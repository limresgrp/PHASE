#!/usr/bin/env python3
"""Compare trajectory ordering by endpoint C-alpha RMSD.

This is deliberately independent of the PHASE package and project metadata.
Each trajectory argument has the form ``label|trajectory|topology|selection|resid_shift``. The
trajectory may be a PDB/mmCIF (topology omitted or repeated) or an XTC/DCD
with a compatible topology.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from typing import Iterable

import numpy as np

try:
    import MDAnalysis as mda
    from MDAnalysis.analysis.rms import rmsd as mda_rmsd
except ImportError as exc:  # pragma: no cover - depends on the user's venv
    raise SystemExit(f"MDAnalysis is required: {exc}") from exc


def _spec(value: str) -> tuple[str, Path, Path, str, int]:
    parts = str(value).split("|", 4)
    if len(parts) not in {3, 4, 5}:
        raise argparse.ArgumentTypeError("expected label|trajectory|topology|selection|resid_shift")
    label, trajectory, topology = (part.strip() for part in parts[:3])
    selection = parts[3].strip() if len(parts) >= 4 else "protein"
    resid_shift = int(parts[4].strip() or "0") if len(parts) == 5 else 0
    if not label or not trajectory or not topology:
        raise argparse.ArgumentTypeError("label, trajectory, and topology are required")
    return label, Path(trajectory).expanduser(), Path(topology).expanduser(), selection or "protein", resid_shift


def _load(spec: tuple[str, Path, Path, str, int]):
    label, trajectory, topology, selection, _ = spec
    if not trajectory.exists():
        raise FileNotFoundError(f"{label}: trajectory not found: {trajectory}")
    if not topology.exists():
        raise FileNotFoundError(f"{label}: topology not found: {topology}")
    universe = mda.Universe(str(topology), str(trajectory))
    atoms = universe.select_atoms(f"({selection}) and name CA")
    if not len(atoms):
        raise ValueError(f"{label}: selection '({selection}) and name CA' contains no atoms")
    return label, universe, atoms


def _residue_keys(atoms, resid_shift: int = 0) -> list[str]:
    return [f"{str(res.resname).strip()}_{int(res.resid) + int(resid_shift)}" for res in atoms.residues]


def _validate_atom_mapping(reference, target, label: str, reference_shift: int = 0, target_shift: int = 0) -> None:
    if len(reference) != len(target):
        raise ValueError(f"{label}: found {len(target)} C-alpha atoms; endpoint has {len(reference)}")
    ref_keys = _residue_keys(reference, reference_shift)
    target_keys = _residue_keys(target, target_shift)
    if ref_keys != target_keys:
        mismatches = [f"{a}!={b}" for a, b in zip(ref_keys, target_keys) if a != b]
        preview = ", ".join(mismatches[:5])
        raise ValueError(f"{label}: C-alpha residue mapping differs from endpoint ({preview})")


def _endpoint_coordinates(spec):
    label, universe, atoms = _load(spec)
    universe.trajectory[0]
    return label, universe, atoms, np.asarray(atoms.positions, dtype=np.float64).copy(), int(spec[4])


def _safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("._") or "trajectory"


def _compute_one(spec, coords_a, coords_b, reference_keys, stride: int, align: bool):
    label = str(spec[0])
    try:
        _, universe, atoms = _load(spec)
        target_keys = _residue_keys(atoms, int(spec[4]))
        if target_keys != reference_keys:
            mismatches = [f"{a}!={b}" for a, b in zip(reference_keys, target_keys) if a != b]
            detail = ", ".join(mismatches[:5]) or f"{len(target_keys)} selected C-alpha atoms; expected {len(reference_keys)}"
            raise ValueError(f"C-alpha residue mapping differs from endpoint ({detail})")
        rows = []
        for frame in range(0, len(universe.trajectory), stride):
            universe.trajectory[frame]
            coords = np.asarray(atoms.positions, dtype=np.float64)
            rmsd_a = float(mda_rmsd(coords, coords_a, center=True, superposition=align))
            rmsd_b = float(mda_rmsd(coords, coords_b, center=True, superposition=align))
            rows.append({
                "trajectory": label,
                "frame": frame,
                "rmsd_to_a": rmsd_a,
                "rmsd_to_b": rmsd_b,
                "coordinate": rmsd_a - rmsd_b,
            })
        return {"label": label, "rows": rows, "error": None}
    except Exception as exc:
        return {"label": label, "rows": [], "error": str(exc)}


def _compute(specs: Iterable[tuple[str, Path, Path, str, int]], endpoint_a, endpoint_b, align: bool, stride: int, workers: int):
    label_a, universe_a, atoms_a, coords_a, shift_a = _endpoint_coordinates(endpoint_a)
    label_b, universe_b, atoms_b, coords_b, shift_b = _endpoint_coordinates(endpoint_b)
    _validate_atom_mapping(atoms_a, atoms_b, label_b, shift_a, shift_b)
    specs = list(specs)
    reference_keys = _residue_keys(atoms_a, shift_a)
    args = [(spec, coords_a, coords_b, reference_keys, stride, align) for spec in specs]
    if workers > 1 and len(args) > 1:
        with ProcessPoolExecutor(max_workers=min(workers, len(args))) as pool:
            result = pool.map(_compute_one, *zip(*args))
            results = list(result)
    else:
        results = [_compute_one(spec, coords_a, coords_b, reference_keys, stride, align) for spec in specs]
    rows = [row for result in results for row in result["rows"]]
    errors = [{"trajectory": result["label"], "error": result["error"]} for result in results if result["error"]]
    return label_a, label_b, rows, errors


def _plot(rows: list[dict], output: Path, title: str) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("[rmsd-ordering] matplotlib unavailable; skipping PNG plot", file=sys.stderr)
        return
    labels = list(dict.fromkeys(row["trajectory"] for row in rows))
    fig, ax = plt.subplots(figsize=(10, 5.8), constrained_layout=True)
    colors = plt.get_cmap("tab10")
    for index, label in enumerate(labels):
        values = np.asarray([row["coordinate"] for row in rows if row["trajectory"] == label], dtype=float)
        if not len(values):
            continue
        ax.hist(values, bins=min(60, max(10, int(np.sqrt(len(values))))), density=True,
                histtype="step", linewidth=1.8, color=colors(index % 10), label=label)
    ax.axvline(0.0, color="black", linewidth=1, linestyle="--", alpha=0.7)
    ax.set_xlabel("C-alpha RMSD to endpoint A - C-alpha RMSD to endpoint B (angstrom)")
    ax.set_ylabel("Density")
    ax.set_title(title)
    ax.grid(True, linestyle=":", alpha=0.35)
    ax.legend(loc="best", frameon=False)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--endpoint-a", required=True, type=_spec, help="label|trajectory|topology|selection|resid_shift")
    parser.add_argument("--endpoint-b", required=True, type=_spec, help="label|trajectory|topology|selection|resid_shift")
    parser.add_argument("--trajectory", action="append", required=True, type=_spec,
                        help="label|trajectory|topology|selection|resid_shift; repeat for each trajectory to evaluate")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--no-alignment", action="store_true", help="Do not remove translation/rotation before RMSD")
    parser.add_argument("--stride", type=int, default=100, help="Evaluate every Nth frame of selected trajectories (default: 100)")
    parser.add_argument("--workers", type=int, default=1, help="Parallel trajectory workers (default: 1)")
    args = parser.parse_args(argv)

    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    if args.stride < 1:
        raise SystemExit("--stride must be at least 1")
    if args.workers < 1:
        raise SystemExit("--workers must be at least 1")
    label_a, label_b, rows, errors = _compute(
        args.trajectory, args.endpoint_a, args.endpoint_b,
        not args.no_alignment, args.stride, args.workers,
    )
    for error in errors:
        print(f"[rmsd-ordering] skipped {error['trajectory']}: {error['error']}", file=sys.stderr)
    if not rows:
        raise SystemExit("No trajectory frames were processed")

    csv_path = output_dir / "rmsd_ordering.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["trajectory", "frame", "rmsd_to_a", "rmsd_to_b", "coordinate"])
        writer.writeheader()
        writer.writerows(rows)
    np.savez_compressed(
        output_dir / "rmsd_ordering.npz",
        trajectory=np.asarray([row["trajectory"] for row in rows]),
        frame=np.asarray([row["frame"] for row in rows], dtype=np.int64),
        rmsd_to_a=np.asarray([row["rmsd_to_a"] for row in rows], dtype=np.float64),
        rmsd_to_b=np.asarray([row["rmsd_to_b"] for row in rows], dtype=np.float64),
        coordinate=np.asarray([row["coordinate"] for row in rows], dtype=np.float64),
    )
    metadata = {
        "endpoint_a": label_a,
        "endpoint_b": label_b,
        "aligned": not args.no_alignment,
        "stride": args.stride,
        "workers": args.workers,
        "skipped_trajectories": errors,
        "selection": "state residue selection and name CA",
        "coordinate": "rmsd_to_a - rmsd_to_b",
        "n_frames": len(rows),
        "trajectories": list(dict.fromkeys(row["trajectory"] for row in rows)),
    }
    (output_dir / "rmsd_ordering_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    _plot(rows, output_dir / "rmsd_ordering.png", f"RMSD ordering: {label_a} vs {label_b}")
    print(f"[rmsd-ordering] frames={len(rows)} trajectories={len(metadata['trajectories'])}")
    print(f"[rmsd-ordering] output={output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
