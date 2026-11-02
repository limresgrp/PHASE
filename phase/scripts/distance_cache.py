"""Build optional cluster distances or retrieve one filtered pair histogram."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

from phase.services.project_store import ProjectStore
from phase.workflows.distances import build_distance_cache, load_distance_cache, pair_distance_histogram, read_cluster_frames
from phase.workflows.distance_analysis import run_distance_js_analysis


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=os.environ.get("PHASE_DATA_ROOT", "data"))
    parser.add_argument("--project-id", required=True)
    parser.add_argument("--system-id", required=True)
    parser.add_argument("--cluster-id", required=True)
    commands = parser.add_subparsers(dest="command", required=True)
    build = commands.add_parser("build", help="Explicitly build/refresh the optional cache")
    build.add_argument("--max-distance-cutoff", type=float, default=12.0)
    build.add_argument("--no-pbc", action="store_true", help="Disable minimum-image distances")
    build.add_argument("--trajectory", action="append", default=[], metavar="STATE_ID=PATH", help="Override an unavailable stored trajectory; repeat per state")
    build.add_argument("--progress", action="store_true")
    commands.add_parser("status", help="Validate and inspect cache metadata, or report the single-bin fallback")
    pair = commands.add_parser("pair", help="Retrieve a compact pair histogram as JSON")
    pair.add_argument("--pair", nargs=2, required=True, metavar=("R", "S"), help="Residue keys (res_10) or zero-based cluster column indices")
    pair.add_argument("--state-ids", default=None, help="Comma-separated macro state IDs")
    pair.add_argument("--metastable-ids", default=None, help="Comma-separated metastable IDs")
    pair.add_argument("--thresholds", default="", help="Strictly increasing cutoffs, e.g. 6,8,10")
    pair.add_argument("--bins", type=int, default=60)
    js = commands.add_parser("js", help="Persist independent distance-conditioned pair JS")
    js.add_argument("--states-a", required=True, help="Comma-separated states in condition A")
    js.add_argument("--states-b", required=True, help="Comma-separated states in condition B")
    js.add_argument("--metastable-a", default="", help="Optional comma-separated metastable IDs in A")
    js.add_argument("--metastable-b", default="", help="Optional comma-separated metastable IDs in B")
    js.add_argument("--thresholds", default="")
    js.add_argument("--min-bin-frames", type=int, default=1)
    js.add_argument("--name", default=None)
    args = parser.parse_args(argv)
    store = ProjectStore(base_dir=Path(args.root).expanduser().resolve() / "projects")
    entry = store.get_cluster_entry(args.project_id, args.system_id, args.cluster_id)
    cluster_path = store.resolve_path(args.project_id, args.system_id, entry["path"])

    if args.command == "build":
        overrides = {}
        for item in args.trajectory:
            state_id, sep, path = item.partition("=")
            if not sep or not state_id.strip() or not path.strip():
                parser.error("--trajectory must be STATE_ID=PATH")
            overrides[state_id.strip()] = str(Path(path).expanduser().resolve())
        bar = None
        def progress(message, current, total):
            nonlocal bar
            from tqdm import tqdm
            if bar is None or bar.desc != message:
                if bar is not None:
                    bar.close()
                bar = tqdm(total=total, desc=message, unit="frame")
            bar.update(current - bar.n)
        try:
            output = build_distance_cache(
                args.project_id, args.system_id, args.cluster_id, store=store,
                max_distance_cutoff=args.max_distance_cutoff, periodic=not args.no_pbc,
                trajectory_overrides=overrides, progress_callback=progress if args.progress else None,
            )
        finally:
            if bar is not None:
                bar.close()
    elif args.command == "js":
        def csv(value):
            return [item.strip() for item in value.split(",") if item.strip()]
        a, b = {"state_ids": csv(args.states_a)}, {"state_ids": csv(args.states_b)}
        if args.metastable_a:
            a["metastable_ids"] = csv(args.metastable_a)
        if args.metastable_b:
            b["metastable_ids"] = csv(args.metastable_b)
        output = run_distance_js_analysis(
            args.project_id, args.system_id, args.cluster_id, store=store,
            condition_a=a, condition_b=b, thresholds=[float(value) for value in csv(args.thresholds)],
            min_bin_frames=args.min_bin_frames, name=args.name,
        )
    elif args.command == "status":
        output = load_distance_cache(cluster_path).metadata
    else:
        keys = read_cluster_frames(cluster_path).residue_keys.tolist()
        indices = [keys.index(value) if value in keys else int(value) for value in args.pair]
        def csv(value):
            return [item.strip() for item in value.split(",") if item.strip()] if value is not None else None
        output = pair_distance_histogram(
            cluster_path, *indices, state_ids=csv(args.state_ids), metastable_ids=csv(args.metastable_ids),
            thresholds=[float(item) for item in csv(args.thresholds)], bins=args.bins,
        )
    print(json.dumps(output, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
