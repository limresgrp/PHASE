from __future__ import annotations

import argparse
import re
import shutil
import uuid
import json
from pathlib import Path

from phase.services.project_store import ProjectStore, SystemMetadata
from phase.workflows.macro_states import add_state as add_macro_state, refresh_system_metadata


def _root_path(value: str) -> Path:
    return Path(value).expanduser().resolve()


def _slugify(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9]+", "_", value.strip().lower())
    cleaned = cleaned.strip("_")
    return cleaned


def _unique_id(base: str, existing: set[str]) -> str:
    candidate = base or str(uuid.uuid4())
    if candidate not in existing:
        return candidate
    idx = 2
    while f"{candidate}-{idx}" in existing:
        idx += 1
    return f"{candidate}-{idx}"


def init_project(
    store: ProjectStore,
    name: str,
    description: str | None,
    *,
    use_slug_ids: bool,
) -> str:
    project_id = None
    if use_slug_ids:
        slug = _slugify(name)
        existing = {p.project_id for p in store.list_projects()}
        project_id = _unique_id(slug, existing)
    project = store.create_project(name=name, description=description, project_id=project_id)
    print(project.project_id)
    return project.project_id


def create_system(
    store: ProjectStore,
    project_id: str,
    name: str | None,
    description: str | None,
    *,
    use_slug_ids: bool,
) -> str:
    system_id = None
    if use_slug_ids:
        slug = _slugify(name or "")
        project_meta = store.get_project(project_id)
        existing = set(project_meta.systems or [])
        system_id = _unique_id(slug, existing)
    system = store.create_system(project_id=project_id, name=name, description=description, system_id=system_id)
    print(system.system_id)
    return system.system_id


def clone_system(
    store: ProjectStore,
    project_id: str,
    source_system_id: str,
    name: str,
    description: str | None,
    *,
    use_slug_ids: bool,
) -> str:
    if not name.strip():
        raise ValueError("A name is required for the cloned system.")
    system_id = None
    if use_slug_ids:
        slug = _slugify(name)
        project_meta = store.get_project(project_id)
        system_id = _unique_id(slug, set(project_meta.systems or []))
    cloned = store.clone_system(
        project_id,
        source_system_id,
        name=name,
        description=description,
        system_id=system_id,
    )
    print(cloned.system_id)
    return cloned.system_id


def add_state(
    store: ProjectStore,
    project_id: str,
    system_id: str,
    state_id: str,
    name: str | None,
    pdb_path: Path,
    traj_path: Path | None,
    residue_selection: str | None,
    copy_traj: bool,
    build_descriptors: bool,
    slice_spec: str | None,
    resid_shift: int,
) -> None:
    try:
        add_macro_state(
            store,
            project_id,
            system_id,
            state_id,
            name,
            pdb_path,
            traj_path,
            residue_selection,
            copy_traj,
            build_descriptors,
            slice_spec,
            resid_shift,
        )
    except Exception as exc:
        print(f"Failed to add state '{state_id}': {exc}")
        raise

    print(state_id)


def delete_state(
    store: ProjectStore,
    project_id: str,
    system_id: str,
    state_id: str,
) -> None:
    system = store.get_system(project_id, system_id)
    state = (system.states or {}).get(state_id)
    if state is None:
        raise FileNotFoundError(f"State '{state_id}' not found in system '{system_id}'.")

    dirs = store.ensure_directories(project_id, system_id)
    system_dir = dirs["system_dir"].resolve()
    states_dir = dirs["states_dir"].resolve()
    states_dir.relative_to(system_dir)
    state_dir = store._state_dir(project_id, system_id, state).absolute()
    if state_dir.parent.resolve() != states_dir or state_dir.name in {"", ".", ".."}:
        raise ValueError(f"Invalid storage folder for state '{state_id}': {state_dir}")

    for field in ("descriptor_file", "descriptor_metadata_file", "trajectory_file", "pdb_file",
                  "metastable_labels_file", "metastable_metadata_file"):
        rel_path = getattr(state, field, None)
        if not rel_path:
            continue
        stored_path = store.resolve_path(project_id, system_id, rel_path).absolute()
        shared = [system.reference_pdb_file] + [other.pdb_file for sid, other in system.states.items() if sid != state_id]
        if any(value and store.resolve_path(project_id, system_id, value).absolute() == stored_path for value in shared):
            continue
        # The folder deletion below owns these paths, including the metadata
        # that hydrates the state list. Do not follow a state-folder symlink.
        if stored_path.is_relative_to(state_dir):
            continue
        try:
            abs_path = stored_path.resolve()
            abs_path.relative_to(system_dir)
        except ValueError:
            # Never delete external/original source files.
            continue
        abs_path.unlink(missing_ok=True)

    if state_dir.is_symlink():
        state_dir.unlink()
    elif state_dir.exists():
        shutil.rmtree(state_dir)

    system.metastable_states = [
        meta for meta in (system.metastable_states or [])
        if str(meta.get("macro_state_id") or "") != str(state_id)
    ]
    if not system.metastable_states:
        system.metastable_locked = False
        system.metastable_model_dir = None
        system.analysis_mode = "macro" if system.macro_locked else None
    system.states.pop(state_id, None)
    refresh_system_metadata(system)
    store.save_system(system)
    print(state_id)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Offline system/project setup for PHASE.")
    ap.add_argument("--root", required=True, help="Offline data root (will create projects/)")
    sub = ap.add_subparsers(dest="cmd", required=True)

    sp = sub.add_parser("init-project")
    sp.add_argument("--name", required=True)
    sp.add_argument("--description")
    sp.add_argument("--use-slug-ids", action="store_true")

    sp = sub.add_parser("create-system")
    sp.add_argument("--project-id", required=True)
    sp.add_argument("--name")
    sp.add_argument("--description")
    sp.add_argument("--use-slug-ids", action="store_true")
    sp.add_argument("--reference")
    sp.add_argument("--chain-id")
    sp.add_argument("--selection", default="protein")

    sp = sub.add_parser("clone-system")
    sp.add_argument("--project-id", required=True)
    sp.add_argument("--source-system-id", required=True)
    sp.add_argument("--name", required=True)
    sp.add_argument("--description")
    sp.add_argument("--use-slug-ids", action="store_true")

    sp = sub.add_parser("add-state")
    sp.add_argument("--project-id", required=True)
    sp.add_argument("--system-id", required=True)
    sp.add_argument("--state-id", required=True)
    sp.add_argument("--name")
    sp.add_argument("--pdb")
    sp.add_argument("--source-state-id")
    sp.add_argument("--use-reference", action="store_true")
    sp.add_argument("--chain-id")
    sp.add_argument("--no-build-descriptors", action="store_true")
    sp.add_argument("--traj")
    sp.add_argument("--residue-selection")
    sp.add_argument("--copy-traj", action="store_true")
    sp.add_argument("--slice-spec")
    sp.add_argument("--resid-shift", type=int, default=0)

    for command in ("reference", "alignment", "alignments", "attach-traj", "add-pdb-ids"):
        sp = sub.add_parser(command)
        sp.add_argument("--project-id", required=True)
        sp.add_argument("--system-id", required=True)
        if command in {"alignment", "attach-traj"}:
            sp.add_argument("--state-id", required=True)
        if command in {"alignment", "alignments", "add-pdb-ids"}:
            sp.add_argument("--json", action="store_true", help="Print full standalone alignment metadata")
        if command == "reference":
            sp.add_argument("--pdb")
            sp.add_argument("--chain-id")
            sp.add_argument("--selection", default="protein")
        if command == "attach-traj":
            sp.add_argument("--traj", required=True)
            sp.add_argument("--slice-spec")
            sp.add_argument("--no-build-descriptors", action="store_true")
        if command == "add-pdb-ids":
            sp.add_argument("--pdb-ids", default="")
            sp.add_argument("--pdb-id-file")
            sp.add_argument("--chain-id")
            sp.add_argument("--build-descriptors", action="store_true")
            sp.add_argument("--preview", action="store_true")

    sp = sub.add_parser("delete-state")
    sp.add_argument("--project-id", required=True)
    sp.add_argument("--system-id", required=True)
    sp.add_argument("--state-id", required=True)

    args = ap.parse_args(argv)
    root = _root_path(args.root) / "projects"
    store = ProjectStore(base_dir=root)

    if args.cmd == "init-project":
        init_project(store, args.name, args.description, use_slug_ids=args.use_slug_ids)
    elif args.cmd == "create-system":
        sid = create_system(store, args.project_id, args.name, args.description, use_slug_ids=args.use_slug_ids)
        if args.reference:
            from phase.workflows.reference_states import set_reference
            try:
                set_reference(store, args.project_id, sid, Path(args.reference).expanduser().resolve(), chain_id=args.chain_id, selection=args.selection)
            except Exception:
                store.delete_system(args.project_id, sid)
                raise
    elif args.cmd == "clone-system":
        clone_system(
            store,
            args.project_id,
            args.source_system_id,
            args.name,
            args.description,
            use_slug_ids=args.use_slug_ids,
        )
    elif args.cmd == "add-state":
        system = store.get_system(args.project_id, args.system_id)
        if system.reference_pdb_file:
            from phase.workflows.reference_states import create_reference_state
            if args.resid_shift:
                raise ValueError("Manual shifts bypass alignment and are allowed only in legacy systems.")
            if not (args.pdb or args.source_state_id or args.use_reference):
                ap.error("Choose --pdb, --source-state-id or --use-reference.")
            state = create_reference_state(store, args.project_id, args.system_id,
                state_id=args.state_id, name=args.name or args.state_id,
                pdb_path=Path(args.pdb).expanduser().resolve() if args.pdb else None,
                source_state_id=args.source_state_id, chain_id=args.chain_id,
                trajectory=Path(args.traj).expanduser().resolve() if args.traj else None,
                copy_traj=args.copy_traj, slice_spec=args.slice_spec, build_descriptors=not args.no_build_descriptors)
            print(state.state_id)
            return 0
        if not args.pdb:
            ap.error("Legacy systems require --pdb.")
        if args.resid_shift:
            print("WARNING: manual residue shifts bypass reference alignment.")
        add_state(
            store,
            args.project_id,
            args.system_id,
            args.state_id,
            args.name,
            Path(args.pdb).expanduser().resolve(),
            Path(args.traj).expanduser().resolve() if args.traj else None,
            args.residue_selection,
            args.copy_traj,
            not args.no_build_descriptors,
            args.slice_spec,
            args.resid_shift,
        )
    elif args.cmd in {"reference", "alignment", "alignments", "attach-traj", "add-pdb-ids"}:
        from phase.workflows.reference_states import set_reference, load_alignment, attach_trajectory, import_pdb_ids
        system = store.get_system(args.project_id, args.system_id)
        if args.cmd == "reference":
            if args.pdb:
                system = set_reference(store, args.project_id, args.system_id, Path(args.pdb).expanduser().resolve(),
                                       chain_id=args.chain_id, selection=args.selection)
            print(json.dumps({"reference": system.reference_pdb_file, "chain": system.reference_chain_id,
                              "sequence_length": len(system.reference_sequence or ""), "selected": system.reference_residue_keys}, indent=2))
        elif args.cmd in {"alignment", "alignments"}:
            states = [system.states[args.state_id]] if args.cmd == "alignment" else system.states.values()
            rows = [{"state": state.name, "state_id": state.state_id, "status": state.alignment_status or "legacy",
                     "alignment": load_alignment(store, system, state)} for state in states]
            if args.json:
                print(json.dumps(rows, indent=2))
            else:
                print("State | Status | Target chain | Identity | Reference coverage | Missing selected | Insertions | Deletions | Quality")
                for row in rows:
                    meta = row["alignment"] or {}
                    percentage = lambda key: f"{100 * meta[key]:.1f}%" if key in meta else "-"
                    print(f"{row['state']} ({row['state_id']}) | {row['status']} | {meta.get('target_chain', '-')} | "
                          f"{percentage('identity')} | {percentage('coverage_reference')} | {meta.get('missing_selected_residue_count', '-')} | "
                          f"{meta.get('insertion_count', '-')} | {meta.get('deletion_count', '-')} | {meta.get('quality', '-')}")
                    if args.cmd == "alignment" and meta:
                        print("Missing selected:", ", ".join(meta.get("missing_selected_residues", [])) or "none")
                        print("Use --json for complete forward/reverse correspondence and provenance.")
        elif args.cmd == "attach-traj":
            attach_trajectory(store, args.project_id, args.system_id, args.state_id, Path(args.traj).expanduser().resolve(),
                              build_descriptors=not args.no_build_descriptors, slice_spec=args.slice_spec)
            print(args.state_id)
        else:
            ids = args.pdb_ids.replace(",", " ").split()
            if args.pdb_id_file:
                ids += Path(args.pdb_id_file).read_text().replace(",", " ").split()
            if not ids:
                ap.error("Provide --pdb-ids or --pdb-id-file.")
            result = import_pdb_ids(store, args.project_id, args.system_id, ids,
                chain_ids={pid: args.chain_id for pid in ids} if args.chain_id is not None else None,
                preview=args.preview, build_descriptors=args.build_descriptors)
            if args.json:
                print(json.dumps(result, indent=2))
            else:
                print("PDB | Status | Identity | Reference coverage | Missing selected | Quality")
                for row in result:
                    meta = row.get("alignment") or {}
                    print(f"{row['pdb_id']} | {row['status']} | {meta.get('identity', '-')} | "
                          f"{meta.get('coverage_reference', '-')} | {meta.get('missing_selected_residue_count', '-')} | "
                          f"{meta.get('quality', '-')} {row.get('error', '')}")
            return int(any(row["status"] == "failed" for row in result))
    elif args.cmd == "delete-state":
        delete_state(
            store,
            args.project_id,
            args.system_id,
            args.state_id,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
