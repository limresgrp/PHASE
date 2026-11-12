"""Resolve display names without confusing topology indices with residue IDs."""
import re


def residue_names_for_keys(keys, *, reference_records=(), pdb_path=None, residue_mapping=None):
    keys = [str(key) for key in keys]
    records = {str(row["residue_key"]): str(row["resname"]).strip()
               for row in reference_records if row.get("residue_key") and row.get("resname")}
    result = {key: records[key] for key in keys if key in records}
    if len(result) == len(keys) or pdb_path is None:
        return result
    import MDAnalysis as mda

    try:
        universe = mda.Universe(str(pdb_path))
    except Exception:
        return result
    try:
        for key in keys:
            if key in result:
                continue
            selection = (residue_mapping or {}).get(key)
            if not selection:
                match = re.fullmatch(r"res_(-?\d+[A-Za-z]?)(?:_[A-Za-z]+)?", key)
                if not match:
                    continue
                selection = f"resid {match[1]}"
            try:
                residues = universe.select_atoms(selection).residues
                if len(residues) == 1:
                    result[key] = str(residues[0].resname).strip()
            except Exception:
                continue
    finally:
        universe.trajectory.close()
    return result
