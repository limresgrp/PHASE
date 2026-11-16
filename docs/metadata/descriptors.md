# Descriptor Metadata (`*_descriptor_metadata.json`)

Location
`projects/<project_id>/systems/<system_id>/descriptors/<state_id>_descriptor_metadata.json`

Purpose
Metadata for descriptor NPZ files produced per macro state.

Fields
- `canonical_reference` (bool): when true, NPZ insertion order and `descriptor_keys` preserve the selected reference order; they must not be lexicographically sorted.
- `missing_residue_keys` (array): canonical residues absent from the structure. Their NPZ columns contain NaN vectors, not zeros or fabricated torsions.
- `missing_descriptor_policy` (string): NaN descriptors predict missing label -2; node fields and couplings touching a missing residue are omitted from evaluated energy. Sample NPZs persist the matching boolean `missing_mask`.
- `descriptor_keys` (array of strings) – residue keys present in the NPZ
- `residue_mapping` (object) – mapping from residue keys to selection strings
- `n_frames` (int)
- `residue_selection` (string) – selection used to build descriptors

Notes
- `system.json` should **not** store global `descriptor_keys`. This file is the source of truth.
