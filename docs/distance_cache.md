# Optional C-alpha Distance Cache

This is the first infrastructure slice in section 29 of
`PHASE_distance_conditioned_JS_infrastructure_spec.md`: construct sparse pair
distances, preserve exact merged-frame identity, reload/validate them, and
retrieve filtered pair histograms. The web page, saved presets, and independent
conditional-JS jobs now build on this infrastructure; see
[Distance Analysis help](distance_analysis.md). Automatic information-gain
ranking remains future work. Potts fitting, sampling, and existing JS results
are unaffected.

## Running

Activate your PHASE environment, then use **Build/inspect CA distance cache** in
`./scripts/phase_console.sh`, or invoke the wrapper directly:

```bash
bash scripts/distance_cache.sh
```

The wrapper selects a project, system, and cluster. It offers build/refresh,
status, a pair histogram, or conditional-JS analysis. Construction uses **all rows already stored in
the selected cluster's merged arrays**. Its 12-angstrom candidate cutoff is a
configurable starting value, independent of any fitted model's contact cutoff.
Construction is explicit; normal clustering and old analyses never require it.

Equivalent command-line use:

```bash
python -m phase.scripts.distance_cache \
  --root /path/to/phase-data --project-id project --system-id system \
  --cluster-id cluster build --max-distance-cutoff 12 --progress

python -m phase.scripts.distance_cache \
  --root /path/to/phase-data --project-id project --system-id system \
  --cluster-id cluster status

python -m phase.scripts.distance_cache \
  --root /path/to/phase-data --project-id project --system-id system \
  --cluster-id cluster pair --pair res_10 res_20 \
  --state-ids active --thresholds 6,8,10 --bins 60
```

Pairs accept the exact stored residue key or a **zero-based column index** in
`residue_keys`. `--state-ids` and `--metastable-ids` are comma-separated IDs;
their intersection filters the saved rows. Omitting a filter selects all rows.
There is no trajectory processing when retrieving histograms.

Build options:

- `--max-distance-cutoff`: select candidate pairs that reach this distance or
  less in at least one merged frame. The cutoff must be finite and positive.
- `--no-pbc`: use Cartesian distances. By default, use minimum-image distances
  if a frame has a valid periodic box; use Cartesian distances otherwise.
- `--trajectory STATE_ID=/path/to/raw.xtc`: override a missing stored trajectory;
  repeat for multiple states. Supply the original unsliced trajectory that
  matches the state's topology and saved slice, not an independently trimmed
  trajectory. Overrides are recorded in the cache, without changing the state.
- `--progress`: display progress separately for candidate discovery and distance
  extraction. A build can run from the CLI independently of the web process.

## Persistent Schema

The optional artifact lives next to `cluster.npz`:

```text
clusters/<cluster_id>/
    cluster.npz
    distance_cache.npz
```

The compressed NPZ contains only numeric/string arrays (no pickle):

| Array | Shape | Meaning |
| --- | --- | --- |
| `residue_keys` | `(N,)` | Exact cluster residue-column order |
| `edge_index` | `(2,E)` | Sorted, unique pairs, always `r < s` |
| `edge_residue_keys` | `(E,2)` | Human-readable pair identities |
| `is_contact_edge` | `(E,)` | Membership in the existing cluster contact graph |
| `frame_state_ids` | `(T,)` | Exact merged state IDs |
| `frame_metastable_ids` | `(T,)` | Exact merged metastable IDs |
| `frame_indices` | `(T,)` | Exact merged descriptor-row indices |
| `distances` | `(T,E)` | Continuous C-alpha distances, float32, angstrom |
| `metadata_json` | scalar string | Authoritative versioned distance metadata |

Metadata identifies schema version 1, `distance_type=ca_ca`, units, candidate
cutoff, frame/residue/edge counts, candidate-selection rule, periodic policy,
contact membership counts, source paths/file signatures, residue-to-atom
mappings, source frame counts, state slices/strides, and residue offsets.
Thresholds are never baked into this artifact. There is no separate competing
metadata sidecar.

## Alignment and Extraction

The workflow reads unfiltered merged labels and all three identity arrays from
`cluster.npz`, respecting its `metadata_json.merged.npz_keys` aliases. A merged
frame index is a descriptor row, so the source frame is obtained by applying
the state's `slice_spec`. The saved stride is used as a fallback only if no
slice is stored; a slice step and stride are never multiplied twice.

State-specific `residue_mapping` selections map logical residue keys to their
physical C-alpha atoms. When these are unavailable, numeric residue keys are
mapped using `physical resid = logical resid - resid_shift`. The state's saved
residue selection is also applied. Missing or ambiguous C-alpha atoms fail
explicitly rather than silently changing residue order.

The first pass finds candidate pairs across every merged frame. The second pass
stores **every candidate pair at every merged frame**, including distances above
the cutoff. Duplicate source frames are read once per pass and copied to all
corresponding merged rows. The cluster-fitting `max_cluster_frames` option is
not applied; invalid/halo labels do not remove geometry rows.

The full distance matrix is assembled in a temporary disk-backed array, then
compressed. No full trajectory coordinate tensor is retained. Temporary disk
space must accommodate roughly `4*T*E` bytes plus the compressed output.

Before replacing a cache, the workflow reloads the cluster and verifies exact
residue, state, metastable, and descriptor-frame identity. A failed build leaves
the previous cache intact. Reading an existing cache also performs these checks;
equal frame counts alone are insufficient. States assigned **after** the cluster
was constructed are not silently appended to its original merged frame axis.

## Missing Distances and Binning

`load_distance_cache(cluster_path)` returns an optional-data object. When the
distance file is absent, `distances=None` and its fallback is one **All distances**
bin per pair. `distance_bin_indices(None, thresholds, frame_count=T)` returns T
zeros, even when thresholds were requested. This preserves unconditional pair
statistics without pretending that distances were measured as zero.

`pair_distance_histogram(...)` reports `distance_available=false`, a single bin
containing all filtered frames, and no physical histogram/quantiles in that case.
Existing global JS consumers continue reading their original arrays unchanged.
A present but corrupt, unsupported, or misaligned cache raises an error instead
of silently switching to legacy behavior.

For available distances, thresholds partition the axis as `[min,max)`: exact
values 6 and 8 belong to the bins starting at 6 and 8, respectively. Empty bins
retain zero counts. Nonfinite values are excluded with bin ID -1 by the binning
helper; pair histogram queries reject missing geometry. Preset validation rejects
unsupported distance types, unordered thresholds, and thresholds above the
candidate cutoff, since a larger threshold would need a larger candidate set.

The histogram helper returns counts, bin edges, threshold occupancy, quantiles,
and selected-frame/state counts, never the full distance matrix. The current
NPZ reader decompresses the matrix server-side; a chunk-addressable container
can be introduced later if repeated interactive queries warrant it.

## Verification

`tests/test_distance_cache.py` exercises known toy distances, full frame identity,
state/metastable filtering, slices and strides, residue shifts, duplicate rows,
halo-label preservation, sparse/empty candidate sets, distances above the cutoff,
periodic boxes, PDB-only inputs, atomic refresh failure, bin boundaries, preset
compatibility, and the equality of legacy single-bin and unconditional pair
statistics.

```bash
python -m pytest -q tests/test_distance_cache.py
```

A representative build on A2A-small (300 residues, 19,000 merged frames from
seven states, 12-angstrom candidate cutoff) produced 6,490 candidate pairs out
of 44,850 possible pairs. The compressed file was 431,178,957 bytes (411.2 MiB),
versus 493,240,000 bytes for its uncompressed float32 distance matrix. The build,
including reload/alignment validation, took 149 seconds on the development host.
The artifact was written to a temporary benchmark directory; existing project
results were not replaced. These figures are system-specific, not a performance
guarantee. Pair-histogram retrieval was also verified on all 19,000 rows.
