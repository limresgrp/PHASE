# Distance Analysis and Conditional Pair JS

## Purpose

Distance Analysis asks whether two ensembles have different residue-pair
microstate statistics *at the same C-alpha distance*. Global pair JS mixes
all distances; conditional JS compares only frames in the same distance bin.

A pair can have low global JS but high conditional JS if its local-state
relationships change with geometry. Conversely, distance conditioning may
remove a global difference caused mainly by different distance occupancies.
This is a descriptive comparison, not evidence of a causal interaction.

## Getting Started

1. Open **Distance analysis** in the system workspace navigation.
2. Select a cluster. Only states with rows in that cluster's merged frame axis
   are offered. Adding a state or assigning an independent MD sample does not
   silently append it to an old merged frame axis.
3. Click **Build / refresh cache** to compute continuous C-alpha pair distances.
   This runs as an RQ background job and works with one worker.
4. Choose inspection trajectories and an optional metastable subset. Click
   **Refresh inspection** after editing the inspection filters.
5. Click a pair-table row to inspect its distance distribution and joint-state
   heatmaps. Choose thresholds or load a saved distance preset to split frames.
6. Click **New conditional-JS analysis**, select conditions A and B, set the
   thresholds/minimum counts, and submit the background job.
7. Select the saved analysis in the sidebar and open **Conditional-JS results**.

## Cache Inputs

Conditional JS needs a residue-pair set. If the cluster has neither a distance
cache nor saved contact edges, first use **Build / refresh cache** and wait for
completion. The analysis button is disabled until pairs exist. The legacy
single-bin fallback only uses saved contact edges; it does not invent an
all-pairs graph. If a built cache contains no pairs, increase the candidate
cutoff and rebuild.

- **Candidate-edge cutoff**, initially 12 angstrom: the union of pairs within
  this cutoff in any merged frame. Every candidate's distances are stored at
  every merged frame, including values above the cutoff.
- **Periodic distances**: use minimum-image distances if the frame has a valid
  simulation box; otherwise use Cartesian distances.
- **Trajectory overrides**: optional raw trajectory paths visible to the
  backend/worker, as state_id=/path/to/trajectory. Stored state topology,
  residue selection/mapping, residue offset, slice and stride still apply.
- All merged rows are included. Cluster-fitting subsampling is not reused.

The candidate cutoff controls which pairs exist in the cache. It is independent
of the existing Potts graph/contact cutoff and the displayed threshold bins.
Extending thresholds beyond the candidate cutoff requires a cache rebuild.
Cache construction never changes the Potts model or the residue-state labels.

## Inspection Controls

The inspection state selector initially chooses one trajectory. Ctrl/Cmd-select
several states, or use **Select all trajectories**. The metastable selector
intersects this selection. Filtering uses saved frame identities and does not
reread the original trajectory.

**Display contact cutoff** defines the pair table's contact fraction: the
fraction of selected frames with distance strictly below that value. It does
not change the candidate set or a fitted model's edge set.

Enable **Only pairs contacting below the cutoff** to actually filter the pair
list, then click **Refresh inspection**. The minimum fraction applies to the
selected inspection trajectories and metastable subset: 0 means at least one
frame below the cutoff, 0.1 means at least 10%, and 1 means every selected frame.
Pairs are filtered before pagination. Their distance distributions still show
all selected frames, including non-contact frames. Without this checkbox the
cutoff only defines the reported contact fraction. Geometry filtering requires
a built distance cache; it is unavailable for the legacy single-bin fallback.

## Finding Distance-Dependent Differences

1. Run a conditional-JS analysis comparing ensembles A and B with nonempty
   thresholds, for example **6,8,10 angstrom** (the web default). Blank thresholds
   produce only global JS: existing one-bin analyses must be rerun to obtain
   distance-resolved results. The distance cache is reusable.
2. Select the saved analysis and open **Conditional-JS results**. The heatmap
   compares the joint cluster distribution of each residue pair in A versus B,
   separately within each distance bin. High JS means these ensembles have
   different conformational distributions even at comparable pair distances;
   it does not merely mean the distance itself is different.
3. Rank by **JS range across bins**, defined as max(bin JS) - min(bin JS).
   For example, normalized JS of 0.05 at short distance and 0.65 at long distance
   gives a range of 0.60. A high, constant JS across all bins has range zero and
   is not evidence of distance-dependent separation. Global JS is displayed
   separately and is not averaged into the bin scores.
4. Require adequate **nA and nB per bin** (display default 20; increase for more
   conservative screening). Only bins meeting both counts and containing a
   saved finite JS contribute. At least two supported bins are needed for a
   range. Missing bins are n/a, never zero. Counts are observations, not estimates
   of independent samples: MD autocorrelation can reduce effective support.
5. Optionally restrict **Use distances below** to a saved threshold. Only whole
   saved bins with upper boundary at or below that threshold participate.
   A cutoff leaving just one bin cannot measure between-bin variation. This
   results filter is separate from the inspection pair-contact filter.
6. Click a pair to see **JS versus distance**, with nA/nB on hover and a dashed
   global-JS reference line. Gaps identify unsupported bins. Occupancies report
   the fraction of pair-valid frames in each condition falling in that bin.
7. Inspect the residue ranking. Its mean and maximum summarize incident
   passing-pair JS ranges, and partner count indicates how many pairs support
   the score. These are edge-derived scores, not single-residue JS values.
8. Open **3D structure**, choose a reference PDB and limit the displayed edges
   (default 40, maximum 200). Colored C-alpha connectors identify the ranked
   pairs. Residue colors average the incident displayed-edge ranges. Green is
   small change; red saturates at a user-set normalized range (default 0.3).
   Gray residues are not scored by the displayed edges. Hover connectors for
   pair identities and range. Stored state residue mappings and numbering
   offsets are applied to align cluster residues with the reference PDB.

The 3D coordinates are those of the reference PDB, not a distance-conditioned
mean structure. A connector need not have the bin's distance in this particular
PDB. The opaque cartoon loads only the reference structure, not the trajectory.

All these rankings and filters are visualization-only: saved JS remains in
nats, while range thresholds and residue scores use normalized 0-1 units.
Changing the raw/normalized display does not change which pairs pass.

**Statistical caution:** JS range is an effect size, not a p-value or a claim of
allostery/causation. Rare bins, different sampling, autocorrelated frames and
screening many pairs can produce apparent changes. Follow up candidates using
independent replicas or a block-resampling analysis with appropriate multiple
testing control. This visualization does not compute confidence intervals or
statistical significance.

The pair table can be searched and sorted by residue keys, contact occupancy,
minimum, median, Q10, Q90, range, frame count or state count. Pagination keeps
large candidate sets manageable.

**Distance thresholds** are positive ascending cutoffs such as 6,8,10.
They define d<6, 6<=d<8, 8<=d<10, and d>=10. Exact boundary values belong to the
bin on their right. Blank thresholds mean one All distances bin.

Save/load/delete named presets from the sidebar. They are JSON UI setup files
under this cluster, shared by users accessing the same PHASE data root.
A preset whose type/range is incompatible with the current cache reports the
reason instead of truncating its bins.

## Distribution Panels

The distance histogram overlays counts by trajectory on common histogram bins.
Dashed vertical lines show applied thresholds. Occupancy labels count all
selected geometry frames in each bin.

Each joint-state heatmap shows p(x_r,x_s | bin), normalized within its own valid
frame subset. The horizontal axis is residue s's cluster and the vertical axis
is residue r's cluster. Cluster IDs have the same meaning in every panel.
An additional All distances heatmap gives the unconditional joint distribution.
A pair with an invalid assignment at either endpoint is excluded from that
pair's heatmap; invalid labels at unrelated residues do not remove the frame.
Empty subsets are explicitly marked as empty.

## Conditional-JS Inputs and Results

- **Conditions A and B** each select one or more trajectories and an optional
  metastable subset. Grouped trajectories are pooled by frame count.
- **Minimum valid frames** applies separately to both conditions in each bin.
  Default 1 reproduces the empirical distributions; increase it to mask bins
  with weak support. Counts below the limit produce n/a, not zero divergence.
- **Analysis name** is optional. Every submission creates a separately stored
  distance_js analysis, listed in the sidebar with a delete action.
- Overlapping conditions are allowed; the result reports overlapping merged
  rows so self-comparisons can be recognized.

The result heatmap compares p_A(x_r,x_s | bin) against p_B(x_r,x_s | bin).
The table also shows global pair JS, per-bin JS, valid counts nA/nB, and
per-bin occupancy for A/B. Occupancy is the bin's pair-valid frame count divided
by all pair-valid frames in its condition. Occupancy differences and conditional
distribution differences answer different questions; inspect both.

Values are stored in nats, with theoretical maximum ln(2). The display defaults
to normalized JS = raw JS / ln(2), on [0,1]. Turn normalization off to show nats.
The results can be ranked by global JS or any conditional bin and are paginated.

Interpret high JS only alongside the counts. Low occupancy can create unstable
empirical differences. These are point estimates; uncertainty intervals and
automatic information-gain ranking are not part of this implementation.

## Existing Clusters Without Distance Data

Without distance_cache.npz, every pair uses one **All distances** bin. Requested
thresholds collapse to this single bin and pair JS equals unconditional pair JS.
The physical histogram is unavailable, rather than displaying invented zero
distances. Conditional-JS results record whether real distances were available.
A present but stale/corrupt cache raises an explicit error and must be rebuilt.

The new distance_js outputs do not overwrite or reinterpret delta_js results.
Existing global JS pages continue to use their original analysis files.

## Files and Console Access

- Continuous geometry: clusters/<cluster>/distance_cache.npz
- Presets: clusters/<cluster>/ui_setups/<preset>.json, setup_type=distance_bins
- Results: clusters/<cluster>/analyses/distance_js/<analysis>/analysis.npz
- Result metadata: analysis_metadata.json beside that analysis NPZ

The NPZ stores edge_index (2,E), edge_global_js (E), edge_bin_js (B,E),
edge_bin_counts_a/b (B,E), edge_valid_counts_a/b (E), and
edge_bin_occupancy_a/b (B,E), along with residue keys and bin/edge labels.
There is no globally padded K-by-K tensor for heterogeneous residue alphabets.
The exact joint matrices shown in inspection are calculated only for the
requested pair from its aligned labels and geometry.

The phase_console distance-cache command also offers **Compute conditional
pair JS**. Scripted equivalent:

python -m phase.scripts.distance_cache --root /path/to/phase-data \
  --project-id PROJECT --system-id SYSTEM --cluster-id CLUSTER \
  js --states-a active --states-b inactive --thresholds 6,8,10 \
  --min-bin-frames 20 --name "Active vs inactive, distance-conditioned"

Both web and console use the same calculation and output schema.

