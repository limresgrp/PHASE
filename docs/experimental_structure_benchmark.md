# Experimental Structure Benchmark

## Purpose and Scope

Compare supplied experimental/static structures using several Potts-model contrasts and user-authored labels. The output is quantitative evidence only. PHASE does not infer biological meaning, efficacy, activation, allostery, causality or mechanism from the results. Labels are annotations supplied by the scientist, not predictions.

A-like and B-like mean only lower and higher ordering along `E_A - E_B` within the selected structure set. They do not assign biological states. Rank coordinates are not probabilities. Potts scores are in the fitted model's energy units, not experimentally measured free energies.

This is an independent analysis. Existing model-pair, Sampling Explorer, hotspot and other results are not modified.

## Preparation

1. Create a PDB-only state for each experimental structure on the System page (or console), then build descriptors. Missing residues are allowed.
2. Assign these states to the same cluster categorical basis. Previously added states are read through their individual assigned MD sample files, not through the original cluster frame table.
3. Use full or combined Potts models from that cluster. Raw delta correction patches are not full endpoints and are rejected.
4. Select at least two static states and one or more ordered contrasts `(model A, model B)`. At least two contrasts are needed for a 2D axis map; three compatible contrasts form a model triangle for the algebra check.

A selected assignment must contain exactly one valid frame. The benchmark never silently selects a trajectory's first frame or averages a trajectory and calls it an experimental structure. Create a separate PDB-only state if necessary. Ordinary unassigned labels (-1) are errors, not masked missing residues. Missing coordinates use the framework's -2 marker.

By default, missing assignments can be generated automatically for descriptor-ready, single-frame states using the existing assignment workflow with one core. This may take time for large clustering models. It does not rebuild missing descriptors or change a trajectory state into a static state. Disable this option to require all assignments in advance. If multiple assignments exist for one state, explicitly select the intended one.

## Web Workflow

Open **Structure benchmark** in the system top navigation. The left sidebar contains the cluster and all existing benchmark results. Click a result to open its dedicated result URL.

Use **Label schemes** to create/load/edit/delete system-level annotations. They are saved under:

```text
systems/<system_id>/analysis_annotations/label_schemes/<scheme_id>.json
```

They do not rename states or modify scientific state metadata. A scheme contains a name, categorical/ordinal kind, categories, and a map of state IDs to labels. Labels can be arbitrary strings. Categorical lists have no implied order. For ordinal schemes, explicitly enter categories from low to high; PHASE never guesses order from their names. Unlabelled structures remain visible but are excluded from label-based statistics.

Click **New benchmark** to choose static states, the scheme, contrasts, AUC groups, fingerprint/distance methods, hotspot mode and numerical options. All run parameters are in the overlay, not the sidebar. Jobs run asynchronously and report progress; a single RQ worker is sufficient, with no nested queue-job waiting.

The saved label scheme is snapshotted when the web job is submitted, including its revision. Editing/deleting the scheme while a job waits does not change the submitted labels. Saved results retain the snapshot permanently unless the result itself is deleted.

Result tabs:

- **Rank maps & scores**: choose any pair of contrasts or generate all unordered axis-pair maps automatically (six maps per page to limit browser load). Coordinates are raw-support or common-support ranks, default common support. Color and marker shape encode only the supplied categories. Tooltips contain identities, labels, raw/common scores, both ranks and node/edge coverage. A full score table is below.
- **Label statistics**: ordinal Spearman and optional two-category AUC, separately for raw and common support.
- **Axis relationships**: complete rank-Spearman and common-score Pearson matrices, plus exact raw-score algebra checks for complete triangles.
- **Local fingerprints / PCoA**: embedding, positive-eigenvalue fractions, negative-eigenvalue mass and both full distance matrices.
- **Global / local differences**: structure-pair scatter and a configurable threshold-selected table. Changing display thresholds does not rewrite the immutable analysis.
- **Hotspots**: node and edge rankings and cross-axis top-k comparisons. The page shows the top 100 per kind/axis; downloads contain all rankings.
- **Coverage QC**: coverage associations, common term counts and raw/common rank agreement.

Use the Info button for provenance and parameters, Delete to remove a result, and **Download CSV / JSON / NPZ bundle** for the complete export.

## Console and Python CLI

Open a cluster in `./scripts/phase_console.sh` and choose **Experimental Structure Benchmark**. The interactive workflow asks for static states, a saved or new label scheme, several model contrasts, optional AUC categories and hotspot references, fingerprint settings and thresholds. New label schemes are persisted at system level just as on the webserver.

The standalone interactive wrapper is:

```bash
./scripts/experimental_structure_benchmark.sh
```

For unattended reproducible runs:

```bash
python -m phase.scripts.experimental_structure_benchmark \
  --root "$PHASE_DATA_ROOT" --project-id PROJECT --system-id SYSTEM \
  --cluster-id CLUSTER --config benchmark.json --progress
```

Example configuration (replace all IDs):

```json
{
  "name": "Supplied experimental structures",
  "state_ids": ["structure_1", "structure_2", "structure_3"],
  "contrasts": [
    {"model_a_id": "model_A", "model_b_id": "model_B"},
    {"model_a_id": "model_A", "model_b_id": "model_C"},
    {"model_a_id": "model_B", "model_b_id": "model_C"}
  ],
  "label_scheme": {
    "name": "Scientist-supplied ordering",
    "kind": "ordinal",
    "categories": ["group_low", "group_middle", "group_high"],
    "labels": {
      "structure_1": "group_low",
      "structure_2": "group_middle",
      "structure_3": "group_high"
    }
  },
  "positive_label": "group_high",
  "negative_label": "group_low",
  "best_oriented": false,
  "fingerprint_mode": "standardized",
  "distance": "euclidean",
  "hotspots_mode": "parameter",
  "top_k": 20,
  "global_threshold": 0.15,
  "local_threshold": null,
  "seed": 0
}
```

Supply exactly one of `label_scheme` (explicit snapshot) or `scheme_id` (saved system annotation). `state_samples` optionally maps state IDs to explicit assigned sample IDs. For ensemble hotspots, each contrast can include `reference_a` and `reference_b`, lists of assigned MD sample IDs. If omitted, PHASE attempts suggestions from the endpoint model training states; empty or overlapping references fail explicitly. Prefer independent held-out reference MDs where possible.

Label management without running an analysis:

```bash
python -m phase.scripts.experimental_structure_benchmark --root "$PHASE_DATA_ROOT" \
  --project-id PROJECT --system-id SYSTEM --list-label-schemes
python -m phase.scripts.experimental_structure_benchmark --root "$PHASE_DATA_ROOT" \
  --project-id PROJECT --system-id SYSTEM --save-label-scheme labels.json
python -m phase.scripts.experimental_structure_benchmark --root "$PHASE_DATA_ROOT" \
  --project-id PROJECT --system-id SYSTEM --delete-label-scheme SCHEME_UUID
```

A scheme JSON can include `scheme_id` to update an existing scheme; revisions increment. Direct `label_scheme` configurations are snapshotted but not automatically registered as reusable annotations. Use the save command or interactive workflow to persist them.

## Mathematical Definitions

### Scores, Masks and Common Support

Both full endpoint models are put into the same zero-sum gauge, following the existing delta-energy and hotspot convention. A contrast is:

```text
d(x) = E_A(x) - E_B(x)
     = sum_r delta_h_r(x_r) + sum_(r,s) delta_J_rs(x_r,x_s)
```

Use the union of edges from all selected models. An interaction absent from a model is mathematically zero; this is different from an unobserved structural term. Per-structure missing node contributions and edges touching missing nodes remain NaN in the NPZ, never observed zero values. Raw totals sum only the axis's observed terms and always report counts and fractions. An entirely unsupported score is unavailable, not zero.

Common support is fixed globally across all selected structures: a node is supported only if observed in every structure; an edge requires both endpoint nodes in every structure. Each axis sums its fields and its union-of-endpoint-model interactions on this support. Absent couplings on other axes add zero, so the algebra check uses precisely the same structural support. Comparisons across different coverage must be interpreted alongside the common-support variant; raw masked totals are not coverage-corrected estimates of a full structure's energy.

If no common terms remain, raw results and coverage are still exported. Common scores, distances and PCoA are unavailable. Undefined values are not imputed.

### Percentiles and Label Statistics

For the finite scores of each contrast and support variant:

```text
percentile = (average_tie_rank - 1) / (n_finite - 1)
```

A singleton or an all-equal finite set has percentile 0.5. Midrank ties mean the extreme observed groups need not lie exactly at 0 or 1. Percentiles depend on the selected structure set, not on an external calibrated distribution. Adding structures can change them. Do not use rank coordinates in additive energy identities.

Spearman compares score ranks to the explicit ordinal category indices. Report rho, p and n. It measures monotonic ordering, not linearity or effect magnitude. Constant variables produce unavailable coefficients; p is unavailable below three supported observations. The p-values use SciPy's standard asymptotic Spearman implementation, not an exact small-sample/permutation test. Optional BH q-values use two prespecified families: all raw-axis tests, and all common-support-axis tests. Missing p-values are excluded from the correction.

AUC uses only the supplied positive and negative categories; all other labels and unsupported scores are excluded. Higher `E_A - E_B` is always the positive-score direction. Ties receive half credit. AUC below 0.5 is retained, not silently flipped. The optional best-oriented value is `max(AUC, 1-AUC)` and is reported separately. Both group sizes and n are reported. Missing either group makes AUC unavailable. This measures separation of the two supplied groups by the PHASE ordering; it does not validate the labels or imply biological efficacy.

### Axis Redundancy and Algebra

Export the complete pairwise Spearman matrix of common-support ranks and the Pearson matrix of common-support raw scores, including supported sample counts. Related axes are not independent evidence.

When three selected contrasts connect the same three models, account for direction and test:

```text
(E_A-E_B) + (E_B-E_C) - (E_A-E_C) = 0
```

Report maximum and mean absolute residual, n, and a scale-aware floating-point tolerance. Only raw, identical-support scores participate, never ranks. A direction-reversed supplied contrast is handled with the appropriate sign. No residual can be tested when support is absent.

### Hotspots

Reuse `phase.potts.hotspots.contrast_models`, `block_statistics` and `separation_statistics`:

- Parameter mode: alphabet-normalized node RMS and edge Frobenius RMS of zero-sum parameter contrasts, over full model terms (not restricted to the selected structures' common support).
- Ensemble mode: `mean_B(delta term) - mean_A(delta term)` with equal reference-trajectory weighting and the existing within-trajectory block bootstrap. Missing observations are excluded. CIs need the existing minimum block support; a single PDB has a mean but not a reliable narrow CI.

Rank nodes and edges by absolute score. Compare top-k identity overlap, Jaccard, and absolute-score Spearman on shared supported identities. All ties at the kth cutoff are included; actual set size can exceed k. Zero/tied scores do not establish a distinctive ranking. Ensemble hotspot references have their own observation support, separate from the common support of the static fingerprints. Full CI, count and block arrays are retained in the NPZ in ensemble mode. Neither mode predicts causal or mutational consequences.

### Local Fingerprints, Distances and PCoA

For each structure, concatenate node and edge contributions for each axis on the fixed common support. Features record axis and canonical term indices. The default is term-wise z-scoring across the selected structures before concatenating axes (equivalently, independently standardizing the concatenated columns), with population standard deviation (`ddof=0`). Constant/numerically constant columns become zero; means and standard deviations are saved. A raw fingerprint option is also supported. Equal fingerprints do not establish identical atomic coordinates. Axes with more supported terms contribute more features; axis dimensionality is not separately normalized. Redundant axes contribute redundant features; use the redundancy panel to decide which contrasts to include.

Export:

- Correlation distance `1 - Pearson(fingerprint_i, fingerprint_j)`, using the selected raw/standardized fingerprint mode. Constant feature vectors have undefined off-diagonal correlation distances, not arbitrary zero distances.
- Euclidean distance on standardized fingerprints, regardless of the raw-fingerprint display/analysis option.

Classical PCoA double-centers squared distances and diagonalizes the resulting Gram matrix. Save all eigenvalues, up to three positive-axis coordinates, their fractions of **total positive eigenvalue mass**, and the negative eigenvalue mass fraction. Negative eigenvalues can occur with non-Euclidean distances and are reported, not hidden. Axes have arbitrary orientation and no intrinsic physical meaning. A zero-distance embedding is explicitly marked degenerate; a matrix with undefined entries cannot be embedded and is not imputed.

### Global/Local Differences and Coverage QC

For each structure pair, global difference is the RMS difference of common-support rank coordinates across the selected contrasts. Local distance is the selected fingerprint metric. Plot both and list pairs satisfying user-selected `global difference <= threshold` and `local distance >= threshold`. These are numerical filters, not biological classifications.

Coverage QC reports raw-score Spearman associations with node, edge and combined term coverage; raw/common rank agreement; per-structure counts/fractions; and common node/edge counts. A methodological warning is flagged at absolute raw-score/combined-coverage rho >= 0.7, even though this is not a significance test. Inspect sample size and unsupported values. Experimental structures may be related or non-independent, so all reported p/q/AUC values need careful study design and independent validation.

## Parameters

- `fingerprint_mode`: `standardized` (default) or `raw`.
- `distance`: `euclidean` (default) or `correlation`, selecting the PCoA and local-discordance distance.
- `hotspots_mode`: `parameter` (default) or `ensemble`.
- `top_k`: hotspot set-comparison cutoff, default 20, including cutoff ties.
- `block_length`: contiguous assigned frames per hotspot bootstrap block, default 50. Choose based on trajectory correlation and test sensitivity; it is not a physical time interval by itself.
- `bootstrap`: hotspot replicates, default 200; 0 disables intervals, otherwise >=50 (web maximum 5000).
- `seed`: reproducible hotspot bootstrap seed, default 0; axis i uses seed+i. Axis order is part of the configuration.
- `positive_label`, `negative_label`: both absent to disable AUC, or two distinct declared categories. `best_oriented` defaults false.
- `bh`: report the two-family Spearman BH correction, default true; applies only to ordinal labels.
- `global_threshold`: RMS rank threshold in [0,1], default 0.15.
- `local_threshold`: nonnegative local distance threshold; default null uses the 75th percentile of supported off-diagonal distances. Both thresholds can also be changed at visualization time without modifying stored results.
- `assign_missing`: default true, allowing single-core assignment of descriptor-ready single-frame states without a saved assignment.

## Artifacts and Reproducibility

Each run creates a new, immutable UUID directory; rerunning never overwrites another result:

```text
clusters/<cluster_id>/analyses/experimental_structure_benchmark/<analysis_id>/
  analysis_metadata.json
  analysis.npz
  summary.json
  statistics.json
  scores.csv
  hotspots.csv
  structure_pairs.csv
  euclidean_distance.csv
  correlation_distance.csv
  rank_spearman_matrix.csv
  common_score_pearson_matrix.csv
  labels.json
  supports.json
  provenance.json
```

Metadata records model SHA256 hashes, assigned sample SHA256 hashes, source PDB hashes where available, state and model IDs/names, the annotation snapshot/revision, support definitions, parameters, seed, PHASE/NumPy/SciPy versions, source commit/dirty flag and analyzer source hash. A clean committed source tree is recommended for reproducibility. Model files are loaded from byte snapshots; sample mutation during reading is rejected.

The NPZ contains assignments, original NaN-masked per-axis term contributions, raw/common scores/ranks, supports, raw/standardized fingerprints, feature mappings and normalization parameters, both distance matrices, correlation matrices, PCoA values, and full hotspot scores/uncertainty arrays. JSON uses null for unavailable values, not numeric zero. The web result endpoint serves compact summaries rather than the full term/fingerprint tensors.

Array conventions (`S` structures, `A` contrasts, `N` residues, `E` global edges, `F` fingerprint features):

- `raw_scores`, `common_scores`, `raw_ranks`, `common_ranks`: `(S, A)`.
- `term_contributions`: `(A, S, N+E)`; nodes first, then edges in `edges: (E,2)` order. Missing observations are NaN, and absent model couplings are zero only where structurally observed.
- `node_observed: (S,N)`, `edge_observed: (S,E)`, `common_support: (N+E,)`, `axis_support: (A,N+E)`.
- `fingerprints_raw`, `fingerprints_standardized`: `(S,F)`; `feature_axis`, `feature_term`, `fingerprint_mean`, `fingerprint_std`: `(F,)`.
- Distance matrices: `(S,S)`; axis correlation matrices: `(A,A)`; PCoA coordinates: `(S,3)` when available, with unused positive dimensions zero-filled and marked through the eigenvalue fractions.
- `hotspot_i_node`, `hotspot_i_edge`: full axis-i score vectors; `hotspot_i_edge_index` maps its local edge order into the global `edges` array. Ensemble `hotspot_i_*` mean, separation, CI, count and block arrays use nodes followed by that axis's local edges, not all global edges.

Atomic publication prevents incomplete artifacts appearing as results. Label-scheme edits/deletion do not change or invalidate saved snapshots. Existing framework cleanup removes a benchmark if a selected state, model, or referenced assignment is deleted; deleting a benchmark removes only its result folder.

## API Overview

System prefix: `/api/v1/projects/{project}/systems/{system}/experimental-benchmarks`.

- `GET /context?cluster_id=...`: models, states, assignments and saved label schemes.
- `GET/POST /label-schemes`, `PUT/DELETE /label-schemes/{scheme_id}`: persistent analysis annotations.
- `POST /analyses`: enqueue a benchmark (configuration plus cluster_id).
- `GET /analyses/{analysis_id}?cluster_id=...`: compact immutable result summary.
- `GET /analyses/{analysis_id}/download?cluster_id=...`: CSV/JSON/NPZ ZIP bundle.

Listing/deleting results uses the existing generic cluster analysis endpoints with analysis type `experimental_structure_benchmark`. No endpoint edits a stored benchmark.
