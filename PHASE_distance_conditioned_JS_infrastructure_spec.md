# PHASE — Distance-Conditioned Pair Statistics: Infrastructure and Visualization Specification

## Purpose

This document is an implementation specification for extending PHASE with **distance-aware pair statistics**.

The immediate goal is **not** to modify the Potts Hamiltonian or the sampler. The first implementation phase should make pairwise Cα distances a persistent, frame-aligned variable and use that variable to analyse pair-state statistics and Jensen–Shannon divergence (JS) conditionally on distance.

A later phase will use the same variable as an explicit additional relational variable in the statistical model and sampling/backmapping pipeline. The data representation created now must therefore be designed for that future use.

The intended implementation is for an LLM/code agent working directly in the PHASE repository. Prefer the existing project/cluster/sample abstractions and file formats over introducing a parallel architecture.

---

## 1. Scientific motivation

### Current PHASE representation

PHASE currently represents each residue by a discrete local conformational state and uses sparse pairwise couplings on a structural contact graph.

Conceptually:

```text
x_r = residue r microstate
E(x) = sum_r h_r(x_r) + sum_(r,s) J_rs(x_r, x_s)
```

The current contact graph is built by asking whether a residue pair ever comes within a Cα cutoff in the selected ensemble.

This creates an important loss of information:

```text
pair r-s is an edge
        |
        +-- d = 4.8 Å
        +-- d = 5.7 Å
        +-- d = 7.1 Å
        +-- d = 8.9 Å
        ...
```

All these frames currently share the same pair coupling/statistical treatment.

### Proposed direction

Treat distance as an **additional relational variable** associated with a residue pair.

For a pair `(r,s)` and frame `t`, store:

```text
d_rs(t) = Cα-Cα distance
```

and derive a categorical distance state from user-defined thresholds:

```text
b_rs(t) = distance bin
```

For example:

```text
bin 0: d < 6 Å
bin 1: 6 Å <= d < 8 Å
bin 2: d >= 8 Å
```

The long-term statistical model should eventually be able to represent something of the form

```text
J_rs(x_r, x_s, b_rs)
```

or an equivalent model in which `b_rs` is an explicit variable that can be sampled and subsequently supplied to the reconstruction/backmapping stage.

**Important:** this document deliberately stops before implementing that distance-conditioned Potts model. The present work establishes the infrastructure required to make the distance variable persistent, inspectable, and usable in conditional JS calculations.

---

# 2. Scope of the implementation

## Phase A — distance infrastructure

Implement:

1. persistent frame-aligned Cα pair distances;
2. persistent metadata describing the distance variable;
3. efficient filtering by pair/state/trajectory/metastable state;
4. distance-bin presets;
5. APIs for retrieving distance distributions and aligned pair/state data;
6. a dedicated distance-analysis page;
7. conditional pair-state distributions and conditional JS divergence;
8. integration of distance information into existing pair visualizations where useful.

## Phase B — information-gain analysis

After Phase A is stable, implement automatic identification of residue pairs for which conditioning on distance adds substantial information.

This phase should rank candidate pairs using statistics such as:

- distance variability;
- occupancy of proposed bins;
- difference between pair-state distributions across distance bins;
- mutual information between distance state and residue-state pair;
- improvement of conditional versus unconditional JS / predictive fit where appropriate.

## Explicitly out of scope for this change

Do **not** yet:

- change `PottsModel`;
- change `fit_potts_pseudolikelihood_torch`;
- change Gibbs/REX/SA sampling;
- alter the QUBO representation;
- make distance bins part of the current Hamiltonian energy;
- infer distances from residue states during sampling;
- implement distance-aware backmapping.

The architecture should nevertheless make those future extensions possible without changing the stored distance representation.

---

# 3. Current repository architecture to preserve

The implementation should build on the existing PHASE workflow.

Relevant concepts/files in the current repository include:

```text
phase/workflows/clustering.py
phase/potts/metrics.py
phase/potts/potts_model.py
phase/potts/analysis_run.py
backend/api/v1/routes/analysis.py
backend/api/v1/routes/clusters.py
backend/tasks.py
frontend/src/pages/DeltaJsEvalPage.jsx
frontend/src/pages/DeltaJsTablePage.jsx
frontend/src/pages/DeltaJs3DPage.jsx
frontend/src/components/layout/SystemWorkspaceNav.jsx
frontend/src/api/projects.js
frontend/src/api/jobs.js
```

The current cluster NPZ already stores merged frame identifiers alongside residue labels, including:

```text
merged__frame_state_ids
merged__frame_metastable_ids
merged__frame_indices
merged__labels
merged__labels_assigned
```

The distance implementation must use the **same merged frame ordering**. Do not create a second implicit frame index.

The current `DescriptorState` metadata already contains the source information needed to access trajectories/topologies, including:

```text
pdb_file
trajectory_file
source_traj
stride
slice_spec
residue_mapping
```

Use those existing fields rather than adding another trajectory registry unless a concrete missing case is found.

---

# 4. Core data-model decision

## 4.1 Distance is a first-class variable

Do not treat distance only as a plotting convenience.

The persisted representation should explicitly identify:

```text
continuous variable:
    d_rs(t)

categorical derived variable:
    b_rs(t)
```

For Phase A, only `d_rs(t)` needs to be persisted permanently. The categorical `b_rs(t)` should normally be computed from a saved threshold preset so that presets can be changed without recomputing trajectories.

This distinction is important:

```text
trajectory data
      |
      +--> continuous distance cache
                  |
                  +--> distance preset A --> bins A
                  +--> distance preset B --> bins B
                  +--> future Potts state variable
```

Do not permanently bake arbitrary user thresholds into the raw distance cache.

---

# 5. Which pairs should be stored?

A dense all-pairs matrix is not acceptable as the default representation.

For a system with `N` residues, all-pair storage has `N(N-1)/2` values per frame. This becomes unnecessarily large even for moderate proteins and long trajectories.

Instead use a sparse **candidate-edge set**.

## Recommended strategy

When building the distance cache, define a configurable maximum analysis cutoff:

```text
max_distance_cutoff
```

For every selected frame, identify residue pairs whose Cα distance is below this maximum cutoff. The union of those pairs across the selected trajectories becomes the stored candidate-edge set.

For each candidate edge, store its distance for every aligned frame in the selected ensemble, including frames where the pair is above the maximum cutoff if those values are needed for the user-selected bins.

Example:

```text
candidate edge: (127, 209)
frame 0: 5.3 Å
frame 1: 5.8 Å
frame 2: 6.4 Å
frame 3: 8.9 Å
frame 4: 11.2 Å
...
```

This supports both:

```text
contact cutoff = 6 Å
```

and later binning such as:

```text
<6
6-8
8-10
>10
```

### Important distinction

Do not confuse:

```text
contact cutoff used to define the PHASE graph
```

with:

```text
maximum distance stored for distance analysis
```

The latter should generally be equal to or larger than the former.

A practical initial default can be a configurable value around 10–12 Å, but the exact default should be exposed rather than hard-coded as a scientific fact.

---

# 6. Recommended on-disk representation

Create a versioned distance artifact associated with the cluster/system.

Recommended layout:

```text
clusters/<cluster_id>/
    cluster.npz
    distance_cache.npz
    distance_metadata.json
```

Alternatively, if the repository convention strongly favors putting metadata inside NPZ, keep one NPZ plus `metadata_json`. Do not duplicate authoritative metadata in multiple places.

## Required conceptual arrays

```text
edge_index           shape (2, E)
edge_residue_keys    shape (E, 2) or reconstructable from residue_keys
frame_state_ids      shape (T,)
frame_metastable_ids shape (T,)
frame_indices        shape (T,)
distances             shape (T, E) or equivalent chunked/sparse representation
```

Where:

- `T` = merged aligned frames used for the cluster/distance analysis;
- `E` = number of candidate residue pairs.

The canonical orientation must be `r < s`.

## Metadata should contain

```json
{
  "schema_version": 1,
  "distance_type": "ca_ca",
  "distance_unit": "angstrom",
  "coordinate_source": "atomistic_trajectory",
  "max_distance_cutoff": 12.0,
  "residue_count": 300,
  "edge_count": 1840,
  "frame_count": 12345,
  "selected_state_ids": [],
  "selected_metastable_ids": [],
  "frame_alignment": {
    "state_ids": "...",
    "metastable_ids": "...",
    "frame_indices": "..."
  }
}
```

The exact field names may be adapted to the repository conventions, but the semantics must remain explicit.

---

# 7. Frame alignment is a critical invariant

The following arrays must have exactly the same first dimension/order:

```text
merged__labels[t, r]
merged__frame_state_ids[t]
merged__frame_metastable_ids[t]
merged__frame_indices[t]
distance_cache.distances[t, e]
```

If frame `t` represents `(state=S, local_frame=1532)`, every artifact must refer to that same frame.

Implement validation helpers that can assert alignment before analysis.

Recommended helper concept:

```python
validate_distance_alignment(
    labels_frame_state_ids,
    labels_frame_indices,
    distance_frame_state_ids,
    distance_frame_indices,
)
```

Validation should fail loudly if state/frame identity differs.

Do not align by array length alone.

---

# 8. How to compute distances

Use the same trajectory/topology loading ecosystem already used by PHASE, preferably MDAnalysis.

The first implementation should support:

```text
CA-CA distance
```

as the canonical and only required distance type.

Keep the API extensible for future definitions such as:

```text
COM-COM
CB-CB
minimum heavy-atom distance
```

but do not implement these unless the existing architecture makes them essentially free.

## Computation

For each selected state/trajectory:

1. load the topology/PDB;
2. load the trajectory using the existing state metadata;
3. apply the existing stride/slice semantics;
4. identify the Cα atom corresponding to each PHASE residue key;
5. compute distances for candidate residue pairs;
6. align the resulting values to the exact frame order used in the cluster;
7. append them to the distance cache.

Use vectorized/chunked distance calculations. Do not construct an `N x N x T` tensor in memory for long trajectories.

Use spatial neighbor search for identifying candidate pairs if necessary.

---

# 9. Interaction with clustering subsampling

The existing cluster pipeline can use `max_cluster_frames` when fitting residue clusters.

Distance analysis should **not** silently inherit the clustering subsample.

The preferred behavior is:

```text
cluster fitting:
    may use max_cluster_frames

stored distance cache:
    uses all selected aligned frames

conditional JS:
    uses all frames unless user explicitly filters/subsamples
```

This avoids making the distance statistics dependent on a clustering-performance option.

If a later performance optimization introduces sampling, it must be explicit in metadata.

---

# 10. Contact/candidate-edge compatibility

The repository already has a concept of contact edges and records `contact_edge_count` in cluster metadata.

The new distance cache should expose enough metadata to answer:

```text
Which stored distance edges are PHASE contact edges?
Which stored pairs are merely additional candidate pairs?
```

Recommended metadata:

```json
{
  "contact_cutoff": 6.0,
  "contact_edges_count": 512,
  "distance_candidate_edges_count": 1840
}
```

Do not overwrite the existing meaning of the PHASE graph edge.

---

# 11. Distance-bin presets

The UI must let the user define and save a reusable threshold preset.

A preset describes a partition of the continuous distance axis.

Example:

```json
{
  "preset_id": "uuid",
  "name": "CA distance 6-8-10",
  "distance_type": "ca_ca",
  "bins": [
    {"label": "< 6 Å", "min": null, "max": 6.0, "include_max": false},
    {"label": "6–8 Å", "min": 6.0, "max": 8.0, "include_max": false},
    {"label": "8–10 Å", "min": 8.0, "max": 10.0, "include_max": false},
    {"label": ">= 10 Å", "min": 10.0, "max": null, "include_max": true}
  ]
}
```

The representation must avoid ambiguous boundary ownership.

## Preset compatibility

A preset can be loaded only when:

```text
preset.distance_type == cache.distance_type
```

and the cache supports the full threshold range.

If a preset requires values above the stored maximum distance cutoff, the UI should report that the preset is incompatible and offer to recompute/extend the cache.

Do not silently truncate bins.

---

# 12. New distance-analysis page

Create a dedicated page rather than overloading the current Delta JS page.

Suggested name:

```text
Distance Analysis
```

Suggested route:

```text
/system/.../sampling/distance
```

or an equivalent route consistent with the current router structure.

Add the page to `SystemWorkspaceNav` near the existing JS/model-analysis tools.

## Page layout

### A. Configuration panel

Controls:

```text
Distance variable:
    Cα-Cα

State / trajectory selection:
    [multi-select]

Metastable selection:
    [all / selected]

Analysis cutoff:
    [ 6.0 Å ]

Distance preset:
    [None / saved preset]

Build / refresh distance cache:
    [button]
```

The page should clearly distinguish:

```text
candidate-edge cutoff
```

from:

```text
current display/bin thresholds
```

### B. Pair table

Show stored candidate/contact pairs with:

```text
Residue A
Residue B
Contact occupancy
Min distance
Median distance
Q10
Q90
Distance range
Frame count
State count
Suggested for conditioning (Phase B)
```

Sort/filter by each numerical field.

### C. Distance distribution plot

When a pair is selected, show:

- histogram/KDE-like histogram of Cα distance;
- optional overlay or small multiples by state/trajectory;
- current preset boundaries as vertical lines;
- counts/fractions in each bin.

Use Plotly, matching the existing frontend plotting stack.

### D. Pair-state distribution

For selected pair `(r,s)`, show the residue-state joint distribution in a grid:

```text
All distances

          x_s
x_r       matrix / heatmap
```

and, when a preset is loaded:

```text
<6 Å
6–8 Å
>8 Å
```

Each panel represents:

```text
p(x_r, x_s | distance bin)
```

Use the same state labels and residue ordering as the existing PHASE pair visualizations.

### E. State/trajectory filtering

The same pair can be inspected as:

```text
all selected data
trajectory A
trajectory B
metastable state M1
metastable state M2
```

Filters must operate on the stored frame metadata; never recompute the trajectory distances for a normal UI filter.

---

# 13. Conditional JS analysis

The first distance-aware JS analysis should compare distributions **within the same distance bin**.

For two conditions/samples A and B:

```text
p_A(x_r, x_s | b)
vs
p_B(x_r, x_s | b)
```

and compute:

```text
JS_b(r,s) = JS(p_A(.|b), p_B(.|b))
```

This should remain separate from the current unconditional edge JS output.

Recommended analysis kind:

```text
distance_js
```

Do not overload existing `delta_js` metadata with incompatible array shapes.

---

# 14. Distance-JS output schema

The result should be versioned and self-describing.

Conceptually:

```text
edge_index                  (E, 2)
edge_labels                 (E,)
distance_bin_labels         (B,)
edge_bin_js                 (B, E)
edge_bin_counts_a           (B, E)
edge_bin_counts_b           (B, E)
edge_bin_occupancy_a        (B, E)
edge_bin_occupancy_b        (B, E)
```

Optional but useful:

```text
pair_joint_a                (B, E, Kmax...)
pair_joint_b                (B, E, Kmax...)
```

Because residue alphabets are heterogeneous, do not force a rectangular global `(K,K)` array if the existing code uses edge-specific matrices. A serialized edge-oriented structure is acceptable.

## Minimum viable result

The API must at least be able to return:

```json
{
  "bin_labels": ["<6", "6-8", ">=8"],
  "edges": [
    {
      "r": 127,
      "s": 209,
      "js": [0.12, 0.31, 0.04],
      "counts_a": [1200, 4300, 300],
      "counts_b": [1100, 4100, 500]
    }
  ]
}
```

Use repository-native analysis metadata/results storage rather than inventing a second job persistence system.

---

# 15. Existing JS page integration

Do not immediately rewrite the existing Delta JS analysis.

Instead:

1. keep unconditional `delta_js` unchanged;
2. add `distance_js` as a separate analysis type;
3. allow the existing Delta JS UI to link to/view the distance-conditioned result where appropriate;
4. later add a distance-preset selector to the existing JS visualization if it can be done without making the page ambiguous.

### Delta JS table

Add an optional column or expandable detail:

```text
Global JS    <6 Å JS    6–8 Å JS    >8 Å JS
```

but only when a compatible `distance_js` analysis exists.

The default view must remain the current global JS.

### Delta JS 3D

The 3D page may eventually support:

```text
JS mode:
    Global
    <6 Å
    6–8 Å
    >8 Å
```

For the first implementation, it is acceptable to add this only after the backend/API schema is stable.

---

# 16. Pair-distribution visualization requirements

The current PHASE visualizations already expose pairwise distributions/JS values. Extend them conceptually as follows:

```text
Pair: RES127 — RES209

              Global
             /       \
         distance-conditioned
          /      |       \
       <6Å     6–8Å     >8Å
```

The viewer should make it obvious that conditional panels are **subsets of frames**, not separate residue-state definitions.

Do not change the residue cluster labels when filtering by distance.

---

# 17. Phase B — automatic information-gain analysis

Once distance-conditioned pair statistics work, implement automatic ranking.

The purpose is not merely to find high-variance distances. The goal is to identify pairs where distance conditioning reveals a different relationship between residue microstates.

## 17.1 Geometry variability

For pair `(r,s)` compute:

```text
Q10(d)
Q50(d)
Q90(d)
IQR(d)
range(d)
contact occupancy
```

This is a screening diagnostic only.

## 17.2 Conditional pair-distribution separation

For bins `b1, b2`, compute:

```text
JS(p(x_r,x_s | b1), p(x_r,x_s | b2))
```

Aggregate over valid bin pairs using occupancy weights.

Call this something explicit such as:

```text
conditional_pair_JS
```

## 17.3 Distance-state information

Define a categorical distance variable `B_rs` from the selected preset.

Compute:

```text
I(B_rs ; (X_r, X_s))
```

where `(X_r, X_s)` is treated as the joint discrete pair state.

This directly measures whether the distance regime carries information about the pair's residue-state configuration.

## 17.4 Conditional-versus-global gain

A useful additional diagnostic is:

```text
gain = conditional_statistic - global_statistic
```

However, define the sign and normalization explicitly before implementation. Do not introduce a metric called "information gain" if it is not dimensionally/interpretable as information.

A safer initial ranking can be multi-criteria:

```text
candidate_score =
    f(distance_variability,
      bin_occupancy,
      conditional_pair_JS,
      mutual_information)
```

with each component reported separately in the UI.

Do not hide the individual measurements behind one opaque score.

---

# 18. Automatic candidate suggestions

The UI should eventually suggest pairs such as:

```text
Potentially useful distance-conditioned edges

1. ASP127 — PHE209
   Distance variability: high
   Bin occupancy: valid
   I(B; Xr,Xs): 0.18 bits
   Conditional JS: 0.21

2. GLY143 — GLN144
   Distance variability: medium
   Bin occupancy: valid
   I(B; Xr,Xs): 0.11 bits
   Conditional JS: 0.14
```

Avoid hard-coded scientific claims such as:

```text
"split if the distance goes from 6 to 8 Å"
```

as the only criterion.

Instead, use those values as an initial UI preset and let the data determine whether the pair is informative.

---

# 19. Distance presets in the UI

Presets should be reusable across analyses on compatible systems.

Recommended controls:

```text
Distance preset
    [default: none]

[+ Add threshold]

< 6 Å
6–8 Å
> 8 Å

[Save preset]
```

For each pair, allow the user to override which bins are active only for analysis/display.

Do not mutate the globally saved preset from a pair-specific override.

Suggested metadata fields:

```text
preset_id
name
distance_type
thresholds
created_at
system_id (optional)
source_analysis_id (optional)
```

---

# 20. API design

Follow the existing FastAPI route structure.

Suggested endpoints (names can be adapted to repository conventions):

```text
POST   /projects/{project}/systems/{system}/clusters/{cluster}/distance-cache
GET    /projects/{project}/systems/{system}/clusters/{cluster}/distance-cache
GET    /projects/{project}/systems/{system}/clusters/{cluster}/distance-pairs
GET    /projects/{project}/systems/{system}/clusters/{cluster}/distance-pairs/{r}/{s}
POST   /projects/{project}/systems/{system}/clusters/{cluster}/distance-presets
GET    /projects/{project}/systems/{system}/clusters/{cluster}/distance-presets
DELETE /projects/{project}/systems/{system}/clusters/{cluster}/distance-presets/{preset}
POST   /projects/{project}/systems/{system}/clusters/{cluster}/analyses/distance-js
```

The exact route naming should follow the current `analysis.py`, `clusters.py`, and frontend API patterns.

The API should expose:

- cache status;
- cache metadata;
- compatible presets;
- pair summaries;
- filtered distance distributions;
- conditional pair-state distributions;
- distance-conditioned JS results.

Avoid returning full `T x E` distance matrices to the browser. The frontend should request aggregated distributions or selected-edge series.

---

# 21. Background jobs

Distance-cache construction may be expensive for long trajectories.

Use the existing RQ job infrastructure.

Recommended stages:

```text
submit job
   ↓
load selected trajectories
   ↓
compute candidate edge set
   ↓
compute/cache aligned distances
   ↓
write metadata
   ↓
job finished
```

Do not make the browser wait on trajectory processing.

Conditional JS over an existing distance cache can be synchronous for small selections or an RQ job for larger analyses, depending on empirical runtime.

The analysis result should use the existing analysis metadata/result conventions.

---

# 22. Performance requirements

The implementation should be designed around long MD trajectories.

Avoid:

```python
all_xyz = np.asarray(all_frames)
all_distances = cdist(all_xyz, all_xyz)
```

for the full trajectory.

Prefer:

```text
trajectory chunks
    ↓
CA coordinates for each chunk
    ↓
spatial neighbor search / pair list
    ↓
selected edge distances
    ↓
append to cache
```

Potential optimizations, if needed:

- spatial cell lists / KDTree;
- chunked vectorized pair differences;
- float32 distance storage;
- memory-mapped arrays;
- compressed NPZ only after the final arrays are assembled.

Do not prematurely optimize the first prototype if the code is already fast for realistic PHASE systems. Measure first.

---

# 23. Validation and tests

Every implementation step must include tests.

## Distance computation unit tests

Use a toy topology with known Cα positions.

Verify:

```text
d(A,B) = expected
```

including unit handling and atom selection.

## Frame alignment test

Construct two trajectories/states with known frame IDs and verify that distance values map to the same `(state_id, frame_index)` as the cluster labels.

## Boundary test

For a preset:

```text
<6
6–8
>=8
```

verify exact boundary assignment at:

```text
5.999
6.000
7.999
8.000
```

## Filtering test

Given known state IDs, verify that trajectory filtering produces the correct subset of distances and labels.

## Conditional pair distribution test

Construct a toy dataset where pair states are deliberately correlated with distance bins.

Verify that:

```text
p(xr,xs | bin 1)
!=
p(xr,xs | bin 2)
```

and that JS is positive.

## Null test

Construct a toy dataset where distance is statistically independent of pair state.

Verify that conditional JS is close to zero within sampling noise.

## Sparse storage test

Verify that only candidate edges are persisted and that edge ordering is canonical.

## Cache compatibility test

Verify that incompatible distance types or threshold ranges are rejected instead of silently producing truncated results.

---

# 24. Acceptance criteria for Phase A

Phase A is complete only when all of the following are true:

### Data

- Cα distances are computed from the actual selected trajectories.
- Distances are aligned exactly to cluster frame identity.
- Distance metadata is versioned.
- Candidate pairs are sparse.
- Raw continuous distance values are preserved.

### UI

- A dedicated Distance Analysis page exists.
- User can select trajectories/states.
- User can inspect pair distance distributions.
- User can define/save/load distance threshold presets.
- Pair-state distributions can be split by distance bin.
- Existing global JS views still work unchanged by default.

### JS

- Conditional pair-state JS can be computed between two selected conditions/samples.
- Results are stored as a separate analysis type.
- Distance-bin occupancy is exposed with JS values.
- Results can be inspected per pair and per bin.

### Tests

- alignment, boundaries, filtering, sparse storage, conditional JS, and null controls all pass.

---

# 25. Acceptance criteria for Phase B

Phase B is complete when the system can automatically identify candidate edges based on distance-conditioned information.

At minimum, the UI should expose:

```text
Distance variability
Bin occupancy
Conditional pair JS
Mutual information I(B ; Xr,Xs)
```

and provide a ranked list of candidate pairs.

The ranking must be accompanied by the underlying metrics so that users can inspect why a pair was suggested.

No candidate should be labeled “important” solely because its distance variance is large.

---

# 26. Future Potts/backmapping architecture

The long-term extension should treat the distance regime as an additional model variable, not merely as a label attached to the edge.

Conceptually:

```text
residue state x_r
residue state x_s
        |
        +-------- distance state b_rs
                         |
                         v
                J_rs(x_r, x_s, b_rs)
```

For actual sampling/backmapping, the model must eventually define how `b_rs` is generated or constrained alongside the residue-state variables.

The continuous distance cache created in Phase A must therefore remain independent of the Potts parameterization.

A later implementation may explore:

```text
J_rs(xr,xs,b)
```

or a baseline-plus-correction form:

```text
J_rs^b = J_rs + ΔJ_rs^b
```

with strong regularization of `ΔJ`.

The latter is especially compatible with the existing delta-Potts infrastructure.

However, do not implement either representation in Phase A.

---

# 27. Important implementation principles

1. **Distance is data, not just a plotting filter.** Persist it as a first-class variable.

2. **Keep continuous and categorical representations separate.** Store distances; derive bins from presets.

3. **Preserve exact frame identity.** State ID + local frame index must remain the alignment key.

4. **Do not change the existing PHASE state clustering.** Distance conditioning should operate on the existing residue-state labels.

5. **Do not modify the current Potts model yet.** Phase A is an analysis/infrastructure extension.

6. **Do not overload `delta_js`.** Introduce a versioned distance-aware analysis type.

7. **Keep global statistics available.** Conditional plots must supplement, not silently replace, the current JS results.

8. **Separate contact cutoff from distance-analysis range.** The latter may extend beyond the former.

9. **Avoid dense all-pairs storage.** Use a sparse candidate-edge set.

10. **Do not equate variance with information gain.** Candidate detection must eventually include state-distance dependence.

11. **Keep future sampling/backmapping in mind.** The distance variable must be representable independently of any specific Potts implementation.

---

# 28. Recommended implementation order

```text
1. Add distance-cache schema and metadata helpers.

2. Implement trajectory → aligned Cα distance extraction.

3. Integrate distance-cache generation into cluster/system workflows.

4. Add cache validation and compatibility checks.

5. Add backend endpoints for cache metadata and selected pair distributions.

6. Add Distance Analysis frontend page.

7. Add threshold/preset persistence.

8. Add conditional pair-state visualization.

9. Add distance-conditioned JS analysis as a separate analysis type.

10. Integrate conditional JS results into existing pair/Delta-JS views.

11. Add automated information-gain metrics.

12. Add candidate ranking/suggestions.

13. Benchmark runtime/storage on a realistic PHASE system.

14. Only after the above passes, design the distance-aware Potts/sampling layer.
```

---

# 29. First implementation task for the coding agent

Before changing the frontend, implement and test the following narrow vertical slice:

```text
selected cluster
    ↓
load exact merged frame identity
    ↓
load original selected trajectories
    ↓
extract Cα coordinates
    ↓
compute sparse candidate residue-pair distances
    ↓
save distance cache + metadata
    ↓
reload cache
    ↓
verify exact frame alignment
    ↓
retrieve distance histogram for one pair
```

Only after this vertical slice is working should the UI and conditional JS work begin.

The coding agent should report:

- which files were changed;
- what the persistent schema is;
- how frame alignment is guaranteed;
- how candidate pairs are selected;
- storage size on a representative test system;
- tests added and their results.

Do not proceed to Potts changes as part of this implementation task.
