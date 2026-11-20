# Model-pair Discriminative Hotspots

## Purpose

Prioritize residues and interactions that distinguish two endpoint Potts Hamiltonians. This complements commitment/frustration: parameter differences describe **where the models differ**, while reference-ensemble contrasts test **which differences discriminate the endpoints actually sampled**.

These are candidate-prioritization scores, not causal allostery measurements or mutation free energies. A large score does not establish that mutating the residue will change function. Validate candidates using independent MD and experiments; include low-score controls and specify hypotheses before inspecting experimental results.

## Run From the Webserver

Open **Model-pair analysis**, then **Discriminative hotspots**. Select the cluster and click **New analysis**. Hyperparameters live in this dialog; the sidebar lists existing analyses and visualization filters. Existing commitment, frustration, and delta-energy results are unchanged.

Select two different **full or combined** models from the same cluster. A raw delta patch is a correction, not a full endpoint, and is rejected. Models must use the same categorical alphabets and residue ordering. Available residue-key/cluster metadata are checked; old model files without keys necessarily assume the selected cluster's categorical basis. Do not compare independently reclustered models even if their alphabet sizes happen to match.

### Advanced: All Model Pairs

In **New analysis → Advanced: compare multiple models**, enable **Run all unique
model pairs** and select two or more models. Selecting A, B and C creates three
independent analyses: A–B, A–C and B–C. For N models there are N(N−1)/2 jobs;
the dialog previews their orientation and count. The earlier selected model is
endpoint A, so delta energy always means E(A) − E(B).

In ensemble mode, check/edit the reference MDs **for each model**. Each pair
uses those two models' respective reference groups; do not reuse one fixed A/B
reference group for every pair. Overlapping/invalid reference groups still fail
the normal validation. Additional ensembles and bootstrap settings are shared.
Jobs are queued independently and work with one worker; a failed pair does not
cancel the others. Batch progress tracks all jobs, and each successful result
is listed separately in Available analyses. Failed submissions are reported.

The interactive console offers the same advanced multi-model prompt before
asking for the remaining parameters, and asks for references once per model.
For a noninteractive run, use `--model-ids a,b,c` instead of `--model-a-id` and
`--model-b-id`, with optional repeated `--model-reference a:md_a1,md_a2`,
`--model-reference b:md_b` and `--model-reference c:md_c` arguments. Omitted
per-model reference selections use the existing model-state inference, which
can be ambiguous; explicit references are preferable. Local pairs run in order,
retain successful results on failure and exit nonzero if any pair failed.

Two analysis modes:

- **Parameter contrast only** needs no trajectories. Use it as a first-pass screen.
- **Ensemble-supported discrimination** also requires one or more assigned MD reference samples for each endpoint. Each reference MD sample must represent one source state (select separate samples for multiple trajectories). References are suggested from model training states where available; inspect and correct these choices. Prefer independent held-out replicas. Samples added after clustering are supported through their independent MD assignment NPZs; they need not occur in the original cluster frame table.

Select optional additional MD or sampled ensembles if you want to color their endpoint preferences later. A/B references are included automatically.

Jobs run asynchronously, without nested queue jobs or waiting for another worker. One web worker is sufficient. Progress reports trajectory evaluation and bootstrap stages. No uploaded trajectory coordinates are required: this analysis consumes existing assigned categorical samples.

## Run From the Console

Open a cluster in `./scripts/phase_console.sh`, then choose **Model-pair discriminative hotspots**. Select models, mode, and endpoint reference MD samples. The interactive wrapper is also independently runnable:

```bash
./scripts/potts_model_pair_hotspots.sh
```

The Python CLI additionally accepts optional preference ensembles:

```bash
python -m phase.scripts.potts_model_pair_hotspots \
  --root "$PHASE_DATA_ROOT" --project-id PROJECT --system-id SYSTEM \
  --cluster-id CLUSTER --model-a-id MODEL_A --model-b-id MODEL_B \
  --reference-a MD_SAMPLE_A --reference-b MD_SAMPLE_B \
  --sample-ids EXTRA_SAMPLE_1,EXTRA_SAMPLE_2 \
  --block-length 50 --bootstrap 200 --seed 0 --progress
```

Reference options accept comma-separated sample IDs. Omitting either option attempts inference from the corresponding model's training state IDs. Ambiguous overlapping references or empty groups fail with an actionable error; explicit choices are preferable. The console wrapper lists MD samples only for reference selection.

Use `--mode parameter` for parameter-only analysis.

## Scores and Signs

Both models are transformed into a common **zero-sum gauge**, which makes the decomposition reproducible. Let `d = E_A - E_B`. Use the **union** of their edge sets; a missing coupling is zero. An edge present in only one model is not discarded.

Parameter scores:

```text
node RMS = sqrt(mean((h_A - h_B)^2))
edge RMS = sqrt(mean((J_A - J_B)^2))
```

These equal L2/sqrt(K) and Frobenius/sqrt(K_r*K_s). Alphabet-size normalization avoids a score increasing simply because more categorical parameters exist. Scores remain in model energy units, not on a 0-1 scale, and may depend on fitting regularization. Unsigned strength has no A/B sign.

For each node or edge term, ensemble separation is:

```text
S = mean_B(d) - mean_A(d)
```

First average observed frames within each reference trajectory, then average trajectories with **equal weight**, independently for A and B. A long trajectory does not overwhelm a short replica. Positive S supports the expected endpoint distinction; negative S is counter-discriminative. Rank by absolute S, then assess sign, coverage and uncertainty. Do not rank by significance alone.

For preference coloring of a selected ensemble:

```text
preference = 0.5 * (mean_A(d) + mean_B(d)) - mean_selected(d)
```

Positive (red) is lower relative delta energy, hence more A-like; negative (blue) is more B-like. This is relative to the balanced reference midpoint, not an absolute thermodynamic classification. Counter-discriminative terms are **not silently sign-flipped**.

## Inputs and Uncertainty

- **Block length**, default 50: number of consecutive assigned frames per nonoverlapping bootstrap block. This is not 50 ps unless your frame interval is 1 ps. Choose a length exceeding relevant correlation times and repeat with different lengths. Gaps in source frame indices and changes of source state break blocks.
- **Bootstrap replicates**, default 200: independent block resampling within each reference trajectory; 95% percentile intervals. Set 0 to disable intervals; otherwise use at least 50 (web maximum 5000).
- **Seed**, default 0: reproducible bootstrap draws for fixed inputs and sample order.
- **References**: nonempty, disjoint assigned MD groups. Overlapping source-state/frame identities, including duplicate reference samples, are rejected. Provenance cannot detect physical dependence between differently named states; verify this yourself.
- **Labels**: assigned labels only. Invalid frames and ordinary unassigned (-1) frames are excluded and counted. Missing coordinates (-2) remain masked: a missing node or an edge touching it is **not an observation**, rather than an observed zero. Unsupported scores and intervals remain unavailable.

Intervals require at least four observed blocks per endpoint and at least two in every contributing trajectory. A single PDB can have a mean score but not a misleading narrow interval. Counts and block support are saved per term. Bootstrap uncertainty is conditional on the **fixed fitted models and fixed reference trajectories**; it does not quantify model-fit uncertainty or replace independent replicas. Intervals are not multiple-testing-adjusted p-values.

## Visualization

The ranked tables and **3D network** share the same filters:

- Rank on parameter RMS or ensemble separation; select top N, minimum absolute effect, and/or the upper portion of the ranking using a percentile cutoff. Node and edge rankings are independent.
- Require an individual-term confidence interval excluding zero, or restrict the separation sign.
- Include or exclude a saved **Structures & selections** residue selection. Excluding the pocket also excludes edges touching it. Including a selection includes its incident edges.
- Choose **unsigned importance** (pale to orange) or an ensemble's **endpoint preference** (red A-like, blue B-like). The displayed saturation value is shared between visible node and edge colors; changing filters can change this scale.
- Optional **edge-weighted residue color**, off by default, uses `(1-alpha)*node + alpha*incident_edges`, with alpha default 0.75. Incident edges are averaged by default; sum is available but favors high-degree nodes. Importance combines magnitudes; preference combines signed values and can cancel. At alpha=1, residue color comes entirely from incident edges. All incident edges contribute before display filtering; isolated residues retain the node value.

Mixing is visualization-only: the NPZ always stores separate node and edge scores. Ranking, effect columns, coverage and confidence filtering remain based on **individual unmixed terms**, not uncertainty of the displayed blend. The table's Displayed column makes this distinction explicit.

Select a PDB in the 3D network. It loads one opaque structure, not a trajectory. Canonical identities are mapped to that state's actual residue numbering; missing coordinates are skipped and reported. Use node-only, edge-only, or combined view. Click a table identity to focus a residue, and hover an edge to inspect its identities and scores. Unselected residues remain gray; gray is not an inferred zero score. Legacy states need an explicit residue mapping to be colored safely.

## Files and Reproducibility

Each run creates:

```text
clusters/<cluster_id>/analyses/model_pair_hotspots/<analysis_id>/
  analysis_metadata.json
  analysis.npz
```

Metadata records model/sample names and IDs, file SHA256 hashes, reference membership, excluded/used frame counts, gauge, edge union, weighting, bootstrap settings and cautions. Atomic directory publication prevents partially written results being listed. Old analyses are not rewritten.

NPZ arrays include residue keys, `(E,2)` edges, separate node/edge parameter RMS; in ensemble mode, A/B means, separation, CI bounds, observed frame and block counts, sample IDs and per-sample preference/coverage. Unsupported values use NaN (JSON null). No all-frame-by-edge energy table or bootstrap draws are sent to the browser. Computation reduces one vectorized frame-to-block term at a time; reference block summaries are retained in memory for bootstrap, so large graphs with very short blocks can still be memory-intensive. Increase block length for scientifically appropriate reasons, or analyze fewer replicas at a time; do not choose blocks solely to optimize speed.

Info buttons expose full metadata; Delete removes the analysis. System health cleanup removes analyses whose endpoint model, reference MD, or preference sample was deleted.

For mutant projection, omit the mutated residue and its incident terms when comparing unchanged residues. This analysis does not implement mutation design or infer changes in the energy function of a mutant.
