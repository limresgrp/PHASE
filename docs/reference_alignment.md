# Reference Structures and Residue Correspondence

## Purpose

A new PHASE system has one reference PDB, one reference protein chain and one
ordered residue selection. Additional PDBs are mapped to that reference by
sequence alignment. Their coordinates and PDB numbering are **not** modified.
This establishes residue correspondence, not structural equivalence: active
and inactive structures may have different conformations despite identical
sequence identity.

Existing systems without a reference remain legacy systems. Their selections,
offsets, descriptors, clusters and analyses are not automatically migrated or
reinterpreted. To adopt reference-based ingestion, create a new system. Changing
the reference or its selection after any states or clusters exist is refused.

## Web Workflow

1. Create a system with exactly one reference PDB. Choose a chain or leave it
   blank to select the longest protein chain. Enter a single MDAnalysis residue
   selection, such as `protein` or `resid 20 to 200`. The selection is evaluated
   on reference PDB numbering and restricted to the chosen chain.
   Chain choice is optional; it distinguishes protein copies for sequence
   alignment, not a replacement for the MDAnalysis selection. Common MD
   protonation/disulfide names (for example `HID`, `HIE`, `HIP`, `CYX`, `CYM`)
   use their parent amino acid for sequence alignment while the original residue
   names and coordinates remain unchanged. Selected residues still require CA.
2. Open the system's **Reference Structure** panel. Inspect the chosen chain,
   sequence length and residue map, including insertion codes and selected keys.
3. Choose **Add state / replica**, provide a name and choose the reference as
   topology. Optionally supply a matching trajectory. **Build descriptors now**
   is enabled by default, including for PDB-ID imports.
   Without a trajectory the PDB is a single frame. With a trajectory, only its
   frames are used, never an extra frame from the PDB.
4. For another experimental structure choose **Upload another PDB**, then
   **Preview alignment(s)**. Inspect identity, coverage, gaps and missing selected
   residues before confirming state creation. If chains tie, enter the desired
   target chain and preview again; `(blank)` explicitly chooses an unnamed chain.
   Manual PDB uploads undergo the same sequence correspondence as PDB-ID imports;
   they need not match the reference atom count, residue numbering or atom order.
   If you choose a PDB in the optional trajectory field while using the reference,
   the web UI switches to the aligned PDB upload workflow and requests a preview.
   The shared ingestion service also aligns such self-contained PDB inputs when
   no explicit uploaded or existing-state topology was chosen. Multi-model PDBs
   use their own frames, without adding a reference frame.
   Trajectories attached to an explicitly uploaded topology or an existing state
   must still match that state's atoms; sequence alignment cannot reorder their
   coordinate arrays. Reopening the creation dialog clears previous file inputs.
5. For another trajectory using the same topology, choose **New replica of ...**.
   The topology is shared, but descriptors and trajectories belong to the new
   state. Deleting the source state does not delete a replica's topology.
6. Existing state cards still support explicit trajectory replacement and
   descriptor building. Prefer creating replicas when the original data should
   remain available. Replacement invalidates that state's descriptors and
   metastable assignments, even when immediate descriptor building is disabled.
   Rebuild/reassign MD samples and rerun downstream analyses when replacing data;
   existing cluster/model/analysis snapshots are not silently rewritten.

### Assigning States to Existing Clusters

During state creation or PDB-ID import, select one or more completed runs under
**Also assign to cluster runs**. Descriptors are built first; only successfully
created, descriptor-ready states are submitted for assignment. Each selected
cluster gets a background job using one compute worker, so one RQ worker is
sufficient. Follow **View assignment progress** to see status and errors.

After uploading, use **System > Reference Structure > Assign states to clusters**.
Choose descriptor-ready states and cluster runs, then **Run assignment**. Results
are MD samples in **Potts Models > Samples**. The older entry point remains
**Potts Models > Samples > New sample > Assign MD states**.

Incomplete structures can be built and assigned without reconstructing loops.
Available residues are extracted normally; absent reference residues retain their
canonical columns and are marked missing. See the mask semantics below.

An aligned state's per-state filter/offset controls are hidden: the system
reference selection is authoritative. Legacy state controls remain available.
Expand **Alignment** in a state card for detailed reference-to-target rows.

## Console Workflow

`./scripts/phase_console.sh` asks for a reference PDB, chain and residue selection
when creating a new system. Its system menu includes reference inspection,
adding a state, PDB-ID import, trajectory attachment/replicas and alignment
inspection. Normal additions no longer ask for offsets in reference systems.

Equivalent commands, using your shared data root:

```bash
python -m phase.scripts.offline_system --root "$PHASE_DATA_ROOT" create-system \
  --project-id PROJECT --name SYSTEM --use-slug-ids \
  --reference /path/reference.pdb --chain-id A --selection 'resid 20 to 200'

python -m phase.scripts.offline_system --root "$PHASE_DATA_ROOT" add-state \
  --project-id PROJECT --system-id SYSTEM --state-id active \
  --name Active --use-reference --traj /path/active.xtc --copy-traj

python -m phase.scripts.offline_system --root "$PHASE_DATA_ROOT" add-state \
  --project-id PROJECT --system-id SYSTEM --state-id experimental \
  --name Experimental --pdb /path/experimental.pdb --chain-id B

python -m phase.scripts.offline_system --root "$PHASE_DATA_ROOT" add-state \
  --project-id PROJECT --system-id SYSTEM --state-id replica2 \
  --name Replica2 --source-state-id active --traj /path/replica2.xtc --copy-traj

python -m phase.scripts.offline_system --root "$PHASE_DATA_ROOT" attach-traj \
  --project-id PROJECT --system-id SYSTEM --state-id active \
  --traj /path/replacement.xtc --slice-spec ::10

python -m phase.scripts.offline_system --root "$PHASE_DATA_ROOT" alignments \
  --project-id PROJECT --system-id SYSTEM

python -m phase.scripts.offline_system --root "$PHASE_DATA_ROOT" alignment \
  --project-id PROJECT --system-id SYSTEM --state-id experimental
```

`add-state` normally builds descriptors. `--no-build-descriptors` only registers
the mapped structure without descriptors; partial structures can also build
descriptors. `attach-traj` copies raw
coordinates into the shared data root by default. `--no-build-descriptors`
defers rebuilding. Always use `--copy-traj` on `add-state` if Docker must access
the trajectory; an external host path is not automatically visible in a container.

Programmatic `create-system` without `--reference` remains available for legacy
clients. `reference --pdb FILE --chain-id A --selection ...` configures an empty
system; omit `--pdb` to inspect its existing reference. Manual `--resid-shift` is
supported only on legacy systems, with a warning.

## Mapping and Scientific Provenance

Protein sequences are extracted from residues present in the first coordinate
model and having a C-alpha atom. SEQRES is not substituted for coordinate
identity. The reference therefore describes its **observed coordinate sequence**,
not an asserted complete UniProt sequence. Missing sequence positions that are
absent from the reference itself cannot be inferred from these coordinates.
Target deletions are relative to this observed reference sequence.

Each residue identity includes chain, integer PDB number, insertion code and
residue name. `A:45`, `A:45A` and `A:46` are distinct. Canonical keys use reference
numbering (`res_45`, `res_45A`, `res_46`) and retain reference coordinate order,
not alphabetical order. The selection is an ordered subset of these keys.

Alignment uses Biopython global pairwise alignment with **BLOSUM62**, gap-open
score **-10**, and gap-extension score **-0.5**. The independent core API accepts
alternative gap penalties. Reference chain selection chooses the longest
observed protein chain (file order breaks equal-length reference ties). Target
selection compares global alignment scores; equal best target chains require
an explicit chain choice. Alternate optimal within-chain mappings are flagged
and classified as warning; inspect correspondence before using them.

Definitions:
- **Identity:** identical residues divided by aligned residue pairs (gaps excluded).
- **Similarity:** positive BLOSUM62 substitutions divided by aligned pairs.
- **Reference coverage:** aligned pairs divided by reference coordinate residues.
- **Target coverage:** aligned pairs divided by target coordinate residues.
- **Insertions/deletions:** zero-based sequence positions on target/reference,
  respectively; extra target identities and per-reference correspondence rows
  identify the affected PDB residues.
- **Missing selected residues:** selected canonical keys with no target residue.
- **Quality score:** identity multiplied by reference coverage, a neutral summary.

Thresholds live centrally in `phase/io/alignment.py`:
- Insufficient correspondence fails below 50% reference coverage or 30% identity.
- Excellent: identity and reference coverage at least 98%, no missing selections.
- Good: identity and reference coverage at least 95%, at most two missing selections.
- Warning: other successful alignments, including alternate optimal mappings.

These labels describe correspondence quality, not functional or structural
equivalence. RMSD-based superposition is not used to determine residue identity.
Sequence substitutions in selected residues are also listed explicitly. Review
them before combining conformational descriptors: residue correspondence does
not make different side-chain chemistries or torsional alphabets equivalent.

Every target stores reference and target SHA-256 fingerprints, scoring settings,
chain choices, identity/coverage, insertions/deletions, forward and reverse
residue maps, selected/missing keys, quality and upload/download provenance.
Detailed sidecars are authoritative; the API exposes alignment information
without requiring clients to locate files. Editing a topology after alignment
causes descriptor generation to fail rather than reuse a stale mapping.

## Missing Residues and Descriptor Ordering

Partial alignments retain the **complete ordered canonical selection**. Absent
residues have null correspondence entries and NaN descriptor vectors; following
columns never shift. State/descriptor metadata list `missing_residue_keys`.
The same workflow is used by the console and webserver.

Cluster assignment predicts only finite descriptor observations. Missing
coordinates receive label **-2**, distinct from halo/unassigned label **-1**.
Sample NPZs persist a `(frames, residues)` boolean `missing_mask`. A partial
frame is retained: it is not an invalid SA frame, a new cluster, or an imputed
conformation. JS marginals use observed frames per residue; edge distributions
use frames where both endpoints exist. No observations means unavailable
(NaN/null), not zero JS or a uniform fabricated distribution.

**Energy semantics:** missing node fields and every coupling touching a missing
residue contribute exactly zero, for both single-model and delta-energy scores.
Energy component NPZs save observation masks; analysis metadata record node/edge
coverage and missing residue indices. Curves are marked partial in the web UI.
This is a partial Hamiltonian in the analysis gauge, **not** marginalization over
unknown conformations. Partial scores depend on coverage and gauge; compare
structures using matched residue coverage when making quantitative claims.
There is no automatic rescaling to full-protein energy.

Commitment/frustration on unavailable nodes/edges is unavailable (NaN), and
normalization excludes those values. Complete-system results keep the same
semantics. Missing coordinates are never treated as an extra categorical state
for fitting or as a warm-start conformation: use complete training/initialization
frames. Specialized analyses that require complete conformations are not made
missing-aware merely by attaching a mask.

For existing incomplete states, click **Build**, then **Assign states to clusters**,
then rerun energies/JS to obtain the new masks. Old analysis snapshots are not
rewritten automatically.

Mapped residues are resolved to exact topology `resindex` selections using full
chain/number/insertion identities. Descriptor NPZ keys, state metadata, clustering
and assignments retain canonical order. Existing absent-torsion conventions for
chemically undefined dihedrals (e.g. side chains without chi2) are unchanged.
For new reference-based descriptors, backbone neighbors follow topology order
within each chain/segment, not integer `resid +/- 1`. A first-frame C--N distance
between 0.5 and 1.9 angstrom is required (minimum-image distance when a usable box
exists), preventing torsions across coordinate gaps and correctly handling
`45 -> 45A -> 46`. Legacy extraction keeps its existing MDAnalysis behavior.
The policy is recorded in descriptor metadata and reused for backmapping indices.

## Importing PDB IDs

The web panel supports up to 20 IDs per request. Preview the alignments, choose
entries individually or select high-confidence mappings (including gaps), then create the
selected states. Errors are reported per ID, not as an all-or-nothing batch.
Downloaded structures are PDB-only states. Web imports build descriptors by
default; uncheck **Build descriptors now** to keep a structure without descriptors.

```bash
python -m phase.scripts.offline_system --root "$PHASE_DATA_ROOT" add-pdb-ids \
  --project-id PROJECT --system-id SYSTEM --pdb-ids 3EML,5G53,8JWY --preview

python -m phase.scripts.offline_system --root "$PHASE_DATA_ROOT" add-pdb-ids \
  --project-id PROJECT --system-id SYSTEM --pdb-id-file ids.txt --chain-id A
```

The CLI defaults to no descriptor calculation for batch imports; add
`--build-descriptors` explicitly if desired. `(blank)` selects an unnamed chain.
Alignment/batch inspection commands print summaries by default; use `--json`
to export full standalone metadata. Web creation rechecks previewed download
fingerprints before confirming imported structures.
ID validation requires a digit followed by three alphanumeric characters.
Downloads use fixed RCSB HTTPS URLs, a timeout and a 50 MiB size limit. PDB is
requested first; a 404 falls back to mmCIF, preserving the original mmCIF and
converting a PDB-compatible topology for the current descriptor/visualization
stack. Structures exceeding PDB chain/number/atom representation limits fail
with an explicit message; identifiers are never silently truncated or renumbered.
Uploaded reference/target files currently use PDB format.
Converted mmCIF coordinates use PDB serialization precision; the unmodified
original mmCIF remains available under the recorded provenance path.

## Visualization Compatibility

The new alignment panels and Distance Analysis 3D viewer use explicit
chain/number/insertion mappings. Older analysis-specific 3D viewers that infer
residue identity from numeric keys have not all been migrated in this pass.
Use the reference topology there for ordinary integer-numbered selections;
arbitrary-target coloring/picking and insertion-coded selections in those
legacy viewers require a follow-up migration. Do not interpret a color on an
unmapped target as a validated residue correspondence.

Alignment and import are shared core services, not duplicated web/console logic.
The importer runs synchronously in a worker thread for bounded small batches;
descriptor jobs for large imports should be scheduled separately. No live
internet access is used by unit tests.

## Files and APIs

```text
reference/reference.pdb
reference/reference_metadata.json
reference/topologies/<sha256>.pdb
reference/originals/<sha256>.cif     # when downloaded mmCIF was converted
states/<state>/alignment_metadata.json
states/<state>/state_metadata.json
states/<state>/descriptors.npz
states/<state>/descriptor_metadata.json
states/<state>/trajectory*.xtc
```

State topology paths may point at the reference or a shared content-addressed
topology. State-owned trajectories and descriptors remain separate. Cloning
copies reference data and preserves relative mappings. Shared files survive
deletion while referenced by another state.

API routes below the existing project/system prefix:
- `POST /reference`: one multipart `pdb`, optional `chain_id`, `selection`.
- `GET /reference`: download reference PDB.
- `POST /reference/preview`: preview target `pdb`, optional `chain_id`.
- `POST /states`: `name`, plus `pdb`, `use_reference=true` or `source_state_id`;
  optional trajectory, chain, slice and descriptor-build flag.
- `POST /states/{state_id}/trajectory`: explicit attachment/replacement.
- `POST /states/from-pdb-ids`: IDs, preview/create flag and optional per-ID chains.
- `GET /alignments` and `GET /states/{state_id}/alignment`: alignment provenance.

The old system-creation route remains; initial multi-PDB ingestion requires
`legacy_import=true`. New UI code uses the separate reference route.

Implementation references:
[Biopython pairwise alignment](https://biopython.org/docs/latest/Tutorial/chapter_pairwise.html),
[RCSB file downloads](https://www.rcsb.org/docs/programmatic-access/file-download-services).
