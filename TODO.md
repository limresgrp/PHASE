# PHASE: Reference PDB alignment and multi-structure state ingestion

## 1. Objective

Extend PHASE so that a system has one explicit **reference PDB** defining:

1. the canonical protein sequence and residue ordering;
2. the canonical residue-selection space used by PHASE;
3. the canonical order of descriptor residues.

All additional structures/states must be aligned to this reference before descriptors are generated.

The user workflow must become:

### Initial system creation

1. User creates a system.
2. User uploads exactly one reference PDB.
3. PHASE parses the reference PDB and stores its sequence/residue identity map.
4. User selects residues on the reference PDB.
5. User optionally uploads one trajectory for that state.
6. The trajectory becomes the first state.

### Adding additional states

A user can either:

A. upload a trajectory only and attach it to an existing state/reference PDB; or

B. upload another PDB:
   - PHASE automatically aligns it to the system reference PDB;
   - PHASE creates the residue mapping automatically;
   - the state uses exactly the same reference residue selection and ordering;
   - insertions/deletions/missing residues are handled explicitly;
   - alignment metadata is stored with the state.

### Future extension

Support:

```text
PDB IDs -> automatic download -> alignment to reference -> state creation
```

This must be implemented as a second layer on top of the same alignment service, not as a separate special-case pathway.

---

# 2. Important design principle

Do NOT continue treating `resid_shift` as the primary mechanism.

The existing implementation has:

```python
DescriptorState.resid_shift
DescriptorState.residue_keys
DescriptorState.residue_mapping
```

and `macro_states.py` currently applies a global integer shift to keys.

That is insufficient for structures with:

- insertions;
- deletions;
- chain numbering changes;
- insertion codes such as `45A`;
- missing residues;
- different PDB numbering schemes;
- multiple alignment blocks.

`resid_shift` should remain for backward compatibility, but automatically generated alignment metadata must become authoritative whenever a structure has been aligned to a reference.

The new model should support arbitrary residue-to-reference mappings.

---

# 3. Data model

Extend `SystemMetadata`.

Current system metadata contains:

```python
residue_selections
residue_selections_mapping
states
```

Add:

```python
reference_pdb_file: Optional[str]
reference_chain_id: Optional[str]
reference_sequence: Optional[str]
reference_residue_keys: List[str]
reference_residue_records: List[Dict[str, Any]]
reference_alignment_metadata_file: Optional[str]
reference_selection: Optional[str]
```

Suggested semantics:

### `reference_pdb_file`

Relative path from the system directory.

Example:

```text
reference/reference.pdb
```

### `reference_chain_id`

The chain used as the canonical protein chain.

### `reference_sequence`

Canonical one-letter amino-acid sequence corresponding to `reference_residue_records`.

### `reference_residue_keys`

Canonical PHASE residue keys in reference order.

Example:

```json
[
  "res_1",
  "res_2",
  "res_3",
  ...
]
```

Do not infer canonical order from lexicographical sorting.

The ordering must always be the order defined by the reference structure.

### `reference_residue_records`

Store sufficient information to recover exactly which PDB residue each key corresponds to.

Example:

```json
[
  {
    "index": 0,
    "residue_key": "res_1",
    "chain_id": "A",
    "resname": "MET",
    "resid": 1,
    "insertion_code": ""
  }
]
```

### `reference_alignment_metadata_file`

Optional system-level alignment metadata.

Example:

```text
reference/reference_metadata.json
```

---

# 4. State-level alignment metadata

Extend `DescriptorState` with:

```python
reference_alignment: Optional[Dict[str, Any]]
alignment_status: Optional[str]
```

Recommended `alignment_status` values:

```text
reference
aligned
partial
failed
```

`reference` means this state is the system reference.

`aligned` means all required selected reference residues were mapped.

`partial` means alignment succeeded, but some reference residues are absent in the target.

`failed` means automatic alignment could not be established.

---

# 5. Alignment metadata schema

For every non-reference PDB state store something like:

```json
{
  "method": "global_sequence_alignment",
  "reference_chain": "A",
  "target_chain": "A",
  "score": 123.4,
  "identity": 0.97,
  "coverage_reference": 0.98,
  "coverage_target": 0.99,

  "reference_length": 309,
  "target_length": 307,

  "matches": 304,
  "mismatches": 2,
  "insertions": [
    {
      "target_start": 145,
      "target_end": 146,
      "reference_position": 145
    }
  ],
  "deletions": [
    {
      "reference_start": 201,
      "reference_end": 202
    }
  ],

  "residue_mapping": {
    "res_1": "A:1",
    "res_2": "A:2",
    "res_3": null
  },

  "target_to_reference": {
    "A:1": "res_1",
    "A:2": "res_2"
  },

  "selected_reference_residues": [
    "res_1",
    "res_2"
  ],

  "mapped_selected_residues": [
    "res_1",
    "res_2"
  ],

  "missing_selected_residues": [],

  "extra_target_residues": [
    "A:145A"
  ]
}
```

The exact fields may be adjusted, but the following information is mandatory:

- alignment method;
- reference chain;
- target chain;
- sequence identity;
- reference coverage;
- target coverage;
- insertion positions;
- deletion positions;
- reference -> target residue mapping;
- target -> reference residue mapping;
- missing selected residues;
- extra target residues.

This metadata is scientific provenance, not merely debugging information. It will later be used to rank structures for the Cryo-EM validation.

---

# 6. Residue identity representation

Do not represent a residue only by integer `resid`.

Use a canonical residue identifier:

```python
@dataclass(frozen=True)
class ResidueRef:
    chain_id: str
    resid: int
    insertion_code: str = ""
    resname: Optional[str] = None
```

Provide conversion:

```python
ResidueRef -> "A:45"
ResidueRef -> "A:45A"
```

This is necessary because PDB insertion codes are common and a global integer offset cannot represent them.

---

# 7. Alignment algorithm

Create a dedicated alignment module, e.g.:

```text
phase/io/alignment.py
```

or:

```text
phase/common/alignment.py
```

The module must be independent of the webserver and console.

Suggested public API:

```python
align_pdb_to_reference(
    reference_pdb: Path,
    target_pdb: Path,
    *,
    reference_chain: Optional[str] = None,
    target_chain: Optional[str] = None,
) -> AlignmentResult
```

and:

```python
build_reference_residue_map(
    reference_pdb: Path,
    *,
    chain_id: Optional[str] = None,
) -> ReferenceStructure
```

and:

```python
map_selection_to_target(
    reference_selection: Sequence[ResidueRef],
    alignment: AlignmentResult,
) -> SelectionMappingResult
```

---

# 8. Chain selection

The algorithm must not silently align the wrong chain.

For a PDB with multiple protein chains:

1. identify protein chains;
2. if the user explicitly specifies a chain, use it;
3. for the reference, if no chain is specified:
   - choose the longest protein chain;
   - store the choice;
4. for a target structure:
   - first compare sequence similarity against the reference sequence;
   - select the best matching protein chain;
   - if two chains are equally plausible, mark the alignment ambiguous and require user confirmation.

Do not choose a chain simply because it is chain A.

---

# 9. Sequence extraction

Use BioPython or the existing project molecular-analysis stack.

The sequence must correspond to the actual PDB residues that can be selected and used by descriptors.

Do not use the SEQRES sequence as the sole basis for the structure mapping.

The residue table must preserve:

```text
PDB chain
PDB resid
insertion code
resname
sequence position
```

because the structure can contain missing residues.

Distinguish:

- sequence residue exists in canonical sequence;
- residue is present in the PDB coordinates.

This distinction is essential for Cryo-EM structures.

---

# 10. Sequence alignment

Use a standard classical sequence alignment.

Do NOT invent a custom offset algorithm.

Preferred default:

Needleman-Wunsch global alignment for the chosen chains.

Use a reasonable protein scoring matrix, e.g. BLOSUM62, with configurable gap penalties.

The implementation should expose:

```python
alignment_score
identity
similarity
coverage_reference
coverage_target
```

The default should align the full protein chain sequence.

For structures with clearly missing coordinate segments, sequence alignment still works; mapped residues can then be marked as absent from the coordinate set.

---

# 11. Mapping reference residues to target residues

The alignment must generate an explicit residue correspondence.

For every reference residue:

```text
reference residue -> target residue or None
```

Example:

```text
reference:
1 2 3 4 5 6 7 8 9

target:
1 2 3 5 6 7 8 9 10
```

then:

```text
ref 1 -> target 1
ref 2 -> target 2
ref 3 -> target 3
ref 4 -> None
ref 5 -> target 5
...
```

Do not renumber the target structure to pretend that the deletion never happened.

The reference coordinate space remains canonical.

---

# 12. Canonical selection

This is the most important behavioral requirement.

After the user selects residues on the reference:

```text
reference selection = [res_20, res_25, res_30, res_40]
```

every target structure must use exactly that ordered list.

For a target with a deletion:

```text
res_20 -> A:20
res_25 -> A:25
res_30 -> None
res_40 -> A:40
```

the target state's selected residue list must remain:

```text
[res_20, res_25, res_30, res_40]
```

with an explicit missing entry for `res_30`.

Do not silently remove `res_30` and do not reorder the remaining residues.

This is necessary because downstream descriptor arrays, clustering, Potts models, and state comparisons depend on consistent residue ordering.

---

# 13. Descriptor behavior for missing residues

The implementation must make a conscious decision here.

Preferred behavior:

- a target residue absent from the PDB cannot produce a valid dihedral descriptor;
- represent it as missing/NaN at descriptor-building time;
- record it in `missing_selected_residues`;
- do NOT fabricate a descriptor;
- do NOT shift subsequent residue indices.

However, any existing downstream clustering/Potts code that assumes a rectangular descriptor matrix must be audited.

The safest initial implementation is:

1. allow alignment to succeed with missing residues;
2. preserve the canonical selection;
3. mark missing residues explicitly;
4. refuse descriptor construction if the selected residue is required by the descriptor pipeline and cannot be represented.

The error must explain exactly which residues are missing.

Do not silently produce a different topology.

---

# 14. Automatic selection transfer

The existing UI allows a user to type a residue selection.

This must change conceptually.

The selection belongs to the **reference system**, not to each individual state.

Move toward:

```python
SystemMetadata.reference_selection
```

instead of requiring users to repeat:

```text
residue_selection
resid_shift
```

for every state.

Backward compatibility:

- retain `DescriptorState.residue_selection`;
- when an old system is loaded, continue honoring it;
- new systems should use reference selection as the source of truth.

---

# 15. Reference-first workflow

Change system creation.

Current code:

```python
create_system_with_descriptors(...)
    pdb_files: List[UploadFile]
```

currently allows multiple independent PDBs to be uploaded directly as states.

New behavior:

```text
Create system
    -> upload reference PDB
    -> establish reference structure
    -> create first state from reference
```

The backend should reject multiple PDBs in the initial reference upload endpoint unless the request explicitly uses the legacy compatibility path.

Prefer a new API rather than overloading the old endpoint.

For example:

```http
POST /projects/{project_id}/systems/{system_id}/reference
```

multipart:

```text
pdb
chain_id (optional)
```

Response:

```json
{
  "reference_pdb_file": "...",
  "reference_chain_id": "A",
  "reference_sequence": "...",
  "reference_residue_records": [...]
}
```

---

# 16. First state creation

After creating the reference, the frontend should offer:

```text
Use reference structure as first state
```

This creates:

```text
state_id = generated UUID
name = user supplied name
pdb_file = reference PDB
alignment_status = "reference"
reference_alignment = null
```

The user can then attach a trajectory.

The first state should not require an independent PDB upload.

---

# 17. Add state with another PDB

Replace the current simple endpoint behavior:

```http
POST /projects/{project_id}/systems/{system_id}/states
```

with support for:

```text
name
pdb
source_state_id (optional)
alignment_mode
chain_id (optional)
```

If `pdb` is supplied:

1. save uploaded PDB to temporary location;
2. load system reference;
3. detect best target chain;
4. run sequence alignment;
5. compute residue mapping;
6. transfer canonical reference selection;
7. write alignment metadata;
8. only then finalize the state.

If alignment fails, do not create a partially registered state.

---

# 18. Add trajectory to an existing state

This should remain separate.

The trajectory is assumed to correspond to the PDB already attached to the state.

Do not realign a trajectory.

However, validate:

- topology compatibility;
- atom count;
- residue count;
- chain/residue correspondence where practical.

The PDB associated with the state is the topology used for trajectory descriptor generation.

---

# 19. Multiple trajectories for the same structural state

A user should be able to do:

```text
reference PDB
   |
   +-- trajectory A -> state "active"
   +-- trajectory B -> state "active_replica_2"
   +-- trajectory C -> state "active_replica_3"
```

There are two acceptable representations:

### Preferred

Each trajectory is a distinct PHASE state, sharing the same PDB identity/alignment metadata.

Example:

```text
state active_replica_1
state active_replica_2
state active_replica_3
```

All have:

```text
same reference alignment
same residue order
same topology
different trajectory
```

This keeps the existing PHASE abstraction intact.

Do not duplicate the physical PDB unnecessarily.

Instead allow state metadata to reference the same topology.

---

# 20. Deduplicate shared PDBs

When several states use the same structural PDB:

```text
states/
  active_replica_1/
  active_replica_2/
```

do not needlessly copy the same PDB multiple times.

A practical option is:

```text
structures/
  reference.pdb
  pdb_<hash>.pdb
```

and each state stores a relative reference such as:

```json
"pdb_file": "structures/pdb_abc123.pdb"
```

The state still owns its descriptor and trajectory directory.

Do not introduce this optimization if it requires a large refactor; correctness comes first.

---

# 21. Alignment metadata visualization

The web UI must display alignment quality per state.

At minimum:

```text
Alignment
✓ aligned

Sequence identity: 97.4%
Reference coverage: 98.7%
Target coverage: 99.0%

Matches: 302
Mismatches: 5
Insertions: 2
Deletions: 3
Missing selected residues: 1
```

Provide expandable details showing:

```text
Reference residue | Target residue | Status
---------------------------------------------
res_20             | A:20           | match
res_21             | A:21           | match
res_22             | —              | deletion
res_23             | A:23           | match
```

This is particularly important for the Cryo-EM analysis.

---

# 22. Alignment-quality classification

Create a deterministic quality classification.

Suggested:

```text
excellent:
    identity >= 0.98
    coverage_reference >= 0.98
    missing_selected == 0

good:
    identity >= 0.95
    coverage_reference >= 0.95
    missing_selected <= 2

warning:
    alignment succeeds but thresholds are worse

failed:
    alignment ambiguous or insufficient coverage
```

Do not hard-code these thresholds throughout the application.

Define them centrally:

```python
ALIGNMENT_THRESHOLDS = ...
```

Make them configurable later.

The metadata should contain both the quantitative metrics and the resulting quality label.

---

# 23. Cryo-EM comparison use case

The metadata must allow later selection/ranking of structures.

For example:

```python
alignment_quality_score
missing_selected_residue_count
identity
coverage_reference
insertion_count
deletion_count
```

A downstream analysis can then say:

```text
best structures to compare:
    identity high
    no missing selected residues
    no unusual insertions/deletions

least safe structures:
    low coverage
    missing selected residues
    large alignment gaps
```

Do NOT hard-code "safe" / "unsafe" semantics into the alignment module.

Expose a neutral `quality` and metrics.

The scientific analysis layer can interpret them.

---

# 24. PDB-ID ingestion: second step

Implement a reusable backend function:

```python
download_pdb_structure(pdb_id: str, destination: Path) -> Path
```

Do not couple it to the web UI.

Recommended API:

```http
POST /projects/{project_id}/systems/{system_id}/states/from-pdb-ids
```

Request:

```json
{
  "pdb_ids": [
    "3EML",
    "5G53",
    "8JWY"
  ],
  "name_template": "PDB_{pdb_id}"
}
```

Backend behavior:

For each ID:

1. validate format `[0-9][A-Za-z0-9]{3}`;
2. download PDB/mmCIF from RCSB;
3. determine protein chain;
4. align against system reference;
5. calculate metadata;
6. create state;
7. save original downloaded structure;
8. optionally calculate descriptors.

Return one result per PDB:

```json
{
  "pdb_id": "5G53",
  "state_id": "...",
  "status": "created",
  "alignment_quality": "excellent",
  "identity": 1.0,
  "missing_selected_residue_count": 0
}
```

Failures must be per-PDB, not all-or-nothing.

Example:

```json
{
  "pdb_id": "XXXX",
  "status": "failed",
  "error": "No suitable protein chain found"
}
```

---

# 25. CLI / phase_console

The current `state_menu` contains:

```text
add-state|Add state (PDB, optional trajectory)
```

Change this menu to:

```text
reference|Set/view reference structure
add-state|Add state
add-pdb|Add state from PDB file
add-pdb-ids|Add states from PDB IDs
attach-traj|Attach trajectory to state
alignment|Inspect alignment metadata
delete-state|Delete state
```

For a new system:

```text
Create system
```

should immediately ask:

```text
Reference PDB path:
Reference chain [auto]:
```

Then:

```text
Reference loaded.
Sequence length: 309
Protein residues with coordinates: 304
```

Then ask for the canonical residue selection.

Example:

```text
Residue selection on reference:
```

This selection becomes the system-level selection.

---

# 26. CLI: adding a PDB state

Command:

```bash
python -m phase.scripts.offline_system add-state \
    --project-id ... \
    --system-id ... \
    --state-id ... \
    --name ... \
    --pdb ...
```

must automatically align to the reference.

Do NOT require:

```text
--resid-shift
```

for normal use.

Retain:

```text
--resid-shift
```

only as a legacy/manual override, with a warning:

```text
WARNING: manual residue shifts bypass reference alignment.
```

---

# 27. CLI: PDB IDs

Add:

```bash
python -m phase.scripts.offline_system add-pdb-ids \
    --project-id PROJECT \
    --system-id SYSTEM \
    --pdb-ids 3EML,5G53,8JWY
```

or accept a text file:

```bash
--pdb-id-file pdb_ids.txt
```

Allow both.

Optional:

```bash
--build-descriptors
```

Default should be false for many PDB IDs because descriptor building can be expensive.

---

# 28. CLI: alignment inspection

Add:

```bash
python -m phase.scripts.offline_system alignment \
    --project-id ...
    --system-id ...
    --state-id ...
```

Output:

```text
State: 5G53
Reference: A2A_reference
Target chain: A

Identity             100.0%
Reference coverage    99.7%
Target coverage       99.7%
Insertions                0
Deletions                 1
Missing selected           0

Quality: excellent
```

Also add a system-wide summary:

```bash
python -m phase.scripts.offline_system alignments \
    --project-id ...
    --system-id ...
```

with columns:

```text
state | pdb | identity | coverage | missing_selected | insertions | deletions | quality
```

This is the table we will ultimately use for Cryo-EM curation.

---

# 29. Web UI: system page

The current `StateCards.jsx` is state-centric and has manual:

```text
Residue filter
Residue shift offset
```

These should no longer be the main controls for new systems.

Add a dedicated:

```text
Reference Structure
```

panel above the states list.

It should show:

```text
Reference PDB: 5MZJ
Chain: A
Sequence length: 309
Selected residues: 87
```

Buttons:

```text
Download reference
Change reference
View residue mapping
```

Changing the reference after states exist must be considered a destructive operation.

For the first implementation, either:

- disallow changing the reference once states/descriptors exist; or
- require explicit re-alignment/rebuild of all states.

Do not silently change the canonical coordinate system.

---

# 30. Web UI: add state

Replace the current simple PDB upload modal/flow with:

```text
Add state

( ) Use existing reference PDB
( ) Upload new PDB

State name:
PDB file:
Chain: Auto
Trajectory: optional
```

If "existing reference PDB":

- create state using reference topology;
- no alignment needed.

If "upload new PDB":

- automatically align;
- display alignment summary before final confirmation.

Example:

```text
Alignment found

Reference chain: A
Target chain: A

Identity: 97.8%
Reference coverage: 98.4%
Target coverage: 99.0%

Insertions: 1
Deletions: 2
Missing selected residues: 0

Quality: GOOD

[Create state]
[Cancel]
```

The state should not be created before successful alignment.

---

# 31. Web UI: attach another trajectory

Each state should have:

```text
Add trajectory
```

This creates a new state if the user chooses:

```text
Create new replica/state from this structure
```

or updates the existing state if explicitly requested.

The cleaner default is:

```text
Add trajectory as new replica
```

because PHASE treats each trajectory as a state/ensemble.

---

# 32. Web UI: PDB-ID batch ingestion

Add a button:

```text
Import PDB IDs
```

Dialog:

```text
Enter PDB IDs
----------------
3EML
5G53
8JWY
8JWZ

[Preview alignments]
```

Backend returns alignment results.

Display a table:

| PDB | Identity | Coverage | Missing selected | Ins | Del | Quality |
|-----|----------|----------|------------------|-----|-----|---------|

Allow:

```text
[Select all safe]
[Create selected states]
```

This is much more useful than blindly creating every structure.

---

# 33. API serialization

Update the system serializer so every state exposes:

```json
{
  "state_id": "...",
  "name": "...",
  "pdb_file": "...",
  "trajectory_file": "...",

  "alignment_status": "aligned",

  "alignment": {
    "method": "global_sequence_alignment",
    "reference_chain": "A",
    "target_chain": "A",
    "identity": 0.978,
    "coverage_reference": 0.984,
    "coverage_target": 0.990,
    "insertions": 1,
    "deletions": 2,
    "missing_selected_residue_count": 0,
    "quality": "good"
  }
}
```

Do not force the frontend to read the raw JSON alignment file.

The JSON sidecar remains the persistent source of truth; the API exposes the useful summary.

---

# 34. Preserve existing state metadata compatibility

Existing projects may contain:

```json
"resid_shift": 2
```

and no alignment metadata.

When loading these:

- continue to work exactly as before;
- do not attempt a retroactive alignment automatically;
- mark alignment status as `"legacy"` or null.

For newly uploaded PDBs:

- use the new alignment system.

---

# 35. Descriptor generation changes

Current `build_state_descriptors()` does:

```python
preprocessor.build_single(...)
_apply_residue_shift(...)
```

Change this to:

```python
alignment = state_meta.reference_alignment

selection = resolve_reference_selection(system_meta)

target_selection = map_reference_selection(
    selection,
    alignment,
)

build_result = preprocessor.build_single(
    ...,
    residue_selection=target_selection,
)
```

However, do not produce a target selection by simply converting the reference integers.

The mapping must preserve the target PDB's actual `chain_id`, `resid`, and insertion code.

---

# 36. Critical descriptor ordering invariant

After descriptor construction:

```python
state_meta.residue_keys
```

must always be in the canonical reference order.

Example:

```text
Reference:
res_10
res_11
res_12
res_13

Target actual residues:
A:10
A:11
A:13
A:14
```

The descriptor keys must be:

```text
res_10
res_11
res_12
res_13
```

and not:

```text
res_10
res_11
res_13
res_14
```

`res_12` must be represented as missing.

This is a hard invariant.

---

# 37. Do not conflate sequence alignment and structural alignment

The initial implementation must use:

```text
sequence alignment
```

to determine residue correspondence.

Then optionally calculate a structural RMSD after correspondence is established.

Do not use a structural RMSD superposition as the primary residue mapping method.

Why:

- structural RMSD can produce ambiguous mappings;
- ligand-induced conformational changes should not alter residue identity;
- insertions/deletions need sequence-level handling;
- the target structures may be active/inactive and therefore structurally different.

The alignment metadata can later include:

```text
backbone RMSD
```

as an additional quality metric.

---

# 38. Optional structural validation

After sequence alignment:

1. extract Cα atoms for matched residues;
2. perform Kabsch alignment;
3. calculate backbone/Cα RMSD.

Store:

```json
"ca_rmsd": 1.42
```

This is optional for state creation, but recommended because it helps identify structures that are chemically aligned yet structurally unusual.

Do not reject structures based solely on RMSD.

---

# 39. Testing requirements

Add unit tests for the alignment module.

At minimum:

### Test 1: identical sequence

Reference:

```text
ABCDE
```

Target:

```text
ABCDE
```

Expect 100% identity and one-to-one mapping.

### Test 2: deletion

Reference:

```text
ABCDE
```

Target:

```text
ABDE
```

Expect:

```text
A -> A
B -> B
C -> None
D -> D
E -> E
```

### Test 3: insertion

Reference:

```text
ABCDE
```

Target:

```text
ABXCDE
```

Expect inserted X represented only on target side.

### Test 4: insertion code

Target:

```text
45
45A
46
```

Ensure all three residues are distinct.

### Test 5: sequence mismatch

Same length but substitutions.

Verify identity is correctly calculated.

### Test 6: wrong chain

Reference A matches target B better than target A.

Verify automatic target chain selection.

### Test 7: ambiguous chains

Two target chains with equally good scores.

Verify alignment returns an explicit ambiguous result rather than silently choosing one.

### Test 8: selection propagation

Reference selection:

```text
[res_2, res_4, res_5]
```

target with deletion at res_4.

Verify canonical ordering and missing-residue metadata.

### Test 9: legacy `resid_shift`

Existing systems still build descriptors unchanged.

### Test 10: state creation rollback

If alignment fails, no state directory or metadata remains.

---

# 40. Integration tests

Add an end-to-end test:

```text
create project
create system with reference PDB
create first reference state
attach trajectory
create second state from modified PDB
```

Verify:

```text
same residue order
same reference keys
alignment metadata persisted
descriptor metadata persisted
```

Also test:

```text
PDB-only target
```

because single experimental structures are central to the intended Cryo-EM workflow.

---

# 41. PDB-ID download tests

Do not use the real internet in unit tests.

Mock:

```python
download_pdb_structure()
```

Test:

- valid PDB ID;
- 404;
- malformed ID;
- invalid PDB;
- no protein chains;
- ambiguous chain.

Add one integration test behind an explicit optional marker if desired.

---

# 42. Error handling

Alignment failures must be actionable.

Bad:

```text
Alignment failed
```

Good:

```text
Could not identify a suitable protein chain in target PDB.

Reference chain A:
  length 309

Target protein chains:
  B: length 72
  C: length 84

No chain exceeds the minimum coverage threshold of 0.50.
```

Similarly:

```text
Alignment succeeded, but 5 selected residues are absent from the target structure:
res_124, res_125, res_167, res_168, res_201
```

---

# 43. Recommended code organization

Do not place alignment code in `backend/api/v1/routes/states.py`.

Recommended structure:

```text
phase/
  io/
    alignment.py
    pdb.py

  services/
    project_store.py

  workflows/
    macro_states.py
```

`backend/` should only orchestrate HTTP requests.

`phase_console` should use the same core functions.

This prevents console/webserver behavior from diverging.

---

# 44. Reference PDB download/import provenance

For every reference/target PDB store:

```json
{
  "source": "upload"
}
```

or:

```json
{
  "source": "rcsb",
  "pdb_id": "5G53",
  "download_url": "...",
  "retrieved_at": "..."
}
```

For uploaded files:

```json
{
  "source": "upload",
  "original_filename": "5G53.pdb"
}
```

Do not depend on filenames for scientific identity.

The PDB ID should be explicit metadata for downloaded structures.

---

# 45. File layout

Use a simple layout compatible with the current project structure.

Recommended:

```text
projects/
  <project>/
    systems/
      <system>/
        system.json

        reference/
          reference.pdb
          reference_metadata.json

        structures/
          5G53.pdb
          8JWY.pdb

        states/
          active/
            state_metadata.json
            descriptors.npz
            descriptor_metadata.json
            trajectory.xtc

          active_replica_2/
            state_metadata.json
            descriptors.npz
            descriptor_metadata.json
            trajectory.xtc
```

Do not put alignment information only in `system.json`.

Store the detailed alignment next to the state that it describes.

---

# 46. Metadata documentation

Update:

```text
docs/metadata/system.md
```

with:

- `reference_pdb_file`
- `reference_chain_id`
- `reference_sequence`
- `reference_residue_records`
- `reference_selection`

Update state metadata documentation with:

- `alignment_status`
- `reference_alignment`

Add a new document:

```text
docs/reference_alignment.md
```

documenting:

- sequence alignment;
- residue mapping;
- insertion/deletion semantics;
- missing residues;
- descriptor ordering;
- quality metrics.

---

# 47. Web/API compatibility

Do not break existing API clients.

Existing:

```http
POST /projects/{project_id}/systems
```

may remain as a legacy endpoint.

For new UI code, use:

```http
POST /projects/{project_id}/systems
POST /projects/{project_id}/systems/{system_id}/reference
POST /projects/{project_id}/systems/{system_id}/states
POST /projects/{project_id}/systems/{system_id}/states/{state_id}/trajectory
POST /projects/{project_id}/systems/{system_id}/states/from-pdb-ids
GET  /projects/{project_id}/systems/{system_id}/alignments
GET  /projects/{project_id}/systems/{system_id}/states/{state_id}/alignment
```

The exact routing can follow existing style, but keep the conceptual separation.

---

# 48. Batch PDB-ID workflow should be asynchronous only when necessary

Downloading and aligning 5–20 PDBs is cheap enough to run as a request in many environments.

For larger batches or descriptor generation, use the existing RQ/background-job architecture.

The first implementation may synchronously:

```text
download
align
create state
```

and defer:

```text
descriptor generation
```

This keeps the API simple.

---

# 49. Alignment metadata should be usable outside PHASE

Write a standalone JSON format that can be consumed by analysis scripts.

For the Cryo-EM paper, we should later be able to run:

```python
metadata = load_alignment_metadata(...)
```

and filter:

```python
safe = (
    metadata["alignment"]["missing_selected_residue_count"] == 0
    and metadata["alignment"]["coverage_reference"] > 0.98
)
```

Do not make analysis scripts reverse-engineer PHASE internals.

---

# 50. Scientific interpretation requirement

Do not call an aligned structure "equivalent" to the reference.

The alignment provides:

```text
residue correspondence
```

not:

```text
structural equivalence
```

This distinction must be preserved in UI and documentation.

An aligned active-state GPCR can still have large conformational changes.

---

# 51. Definition of done

The implementation is complete when the following workflow works from both CLI and web UI:

```text
1. Create system
2. Upload reference PDB
3. Select residues on reference
4. Create reference state
5. Attach trajectory
6. Add second trajectory as another state
7. Upload a different PDB
8. Automatic sequence alignment
9. Automatic reference-selection mapping
10. Alignment metadata visible
11. Descriptor generation uses canonical reference order
12. Missing/deleted residues are explicit
13. Insertions do not change canonical ordering
14. Add PDB IDs
15. PDBs download automatically
16. Each PDB is aligned to reference
17. Batch result reports alignment quality
18. Existing legacy projects continue to load
```

The key invariant is:

> **Every state in a system shares the same canonical reference residue coordinate system, regardless of its own PDB numbering.**

This invariant should be enforced throughout the descriptor-building and downstream-analysis pipeline.