import { useState } from 'react';
import { requestJSON } from '../../api/client';
import HelpDrawer from '../common/HelpDrawer';

const field = 'w-full rounded border border-gray-600 bg-gray-900 p-2 text-sm text-gray-100';
const button = 'rounded border border-gray-600 px-3 py-2 text-sm text-gray-100 disabled:opacity-40';
const pct = (value) => value == null ? '—' : `${(100 * value).toFixed(1)}%`;
const highConfidenceImport = (row) => row.status === 'preview' && (['good', 'excellent'].includes(row.alignment?.quality) || (row.alignment?.identity >= 0.95 && row.alignment?.coverage_reference >= 0.95));

export function AlignmentDetails({ alignment, status }) {
  if (!alignment) return <p className="text-xs text-gray-400">Alignment: {status || 'legacy'}{status === 'reference' ? ' · system reference topology' : ' · original selections and shifts retained'}</p>;
  return <details className="rounded border border-gray-700 p-2 text-xs text-gray-300">
    <summary className="cursor-pointer">Alignment: {status || alignment.alignment_status} · {alignment.quality} · identity {pct(alignment.identity)} · coverage {pct(alignment.coverage_reference)}</summary>
    {alignment.error && <p role="alert" className="text-red-400">{alignment.error}</p>}
    <p className="mt-2">Reference chain {alignment.reference_chain || '(blank)'} → target chain {alignment.target_chain || '(blank)'} · target coverage {pct(alignment.coverage_target)}</p>
    <p>Matches {alignment.matches} · mismatches {alignment.mismatches} · inserted residues {alignment.insertion_count} · deleted residues {alignment.deletion_count} · missing selected {alignment.missing_selected_residue_count}</p>
    {alignment.alternate_optimal_mapping && <p className="text-amber-300">Multiple optimal sequence mappings exist. Inspect the correspondence before using this structure.</p>}
    {!!alignment.missing_selected_residues?.length && <p className="text-amber-300">Missing coordinates (masked in descriptors and energies): {alignment.missing_selected_residues.join(', ')}</p>}
    {!!alignment.mismatched_selected_residues?.length && <p className="text-amber-300">Selected sequence substitutions: {alignment.mismatched_selected_residues.join(', ')}. Residue correspondence does not guarantee comparable side-chain chemistry.</p>}
    <div className="max-h-64 overflow-auto"><table className="w-full text-left"><thead><tr><th>Reference residue</th><th>Target residue</th><th>Status</th></tr></thead><tbody>{alignment.correspondence?.map((row) => <tr key={row.reference_residue}><td>{row.reference_residue}</td><td>{row.target_residue || '—'}</td><td>{row.status}</td></tr>)}</tbody></table></div>
    <p className="mt-2">This establishes residue correspondence, not structural equivalence.</p>
  </details>;
}

export default function ReferenceStructurePanel({ projectId, systemId, system, onChanged }) {
  const [modal, setModal] = useState('');
  const [pdb, setPdb] = useState(null);
  const [chain, setChain] = useState('');
  const [selection, setSelection] = useState('protein');
  const [name, setName] = useState('');
  const [source, setSource] = useState('reference');
  const [trajectory, setTrajectory] = useState(null);
  const [slice, setSlice] = useState('');
  const [build, setBuild] = useState(true);
  const [clusterIds, setClusterIds] = useState([]);
  const [assignmentStateIds, setAssignmentStateIds] = useState([]);
  const [assignmentJobs, setAssignmentJobs] = useState([]);
  const [assignmentErrors, setAssignmentErrors] = useState([]);
  const [preview, setPreview] = useState(null);
  const [ids, setIds] = useState('');
  const [batch, setBatch] = useState([]);
  const [selectedIds, setSelectedIds] = useState([]);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState('');
  const [help, setHelp] = useState(false);
  const states = Object.values(system.states || {});
  const clusters = (system.metastable_clusters || []).filter((run) => run.path && !['failed', 'queued', 'running', 'started'].includes(run.status));
  const base = `/projects/${encodeURIComponent(projectId)}/systems/${encodeURIComponent(systemId)}`;
  const open = (kind) => { setModal(kind); setPdb(null); setTrajectory(null); setSlice(''); setSource('reference'); setPreview(null); setError(''); setBatch([]); setChain(''); setSelectedIds([]); setBuild(true); setClusterIds([]); setAssignmentStateIds([]); };
  const selectTrajectory = (file) => {
    if (source === 'reference' && file && /\.(pdb|ent)$/i.test(file.name)) {
      setSource('upload'); setPdb(file); setTrajectory(null); setPreview(null);
    } else {
      setTrajectory(file);
    }
  };
  const queueAssignments = async (stateIds) => {
    if (!stateIds.length || !clusterIds.length) return;
    const jobs = [], errors = [];
    for (const clusterId of clusterIds) {
      try {
        const job = await requestJSON('/submit/md_samples_refresh', {
          method: 'POST', body: { project_id: projectId, system_id: systemId, cluster_id: clusterId,
            state_ids: stateIds, overwrite: true, cleanup: true },
        });
        jobs.push({ ...job, clusterId });
      } catch (err) { errors.push(`${clusterId}: ${err.message}`); }
    }
    setAssignmentJobs(jobs); setAssignmentErrors(errors);
  };
  const perform = async (action) => {
    setBusy(true); setError('');
    try {
      if (modal === 'assign') {
        if (!assignmentStateIds.length || !clusterIds.length) throw new Error('Select states and at least one cluster run.');
        await queueAssignments(assignmentStateIds);
        setModal(''); return;
      }
      if (modal === 'ids') {
        const list = action === 'preview' ? [...new Set(ids.split(/[\s,]+/).filter(Boolean))] : selectedIds;
        const payload = { pdb_ids: list, preview: action === 'preview', chain_ids: Object.fromEntries(list.filter(() => chain.trim()).map((id) => [id, chain.trim()])) };
        if (action !== 'preview') payload.expected_hashes = Object.fromEntries(batch.filter((row) => list.includes(row.pdb_id)).map((row) => [row.pdb_id, row.alignment.target_sha256]));
        payload.build_descriptors = build;
        const result = await requestJSON(`${base}/states/from-pdb-ids`, { method: 'POST', body: payload });
        setBatch(result.results);
        setSelectedIds(result.results.filter(highConfidenceImport).map((row) => row.pdb_id));
        if (action !== 'preview') {
          await onChanged();
          await queueAssignments(result.results.filter((row) => row.status === 'created' && build).map((row) => row.state_id));
        }
        return;
      }
      const body = new FormData();
      if (pdb) body.append('pdb', pdb);
      if (chain.trim()) body.append('chain_id', chain.trim());
      if (modal === 'reference') {
        if (!pdb) throw new Error('Select a reference PDB.');
        body.append('selection', selection);
        await requestJSON(`${base}/reference`, { method: 'POST', body });
      } else {
        if (source === 'upload' && !pdb) throw new Error('Select a target PDB.');
        if (action === 'preview') {
          setPreview(await requestJSON(`${base}/reference/preview`, { method: 'POST', body })); return;
        }
        if (source === 'upload' && !preview) throw new Error('Preview the alignment before creating the state.');
        body.append('name', name.trim());
        if (source === 'reference') body.append('use_reference', 'true');
        else if (source !== 'upload') body.append('source_state_id', source);
        if (trajectory) body.append('trajectory', trajectory);
        if (slice.trim()) body.append('slice_spec', slice.trim());
        body.append('build_descriptors_after_upload', String(build));
        const created = await requestJSON(`${base}/states`, { method: 'POST', body });
        await queueAssignments(Object.values(created.states || {}).filter((state) => !system.states?.[state.state_id] && state.descriptor_file).map((state) => state.state_id));
      }
      await onChanged(); setModal(''); setTrajectory(null); setName('');
    } catch (err) { setError(err.message); }
    finally { setBusy(false); }
  };

  return <section className="space-y-3 rounded-lg border border-gray-700 bg-gray-800 p-4">
    <HelpDrawer open={help} onClose={() => setHelp(false)} title="Reference structure and sequence correspondence" docPath="/docs/reference_alignment.md" />
    <div className="flex flex-wrap items-center justify-between gap-3"><h2 className="text-lg font-semibold">Reference Structure</h2><button className={button} onClick={() => setHelp(true)}>Help</button></div>
    {!!assignmentJobs.length && <div role="status" className="text-sm text-cyan-300">Cluster assignment queued (one compute worker per job). {assignmentJobs.map((job) => <a key={job.job_id} className="mr-3 underline" href={`/jobs/${job.job_id}`}>View assignment progress</a>)}</div>}
    {!!assignmentErrors.length && <div role="alert" className="text-sm text-red-400">Assignment could not be queued: {assignmentErrors.join('; ')}. Uploaded states are retained. Retry using Assign states to clusters.</div>}
    {system.reference_pdb_file && <button className={button} onClick={() => open('assign')}>Assign states to clusters</button>}
    {system.reference_pdb_file ? <>
      <p className="text-sm">Reference PDB: {system.reference_pdb_file} · chain {system.reference_chain_id || '(blank)'} · coordinate sequence length {system.reference_sequence?.length} · selected residues {system.reference_residue_keys?.length}</p>
      <div className="flex flex-wrap gap-2"><a className={button} href={`/api/v1${base}/reference`} download>Download reference</a><button className={button} onClick={() => setModal('map')}>View residue mapping</button><button className={button} disabled={states.length > 0} title="Reference and selection are immutable once states exist" onClick={() => open('reference')}>Change reference / selection</button><button className="rounded bg-blue-700 px-3 py-2 text-sm font-semibold text-white" onClick={() => open('state')}>+ Add state / replica</button><button className={button} onClick={() => open('ids')}>Import PDB IDs</button></div>
    </> : <><p className="text-sm text-gray-400">{states.length ? 'Legacy system: existing selections and residue shifts remain unchanged. Create a new system for reference-based ingestion.' : 'Upload a reference PDB and select its canonical residues before adding states.'}</p><button className={button} disabled={states.length > 0} onClick={() => open('reference')}>Set reference PDB</button></>}
    {modal && <div className="fixed inset-0 z-50 flex items-center justify-center bg-black/70 p-4"><section role="dialog" aria-modal="true" aria-label="Reference structure workflow" className="max-h-[90vh] w-full max-w-3xl space-y-3 overflow-auto rounded-lg border border-gray-600 bg-gray-800 p-5 text-gray-100">
      <div className="flex justify-between"><h3 className="font-semibold">{modal === 'reference' ? 'Set reference structure' : modal === 'ids' ? 'Import PDB IDs' : modal === 'assign' ? 'Assign states to clusters' : modal === 'map' ? 'Canonical reference residue map' : 'Add state / replica'}</h3><button disabled={busy} onClick={() => setModal('')}>Close</button></div>
      {error && <p role="alert" className="text-sm text-red-400">{error}</p>}
      {modal === 'map' ? <div className="overflow-auto"><table className="w-full text-left text-xs"><thead><tr>{['Key', 'Chain', 'PDB residue', 'Residue name', 'Selected'].map((label) => <th key={label}>{label}</th>)}</tr></thead><tbody>{system.reference_residue_records?.map((r) => <tr key={r.residue_key}><td>{r.residue_key}</td><td>{r.chain_id}</td><td>{r.resid}{r.insertion_code}</td><td>{r.resname}</td><td>{system.reference_residue_keys.includes(r.residue_key) ? 'Yes' : 'No'}</td></tr>)}</tbody></table></div> : <>
        {modal === 'state' && <>
          <label className="block text-xs">State name<input className={field} value={name} onChange={(e) => setName(e.target.value)} /></label>
          <label className="block text-xs">Topology<select className={field} value={source} onChange={(e) => { setSource(e.target.value); setPdb(null); setPreview(null); }}><option value="reference">Use reference structure as state</option><option value="upload">Upload another PDB (sequence alignment)</option>{states.map((state) => <option key={state.state_id} value={state.state_id}>New replica of {state.name}</option>)}</select></label>
        </>}
        {(modal === 'reference' || (modal === 'state' && source === 'upload')) && <label className="block text-xs">PDB file<input className={field} type="file" accept=".pdb,.ent" onChange={(e) => { setPdb(e.target.files?.[0]); setPreview(null); }} /></label>}
        {modal === 'state' && source === 'upload' && pdb && <p className="text-xs">Selected PDB: {pdb.name}</p>}
        {(modal === 'reference' || modal === 'ids' || (modal === 'state' && source === 'upload')) && <label className="block text-xs">{modal === 'reference' ? 'Reference chain (blank = longest protein chain)' : 'Target chain (blank = best sequence match; specify to resolve ties)'}<input className={field} value={chain} onChange={(e) => { setChain(e.target.value); setPreview(null); setBatch([]); setSelectedIds([]); }} /></label>}
        {modal === 'reference' && <label className="block text-xs">Residue selection on reference (MDAnalysis)<input className={field} value={selection} onChange={(e) => setSelection(e.target.value)} placeholder="protein or resid 20 to 200" /></label>}
        {modal === 'state' && <>
          <label className="block text-xs">Optional trajectory (must match this topology)<input className={field} type="file" onChange={(e) => selectTrajectory(e.target.files?.[0] || null)} /></label>
          <p className="text-xs text-gray-400">Independent PDBs use Upload another PDB and sequence alignment. A PDB chosen here with the reference topology automatically switches to that workflow; preview its alignment before creating the state. Coordinate-only trajectories must match the selected topology.</p>
          <label className="block text-xs">Frame slice<input className={field} value={slice} onChange={(e) => setSlice(e.target.value)} placeholder="blank = all frames; ::10 = stride 10" /></label>
          {preview && <AlignmentDetails alignment={preview} status={preview.alignment_status} />}
        </>}
        {modal === 'ids' && <>
          <label className="block text-xs">PDB IDs (up to 20; comma or newline separated)<textarea className={field} value={ids} onChange={(e) => { setIds(e.target.value); setBatch([]); setSelectedIds([]); }} /></label>
          {!!batch.length && <><button className={button} onClick={() => setSelectedIds(batch.filter(highConfidenceImport).map((r) => r.pdb_id))}>Select high-confidence mappings (including gaps)</button><div className="space-y-2">{batch.map((row, i) => <div key={`${row.pdb_id}-${i}`} className="rounded border border-gray-600 p-2"><label className="flex gap-2 text-xs"><input type="checkbox" disabled={row.status !== 'preview'} checked={selectedIds.includes(row.pdb_id)} onChange={(e) => setSelectedIds((current) => e.target.checked ? [...current, row.pdb_id] : current.filter((id) => id !== row.pdb_id))} />{row.pdb_id} · {row.status} {row.error || ''}</label>{row.alignment && <AlignmentDetails alignment={row.alignment} status={row.alignment.alignment_status} />}</div>)}</div></>}
        </>}
        {['state', 'ids'].includes(modal) && <label className="flex gap-2 text-xs"><input type="checkbox" checked={build} onChange={(e) => { setBuild(e.target.checked); if (!e.target.checked) setClusterIds([]); }} /> Build descriptors now (PDB-only = one frame)</label>}
        {modal === 'state' && build && !!preview?.missing_selected_residues?.length && <p role="status" className="text-sm text-amber-300">Available residues will be built and assigned. Missing reference residues retain a missing mask; their nodes and incident edges contribute zero energy.</p>}
        {modal === 'assign' && <fieldset className="space-y-2"><legend className="text-sm">Descriptor-ready states</legend>{states.filter((state) => state.descriptor_file).map((state) => <label className="flex gap-2 text-xs" key={state.state_id}><input type="checkbox" checked={assignmentStateIds.includes(state.state_id)} onChange={(e) => setAssignmentStateIds((current) => e.target.checked ? [...current, state.state_id] : current.filter((id) => id !== state.state_id))} />{state.name || state.state_id}</label>)}{!states.some((state) => state.descriptor_file) && <p className="text-xs text-gray-400">Build state descriptors first in the States panel.</p>}</fieldset>}
        {(['state', 'ids'].includes(modal) || modal === 'assign') && <fieldset disabled={busy || (modal !== 'assign' && !build)} className="space-y-2 rounded border border-gray-600 p-3"><legend className="text-sm">{modal === 'assign' ? 'Cluster runs' : 'Also assign to cluster runs (optional)'}</legend>{clusters.map((run) => <label className="flex gap-2 text-xs" key={run.cluster_id}><input type="checkbox" checked={clusterIds.includes(run.cluster_id)} onChange={(e) => setClusterIds((current) => e.target.checked ? [...current, run.cluster_id] : current.filter((id) => id !== run.cluster_id))} />{run.name || run.cluster_name || run.cluster_id}</label>)}{!clusters.length && <p className="text-xs text-gray-400">No completed cluster runs available.</p>}<p className="text-xs text-gray-400">Assignments run as background jobs, one compute worker each. Open the progress link after submission; resulting MD samples appear in Potts Models → Samples.</p></fieldset>}
        <div className="flex gap-2">{(modal === 'ids' || (modal === 'state' && source === 'upload')) && <button className={button} disabled={busy} onClick={() => perform('preview')}>{busy ? 'Working…' : 'Preview alignment(s)'}</button>}<button className="rounded bg-blue-700 px-3 py-2 font-semibold text-white disabled:opacity-40" disabled={busy || (modal === 'state' && (!name.trim() || (source === 'upload' && !preview))) || (modal === 'ids' && !selectedIds.length) || (modal === 'assign' && (!clusterIds.length || !assignmentStateIds.length))} onClick={() => perform('create')}>{busy ? 'Working…' : modal === 'reference' ? 'Save reference' : modal === 'assign' ? 'Run assignment' : 'Create selected state(s)'}</button></div>
        <p className="text-xs text-gray-400">Sequence alignment establishes residue correspondence, not structural equivalence. Missing selected residues are masked, never assigned an ordinary cluster. Their nodes and incident edges contribute zero energy. Trajectories are copied into the shared data root.</p>
      </>}
    </section></div>}
  </section>;
}
