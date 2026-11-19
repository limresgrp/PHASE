import { useCallback, useEffect, useMemo, useState } from 'react';
import { useLocation, useNavigate, useParams } from 'react-router-dom';
import { CircleHelp, Download, Info, Plus, Tags, Trash2, X } from 'lucide-react';
import { requestJSON } from '../api/client';
import { fetchSystem, fetchClusterAnalyses, deleteClusterAnalysis } from '../api/projects';
import { fetchJobStatus } from '../api/jobs';
import HelpDrawer from '../components/common/HelpDrawer';
import ExperimentalBenchmarkResults from '../components/visualization/ExperimentalBenchmarkResults';

const TYPE = 'experimental_structure_benchmark';
const input = 'rounded border border-gray-700 bg-gray-900 px-2 py-1 text-sm';
const panel = 'rounded-lg border border-gray-700 bg-gray-800 p-4';
const empty = { models: [], states: [], md_samples: [], label_schemes: [] };
const blankScheme = { name: '', kind: 'categorical', categoriesText: '', labels: {}, scheme_id: '' };
const defaults = { state_ids: [], state_samples: {}, contrasts: [], scheme_id: '', name: '', assign_missing: true,
  hotspots_mode: 'parameter', top_k: 20, block_length: 50, bootstrap: 200, seed: 0,
  fingerprint_mode: 'standardized', distance: 'euclidean', positive_label: '', negative_label: '', best_oriented: false,
  bh: true, global_threshold: .15, local_threshold: '' };

function Overlay({ title, onClose, children }) {
  return <div role="dialog" aria-modal="true" aria-label={title} className="fixed inset-0 z-50 flex items-center justify-center bg-black/60 p-4"><section className={`${panel} max-h-[90vh] w-full max-w-4xl overflow-auto`}><div className="mb-4 flex items-center justify-between gap-3"><h2 className="text-lg font-semibold">{title}</h2><button aria-label={`Close ${title}`} onClick={onClose}><X /></button></div>{children}</section></div>;
}

export default function ExperimentalBenchmarkPage() {
  const { projectId, systemId, analysisId } = useParams();
  const location = useLocation(), navigate = useNavigate();
  const [system, setSystem] = useState(null);
  const [chosenCluster, setChosenCluster] = useState(new URLSearchParams(location.search).get('cluster_id') || '');
  const [context, setContext] = useState(empty), [analyses, setAnalyses] = useState([]);
  const [result, setResult] = useState(null), [info, setInfo] = useState(null);
  const [error, setError] = useState(''), [loading, setLoading] = useState(false);
  const [help, setHelp] = useState(false), [modal, setModal] = useState(false), [editor, setEditor] = useState(false);
  const [scheme, setScheme] = useState(blankScheme), [form, setForm] = useState(defaults);
  const [busy, setBusy] = useState(false), [revision, setRevision] = useState(0);
  const [job, setJob] = useState(null), [status, setStatus] = useState(null);
  const clusters = useMemo(() => (system?.metastable_clusters || []).filter((c) => c.path), [system]);
  const clusterId = chosenCluster || clusters[0]?.cluster_id || '';
  const page = `/projects/${encodeURIComponent(projectId)}/systems/${encodeURIComponent(systemId)}/experimental_benchmarks`;
  const api = `/projects/${encodeURIComponent(projectId)}/systems/${encodeURIComponent(systemId)}/experimental-benchmarks`;
  const update = (key, value) => setForm((prev) => ({ ...prev, [key]: value }));
  const refresh = useCallback(() => setRevision((v) => v + 1), []);
  useEffect(() => {
    let active = true;
    fetchSystem(projectId, systemId).then((s) => { if (active) setSystem(s); }).catch((e) => { if (active) setError(e.message); });
    return () => { active = false; };
  }, [projectId, systemId]);
  useEffect(() => {
    if (!clusterId) return undefined;
    let active = true;
    Promise.all([requestJSON(`${api}/context?cluster_id=${encodeURIComponent(clusterId)}`), fetchClusterAnalyses(projectId, systemId, clusterId, { analysisType: TYPE })])
      .then(([choices, existing]) => { if (active) { setContext(choices); setAnalyses(existing.analyses || []); } })
      .catch((e) => { if (active) setError(e.message); });
    return () => { active = false; };
  }, [api, projectId, systemId, clusterId, revision]);
  useEffect(() => {
    setResult(null);
    if (!analysisId || !clusterId) return undefined;
    let active = true; setLoading(true);
    requestJSON(`${api}/analyses/${encodeURIComponent(analysisId)}?cluster_id=${encodeURIComponent(clusterId)}`)
      .then((value) => { if (active) setResult(value); }).catch((e) => { if (active) setError(e.message); })
      .finally(() => { if (active) setLoading(false); });
    return () => { active = false; };
  }, [api, analysisId, clusterId]);
  useEffect(() => {
    if (!job) return undefined;
    let active = true, timer;
    const poll = async () => {
      try {
        const state = await fetchJobStatus(job); if (!active) return;
        setStatus(state);
        if (state.status === 'finished') {
          setJob(null); refresh();
          const id = state.result?.results?.analysis_id;
          if (id) navigate(`${page}/${id}?cluster_id=${encodeURIComponent(clusterId)}`);
          return;
        }
        if (['failed', 'canceled', 'stopped'].includes(state.status)) { setError(state.meta?.error || state.result?.error || state.meta?.status || 'Benchmark job failed.'); setJob(null); return; }
      } catch (e) { if (active) setError(e.message); }
      if (active) timer = setTimeout(poll, 2000);
    };
    poll(); return () => { active = false; clearTimeout(timer); };
  }, [job, page, clusterId, navigate, refresh]);
  const openNew = () => {
    setError('');
    setForm({ ...defaults, state_ids: context.states.filter((s) => s.n_frames <= 1).slice(0, 2).map((s) => s.state_id),
      scheme_id: context.label_schemes[0]?.scheme_id || '',
      contrasts: context.models.length > 1 ? [{ model_a_id: context.models[0].model_id, model_b_id: context.models[1].model_id }] : [] });
    setModal(true);
  };
  const submit = async (e) => {
    e.preventDefault(); setBusy(true); setError('');
    try {
      const params = { ...form, cluster_id: clusterId, positive_label: form.positive_label || null, negative_label: form.negative_label || null,
        local_threshold: form.local_threshold === '' ? null : Number(form.local_threshold), name: form.name || null };
      const response = await requestJSON(`${api}/analyses`, { method: 'POST', body: params });
      setJob(response.job_id); setStatus(null); setModal(false);
    } catch (err) { setError(err.message); } finally { setBusy(false); }
  };
  const loadScheme = (id) => {
    const saved = context.label_schemes.find((s) => s.scheme_id === id);
    setScheme(saved ? { ...saved, labels: Object.fromEntries(Object.entries(saved.labels).filter(([id]) => context.states.some((state) => state.state_id === id))), categoriesText: saved.categories.join('\n') } : blankScheme);
  };
  const saveScheme = async (copy = false) => {
    const categories = scheme.categoriesText.split('\n').map((v) => v.trim()).filter(Boolean);
    if (!scheme.name.trim() || !categories.length || new Set(categories).size !== categories.length) {
      setError('Enter a scheme name and unique, nonempty categories.'); return;
    }
    setBusy(true); setError('');
    try {
      const id = copy ? '' : scheme.scheme_id;
      await requestJSON(`${api}/label-schemes${id ? `/${encodeURIComponent(id)}` : ''}`, { method: id ? 'PUT' : 'POST',
        body: { name: scheme.name, kind: scheme.kind, categories: scheme.categoriesText.split('\n').map((v) => v.trim()).filter(Boolean), labels: scheme.labels } });
      setEditor(false); refresh();
    } catch (err) { setError(err.message); } finally { setBusy(false); }
  };
  const deleteScheme = async () => {
    if (!scheme.scheme_id || !window.confirm('Delete this label scheme? Existing benchmark snapshots are preserved.')) return;
    try { await requestJSON(`${api}/label-schemes/${encodeURIComponent(scheme.scheme_id)}`, { method: 'DELETE' }); setScheme(blankScheme); refresh(); }
    catch (err) { setError(err.message); }
  };
  const remove = async (analysis) => {
    if (!window.confirm(`Delete benchmark ${analysis.name}?`)) return;
    try { await deleteClusterAnalysis(projectId, systemId, clusterId, TYPE, analysis.analysis_id); if (analysis.analysis_id === analysisId) navigate(`${page}?cluster_id=${clusterId}`); refresh(); }
    catch (err) { setError(err.message); }
  };
  const changeContrast = (index, key, value) => setForm((prev) => ({ ...prev, contrasts: prev.contrasts.map((c, i) => i === index ? { ...c, [key]: value } : c) }));
  const selectedScheme = context.label_schemes.find((s) => s.scheme_id === form.scheme_id);
  const categories = scheme.categoriesText.split('\n').map((v) => v.trim()).filter(Boolean);
  return <div className="space-y-4 p-4 text-gray-100">
    <div className="flex flex-wrap items-center justify-between gap-3"><div><h1 className="text-xl font-semibold">Experimental Structure Benchmark</h1><p className="text-sm text-gray-400">Compare supplied structures, models and user labels without inferred biological meaning.</p></div><div className="flex flex-wrap gap-2"><button className={input} onClick={() => setHelp(true)}><CircleHelp className="mr-1 inline h-4 w-4" />Help</button><button className={input} onClick={() => { setScheme(blankScheme); setEditor(true); }}><Tags className="mr-1 inline h-4 w-4" />Label schemes</button><button className="rounded bg-blue-700 px-3 py-2 font-semibold text-white disabled:opacity-50" disabled={!clusterId || !!job} onClick={openNew}><Plus className="inline h-4 w-4" /> New benchmark</button></div></div>
    {error && <p role="alert" className="rounded border border-red-500 p-3 text-red-400">{error}</p>}
    {(job || status) && <section className={panel}><p>{status?.meta?.status || status?.status || 'Queued'} ({status?.meta?.progress || 0}%)</p><progress max="100" value={status?.meta?.progress || 0} className="w-full" /></section>}
    <div className="grid items-start gap-4 lg:grid-cols-[260px_minmax(0,1fr)]"><aside className={`${panel} space-y-4`}><label className="block text-sm">Cluster<select className={`${input} mt-1 w-full`} value={clusterId} onChange={(e) => { setChosenCluster(e.target.value); navigate(`${page}?cluster_id=${e.target.value}`); }}><option value="">Choose cluster</option>{clusters.map((c) => <option key={c.cluster_id} value={c.cluster_id}>{c.name || c.cluster_id}</option>)}</select></label><h2 className="font-semibold">Immutable benchmark results</h2>
      {analyses.map((a) => <div key={a.analysis_id} className={`rounded border p-2 ${analysisId === a.analysis_id ? 'border-cyan-500' : 'border-gray-700'}`}><button className="w-full break-words text-left text-sm" onClick={() => navigate(`${page}/${a.analysis_id}?cluster_id=${clusterId}`)}>{a.name}<span className="block text-xs text-gray-400">{a.state_ids?.length} structures · {a.parameters?.contrasts?.length} contrasts<br />{a.created_at?.slice(0, 10)}</span></button><div className="mt-2 flex justify-end gap-3"><button title="Benchmark metadata" onClick={() => setInfo(a)}><Info size={16} /></button><button title="Delete benchmark" onClick={() => remove(a)}><Trash2 size={16} /></button></div></div>)}
      {!analyses.length && <p className="text-sm text-gray-400">No benchmark results yet.</p>}
    </aside><main className="min-w-0 space-y-4">{loading && <p>Loading benchmark summaries...</p>}{result ? <><a className="inline-flex items-center gap-2 rounded bg-blue-700 px-3 py-2 text-sm text-white" href={`/api/v1${api}/analyses/${analysisId}/download?cluster_id=${encodeURIComponent(clusterId)}`}><Download size={16} /> Download CSV / JSON / NPZ bundle</a><ExperimentalBenchmarkResults key={analysisId} result={result} /></> : !loading && <section className={panel}>Select an existing result or create a benchmark. Use PDB-only states with a single assigned categorical frame. Missing residues remain masked. Label schemes are independent system-level analysis annotations.</section>}</main></div>
    <HelpDrawer open={help} onClose={() => setHelp(false)} title="Experimental Structure Benchmark" docPath="/docs/experimental_structure_benchmark.md" />
    {info && <Overlay title="Benchmark metadata" onClose={() => setInfo(null)}><pre className="whitespace-pre-wrap break-words text-xs">{JSON.stringify(info, null, 2)}</pre></Overlay>}
    {editor && <Overlay title="System-level label schemes" onClose={() => setEditor(false)}><div className="space-y-4"><p className="text-sm text-gray-400">Analysis annotations only: no state names or identities are modified. Changing/deleting a scheme does not change saved benchmarks. For ordinal schemes, category order is explicitly the order entered, never inferred from category names.</p>
      {error && <p role="alert" className="text-red-400">{error}</p>}
      <label className="block">Load scheme<select className={`${input} w-full`} value={scheme.scheme_id || ''} onChange={(e) => loadScheme(e.target.value)}><option value="">New scheme</option>{context.label_schemes.map((s) => <option key={s.scheme_id} value={s.scheme_id}>{s.name} ({s.kind}, revision {s.revision})</option>)}</select></label>
      <label className="block">Name<input className={`${input} w-full`} value={scheme.name} onChange={(e) => setScheme((s) => ({ ...s, name: e.target.value }))} /></label><label className="block">Kind<select className={`${input} w-full`} value={scheme.kind} onChange={(e) => setScheme((s) => ({ ...s, kind: e.target.value }))}><option value="categorical">Categorical (no order)</option><option value="ordinal">Ordinal (explicit low-to-high order)</option></select></label>
      <label className="block">Categories, one per line (ordinal: low-to-high)<textarea rows="4" className={`${input} w-full`} value={scheme.categoriesText} onChange={(e) => setScheme((s) => ({ ...s, categoriesText: e.target.value }))} /></label>
      <div className="max-h-80 space-y-2 overflow-auto">{context.states.map((state) => <label className="flex items-center justify-between gap-3 text-sm" key={state.state_id}><span>{state.name}</span><select className={input} value={scheme.labels[state.state_id] || ''} onChange={(e) => setScheme((s) => { const labels = { ...s.labels }; if (e.target.value) labels[state.state_id] = e.target.value; else delete labels[state.state_id]; return { ...s, labels }; })}><option value="">Unlabelled</option>{categories.map((c, i) => <option key={`${i}:${c}`} value={c}>{c}</option>)}</select></label>)}</div>
      <div className="flex flex-wrap gap-2"><button className="rounded bg-blue-700 px-3 py-2 text-white disabled:opacity-50" disabled={busy} onClick={() => saveScheme(false)}>Save scheme</button>{scheme.scheme_id && <><button className={input} disabled={busy} onClick={() => saveScheme(true)}>Save as new scheme</button><button className={input} onClick={deleteScheme}>Delete scheme</button></>}</div></div></Overlay>}
    {modal && <Overlay title="New Experimental Structure Benchmark" onClose={() => setModal(false)}><form onSubmit={submit} className="space-y-4">
      {error && <p role="alert" className="text-red-400">{error}</p>}
      <p className="text-sm text-gray-400">Scores are E_A − E_B in zero-sum gauge. Select static structures, explicit labels and independent model contrasts. Missing coordinates are masked; all comparisons use the same common-support definition.</p>
      <label className="block">Benchmark name<input className={`${input} w-full`} value={form.name} onChange={(e) => update('name', e.target.value)} placeholder="Experimental structure benchmark" /></label>
      <label className="block">Label scheme<select required className={`${input} w-full`} value={form.scheme_id} onChange={(e) => update('scheme_id', e.target.value)}><option value="">Choose an existing scheme</option>{context.label_schemes.map((s) => <option key={s.scheme_id} value={s.scheme_id}>{s.name} ({s.kind})</option>)}</select></label>{!context.label_schemes.length && <p className="text-sm">Create a scheme first using the Label schemes button.</p>}
      <h3 className="font-semibold">Static states (at least two)</h3><div className="max-h-64 space-y-2 overflow-auto">{context.states.map((state) => <div className="flex flex-wrap items-center justify-between gap-2 text-sm" key={state.state_id}><label><input type="checkbox" disabled={state.n_frames > 1} checked={form.state_ids.includes(state.state_id)} onChange={(e) => update('state_ids', e.target.checked ? [...form.state_ids, state.state_id] : form.state_ids.filter((id) => id !== state.state_id))} /> {state.name} <span className="text-xs text-gray-400">({state.n_frames} frames; {selectedScheme?.labels[state.state_id] || 'unlabelled'})</span></label>
        {state.assigned_samples.length > 1 && <select className={input} value={form.state_samples[state.state_id] || ''} onChange={(e) => update('state_samples', { ...form.state_samples, [state.state_id]: e.target.value })}><option value="">Choose one-frame assignment</option>{state.assigned_samples.map((s) => <option key={s.sample_id} value={s.sample_id}>{s.name}</option>)}</select>}
      </div>)}</div><label className="block text-xs"><input type="checkbox" checked={form.assign_missing} onChange={(e) => update('assign_missing', e.target.checked)} /> Assign missing single-frame descriptor-ready states automatically (one core)</label>
      <h3 className="font-semibold">Model contrasts</h3>{form.contrasts.map((contrast, i) => <section className="space-y-2 rounded border border-gray-700 p-3" key={i}><div className="flex flex-wrap items-center gap-2"><span>Axis {i + 1}</span>{['a', 'b'].map((side) => <label key={side}>Model {side.toUpperCase()} <select required className={input} value={contrast[`model_${side}_id`]} onChange={(e) => changeContrast(i, `model_${side}_id`, e.target.value)}>{context.models.map((m) => <option key={m.model_id} value={m.model_id}>{m.name}</option>)}</select></label>)}<button type="button" className={input} onClick={() => update('contrasts', form.contrasts.filter((_, index) => index !== i))}>Remove</button></div>
        {form.hotspots_mode === 'ensemble' && <div className="grid gap-3 md:grid-cols-2">{['a', 'b'].map((side) => <label className="block text-xs" key={side}>Endpoint {side.toUpperCase()} reference MD samples<select multiple size="4" className={`${input} w-full`} value={contrast[`reference_${side}`] || context.models.find((m) => m.model_id === contrast[`model_${side}_id`])?.suggested_references || []} onChange={(e) => changeContrast(i, `reference_${side}`, Array.from(e.target.selectedOptions, (o) => o.value))}>{context.md_samples.map((s) => <option key={s.sample_id} value={s.sample_id}>{s.name}</option>)}</select></label>)}</div>}
      </section>)}<button type="button" className={input} disabled={context.models.length < 2} onClick={() => update('contrasts', [...form.contrasts, { model_a_id: context.models[0].model_id, model_b_id: context.models[1].model_id }])}>Add contrast</button>
      <div className="grid gap-3 md:grid-cols-2"><label>Fingerprint mode<select className={`${input} w-full`} value={form.fingerprint_mode} onChange={(e) => update('fingerprint_mode', e.target.value)}><option value="standardized">Term-wise z-scored (default)</option><option value="raw">Raw contributions</option></select></label><label>PCoA / local distance<select className={`${input} w-full`} value={form.distance} onChange={(e) => update('distance', e.target.value)}><option value="euclidean">Standardized Euclidean</option><option value="correlation">Correlation distance</option></select></label>
        <label>Hotspots<select className={`${input} w-full`} value={form.hotspots_mode} onChange={(e) => update('hotspots_mode', e.target.value)}><option value="parameter">Parameter-only RMS contrasts</option><option value="ensemble">Ensemble-supported discrimination</option></select></label>
        {[['top_k', 'Hotspot top-k', 1], ['global_threshold', 'Global rank RMS threshold', 0], ['seed', 'Seed', 0], ...(form.hotspots_mode === 'ensemble' ? [['block_length', 'Bootstrap block length', 1], ['bootstrap', 'Bootstrap replicates (0 or ≥50)', 0]] : [])].map(([key, label, min]) => <label key={key}>{label}<input type="number" required step={key === 'global_threshold' ? '.01' : '1'} min={min} max={key === 'global_threshold' ? 1 : undefined} className={`${input} w-full`} value={form[key]} onChange={(e) => update(key, Number(e.target.value))} /></label>)}
        <label>Local distance threshold (blank = 75th percentile)<input type="number" min="0" step=".1" className={`${input} w-full`} value={form.local_threshold} onChange={(e) => update('local_threshold', e.target.value)} /></label>
        {['positive_label', 'negative_label'].map((key) => <label key={key}>AUC {key === 'positive_label' ? 'positive' : 'negative'} category<select className={`${input} w-full`} value={form[key]} onChange={(e) => update(key, e.target.value)}><option value="">No AUC</option>{(selectedScheme?.categories || []).map((c) => <option key={c} value={c}>{c}</option>)}</select></label>)}
      </div><label className="block text-xs"><input type="checkbox" checked={form.best_oriented} onChange={(e) => update('best_oriented', e.target.checked)} /> Also report best-oriented AUC separately (preserve original AUC)</label><label className="block text-xs"><input type="checkbox" checked={form.bh} onChange={(e) => update('bh', e.target.checked)} /> Benjamini–Hochberg q for ordinal Spearman tests (raw/common families separately)</label>
      <p className="text-xs text-gray-400">Single PDB references cannot provide reliable bootstrap intervals. Label names have no software-inferred biological ordering or meaning. Constant/unsupported quantities remain unavailable.</p>
      <button disabled={busy || form.state_ids.length < 2 || !form.contrasts.length || !form.scheme_id} className="rounded bg-blue-700 px-4 py-2 font-semibold text-white disabled:opacity-50">{busy ? 'Submitting...' : 'Run benchmark'}</button>
    </form></Overlay>}
  </div>;
}
