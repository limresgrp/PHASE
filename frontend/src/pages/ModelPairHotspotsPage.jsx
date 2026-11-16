import { useCallback, useEffect, useMemo, useState } from 'react';
import { Link, useLocation, useParams } from 'react-router-dom';
import { CircleHelp, Info, Plus, Trash2, X } from 'lucide-react';
import HelpDrawer from '../components/common/HelpDrawer';
import HotspotStructure from '../components/visualization/HotspotStructure';
import { requestJSON } from '../api/client';
import { fetchSystem, fetchClusterAnalyses, fetchClusterAnalysisData, deleteClusterAnalysis, fetchClusterUiSetups } from '../api/projects';
import { fetchJobStatus } from '../api/jobs';
import { hotspotRows, hotspotColor } from '../utils/hotspotMetrics';

const TYPE = 'model_pair_hotspots';
const input = 'rounded border border-gray-700 bg-gray-900 px-2 py-1 text-sm';
const panel = 'rounded-lg border border-gray-700 bg-gray-800 p-4';
const emptyContext = { models: [], samples: [] };
const format = (value) => Number.isFinite(value) ? value.toFixed(4) : 'unavailable';
const colorCss = (value) => `#${value.toString(16).padStart(6, '0')}`;

function MultiSelect({ title, options, value, onChange }) {
  return <label className="block space-y-1 text-sm"><span>{title}</span><select multiple size={Math.min(6, Math.max(2, options.length))} className={`${input} w-full`} value={value}
    onChange={(e) => onChange(Array.from(e.target.selectedOptions, (option) => option.value))}>
    {options.map((item) => <option key={item.sample_id} value={item.sample_id}>{item.name}</option>)}
  </select><span className="block text-xs text-gray-400">Ctrl/Cmd-click to select multiple trajectories.</span></label>;
}

export default function ModelPairHotspotsPage() {
  const { projectId, systemId } = useParams();
  const location = useLocation();
  const [system, setSystem] = useState(null);
  const [chosenCluster, setChosenCluster] = useState(new URLSearchParams(location.search).get('cluster_id') || '');
  const [context, setContext] = useState(emptyContext);
  const [analyses, setAnalyses] = useState([]);
  const [analysisId, setAnalysisId] = useState('');
  const [payload, setPayload] = useState(null);
  const [selections, setSelections] = useState([]);
  const [selectionId, setSelectionId] = useState('');
  const [selectionMode, setSelectionMode] = useState('exclude');
  const [error, setError] = useState('');
  const [loading, setLoading] = useState(false);
  const [help, setHelp] = useState(false);
  const [modal, setModal] = useState(false);
  const [info, setInfo] = useState(null);
  const [job, setJob] = useState(null);
  const [jobStatus, setJobStatus] = useState(null);
  const [submitting, setSubmitting] = useState(false);
  const [revision, setRevision] = useState(0);
  const [form, setForm] = useState({ mode: 'ensemble', model_a_id: '', model_b_id: '', reference_a: [], reference_b: [], sample_ids: [], block_length: 50, bootstrap: 200, seed: 0 });
  const [tab, setTab] = useState('table');
  const [metric, setMetric] = useState('separation');
  const [view, setView] = useState('importance');
  const [sampleIndex, setSampleIndex] = useState(0);
  const [minimum, setMinimum] = useState(0);
  const [top, setTop] = useState(30);
  const [percentile, setPercentile] = useState(0);
  const [reliable, setReliable] = useState(false);
  const [direction, setDirection] = useState('all');
  const [blend, setBlend] = useState(false);
  const [alpha, setAlpha] = useState(.75);
  const [aggregation, setAggregation] = useState('mean');
  const [network, setNetwork] = useState('combined');
  const [focusKey, setFocusKey] = useState('');
  const clusters = useMemo(() => (system?.metastable_clusters || []).filter((c) => c.path), [system]);
  const clusterId = chosenCluster || clusters[0]?.cluster_id || '';
  const base = `/projects/${encodeURIComponent(projectId)}/systems/${encodeURIComponent(systemId)}`;
  const hotspotBase = `${base}/metastable/clusters/${encodeURIComponent(clusterId)}/hotspots`;

  useEffect(() => {
    let active = true;
    fetchSystem(projectId, systemId).then((value) => { if (active) setSystem(value); }).catch((e) => { if (active) setError(e.message); });
    return () => { active = false; };
  }, [projectId, systemId]);

  useEffect(() => {
    setPayload(null); setAnalyses([]); setAnalysisId(''); setContext(emptyContext); setSelectionId(''); setJob(null); setJobStatus(null);
    if (!clusterId) return undefined;
    let active = true;
    Promise.all([requestJSON(`${hotspotBase}/context`), fetchClusterUiSetups(projectId, systemId, clusterId, { setupType: 'residue_selection' })])
      .then(([choices, setups]) => { if (!active) return; setContext(choices); setSelections(setups.setups || []);
        setForm((prev) => ({ ...prev, model_a_id: choices.models[0]?.model_id || '', model_b_id: choices.models[1]?.model_id || '',
          reference_a: choices.models[0]?.suggested_references || [], reference_b: choices.models[1]?.suggested_references || [], sample_ids: [] })); })
      .catch((e) => { if (active) setError(e.message); });
    return () => { active = false; };
  }, [projectId, systemId, clusterId, hotspotBase]);

  useEffect(() => {
    if (!clusterId) return undefined;
    let active = true;
    fetchClusterAnalyses(projectId, systemId, clusterId, { analysisType: TYPE }).then((res) => {
      if (!active) return;
      const list = res.analyses || []; setAnalyses(list);
      setAnalysisId((id) => list.some((a) => a.analysis_id === id) ? id : list[0]?.analysis_id || '');
    }).catch((e) => { if (active) setError(e.message); });
    return () => { active = false; };
  }, [projectId, systemId, clusterId, revision]);

  useEffect(() => {
    setPayload(null); setFocusKey(''); setSampleIndex(0);
    if (!analysisId || !clusterId) return undefined;
    let active = true; setLoading(true);
    fetchClusterAnalysisData(projectId, systemId, clusterId, TYPE, analysisId).then((res) => {
      if (!active) return;
      setPayload(res); setMetric(res.metadata.mode === 'parameter' ? 'parameter_rms' : 'separation');
      setView('importance'); setReliable(false);
    }).catch((e) => { if (active) setError(e.message); }).finally(() => { if (active) setLoading(false); });
    return () => { active = false; };
  }, [projectId, systemId, clusterId, analysisId]);

  useEffect(() => {
    if (!job) return undefined;
    let active = true, timer;
    const poll = async () => {
      try {
        const status = await fetchJobStatus(job); if (!active) return;
        setJobStatus(status);
        if (status.status === 'finished') {
          setRevision((v) => v + 1); setAnalysisId(status.result?.results?.analysis_id || ''); setJob(null); return;
        }
        if (['failed', 'stopped', 'canceled'].includes(status.status)) {
          setError(status.meta?.error || status.result?.error || status.meta?.status || 'Analysis failed.'); setJob(null); return;
        }
      } catch (e) { if (active) setError(e.message); }
      if (active) timer = setTimeout(poll, 2000);
    };
    poll(); return () => { active = false; clearTimeout(timer); };
  }, [job]);

  const labels = useMemo(() => {
    const mapping = Object.fromEntries((system?.reference_residue_records || []).map((r) => [r.residue_key, `${r.resname}${r.resid}${r.insertion_code || ''}`]));
    for (const key of payload?.data?.residue_keys || []) {
      if (!mapping[key]) { const match = key.match(/^res_(-?\d+)_(\w+)$/); mapping[key] = match ? `${match[2]}${match[1]}` : key; }
    }
    return mapping;
  }, [system, payload]);
  const selection = useMemo(() => selections.find((s) => s.setup_id === selectionId)?.payload?.selected_residue_indices || [], [selections, selectionId]);
  const rows = useMemo(() => hotspotRows(payload?.data, { metric, preferenceSample: sampleIndex, view, blend, alpha, aggregation,
    minimum, top, percentile, reliable: reliable && metric === 'separation', direction, selection, selectionMode, labels }),
  [payload, metric, sampleIndex, view, blend, alpha, aggregation, minimum, top, percentile, reliable, direction, selection, selectionMode, labels]);
  const maximum = useMemo(() => [...rows.nodes, ...rows.edges].reduce((v, row) => Math.max(v, Math.abs(row.color)), 0) || 1, [rows]);
  const update = (key, value) => setForm((prev) => ({ ...prev, [key]: value }));
  const changeModel = (side, value) => setForm((prev) => ({ ...prev, [`model_${side}_id`]: value,
    [`reference_${side}`]: context.models.find((m) => m.model_id === value)?.suggested_references || [] }));
  const submit = async (event) => {
    event.preventDefault(); setError(''); setSubmitting(true);
    try {
      const result = await requestJSON(`${hotspotBase}/analyses`, { method: 'POST', body: form });
      setJob(result.job_id); setJobStatus(null); setModal(false);
    } catch (e) { setError(e.message); } finally { setSubmitting(false); }
  };
  const remove = useCallback(async (analysis) => {
    if (!window.confirm(`Delete ${analysis.name || analysis.analysis_id}?`)) return;
    try { await deleteClusterAnalysis(projectId, systemId, clusterId, TYPE, analysis.analysis_id); setRevision((v) => v + 1); }
    catch (e) { setError(e.message); }
  }, [projectId, systemId, clusterId]);
  const referenceOptions = context.samples.filter((s) => s.type === 'md_eval');
  const ensemble = payload?.metadata?.mode === 'ensemble';
  const table = (items, title) => <section className={panel}><h2 className="mb-2 font-semibold">{title}</h2>
    <p className="mb-2 text-xs text-gray-400">Ranked by absolute {metric === 'separation' ? 'endpoint separation' : 'parameter RMS'}. Color represents {view}. Click an identity to focus in 3D.</p>
    <div className="overflow-auto"><table className="w-full text-left text-sm"><thead><tr><th>Identity</th><th>Effect</th><th>95% CI</th><th>Observed A / B</th><th>Displayed</th></tr></thead>
      <tbody>{items.map((row) => <tr className="border-t border-gray-700" key={row.index}><td className="py-2"><button className="text-cyan-400 underline" onClick={() => { setFocusKey(payload.data.residue_keys[row.indices[0]]); setTab('structure'); }}>{row.label}</button></td>
        <td>{format(row.effect)}</td><td>{metric === 'separation' ? `${format(row.low)} / ${format(row.high)}` : '-'}</td><td>{metric === 'separation' ? `${row.countA ?? '-'} / ${row.countB ?? '-'}` : '-'}</td>
        <td><span className="inline-block rounded px-2 py-1 text-black" style={{ background: colorCss(hotspotColor(row.color, maximum, view === 'preference')) }}>{format(row.color)}</span></td></tr>)}</tbody>
    </table>{!items.length && <p className="p-3 text-sm">No terms pass the current filters.</p>}</div></section>;

  return <div className="space-y-4 p-4 text-gray-100">
    <div className="flex flex-wrap items-center justify-between gap-3"><div><h1 className="text-xl font-semibold">Model-pair hotspots</h1><p className="text-sm text-gray-400">Discriminative residues and interactions, not causal mutation predictions.</p></div>
      <div className="flex gap-2"><button className={input} onClick={() => setHelp(true)}><CircleHelp className="mr-1 inline h-4 w-4" />Help</button><button disabled={!clusterId || !!job} className="rounded bg-blue-700 px-3 py-2 font-semibold text-white disabled:opacity-50" onClick={() => setModal(true)}><Plus className="inline h-4 w-4" /> New analysis</button></div></div>
    <div className="flex gap-2 text-sm"><Link className={input} to={`${base}/sampling/delta_eval?cluster_id=${clusterId}`}>Energies, commitment & frustration</Link><span className="rounded bg-cyan-700 px-3 py-1 text-white">Discriminative hotspots</span></div>
    {error && <p role="alert" className="rounded border border-red-500 p-3 text-red-400">{error}</p>}
    {(job || jobStatus) && <section className={panel}><p>{jobStatus?.meta?.status || jobStatus?.status || 'Queued'} {jobStatus?.meta?.progress ?? 0}%</p><progress className="w-full" value={jobStatus?.meta?.progress || 0} max="100" /></section>}
    <div className="grid items-start gap-4 lg:grid-cols-[280px_minmax(0,1fr)]">
      <aside className={`${panel} space-y-4`}><label className="block text-sm">Cluster<select className={`${input} mt-1 w-full`} value={clusterId} onChange={(e) => setChosenCluster(e.target.value)}>{clusters.map((c) => <option key={c.cluster_id} value={c.cluster_id}>{c.name || c.cluster_id}</option>)}</select></label>
        <h2 className="font-semibold">Available analyses</h2>{analyses.map((a) => <div key={a.analysis_id} className={`rounded border p-2 ${analysisId === a.analysis_id ? 'border-cyan-500' : 'border-gray-700'}`}>
          <button className="w-full break-words text-left text-sm" onClick={() => setAnalysisId(a.analysis_id)}>{a.name || a.analysis_id}<span className="block text-xs text-gray-400">{a.mode} · {a.created_at?.slice(0, 10)}</span></button>
          <div className="mt-2 flex justify-end gap-3"><button title="Analysis metadata" onClick={() => setInfo(a)}><Info size={16} /></button><button title="Delete analysis" onClick={() => remove(a)}><Trash2 size={16} /></button></div></div>)}
        {!analyses.length && <p className="text-sm text-gray-400">No hotspot analyses yet.</p>}
        {payload && <><h2 className="font-semibold">Display filters</h2>
          <label className="block text-xs">Ranking metric<select className={`${input} w-full`} value={metric} onChange={(e) => setMetric(e.target.value)}><option value="parameter_rms">Parameter contrast (RMS)</option>{ensemble && <option value="separation">Ensemble separation</option>}</select></label>
          <label className="block text-xs">Color<select className={`${input} w-full`} value={view} onChange={(e) => setView(e.target.value)}><option value="importance">Unsigned importance</option>{ensemble && <option value="preference">Endpoint preference</option>}</select></label>
          {view === 'preference' && <label className="block text-xs">Ensemble<select className={`${input} w-full`} value={sampleIndex} onChange={(e) => setSampleIndex(Number(e.target.value))}>{(payload.metadata.samples || []).map((s, i) => <option key={s.sample_id} value={i}>{s.name}</option>)}</select></label>}
          <label className="block text-xs">Top terms (per table)<input className={`${input} w-full`} type="number" min="1" max="10000" value={top} onChange={(e) => setTop(Math.max(1, Number(e.target.value)))} /></label>
          <label className="block text-xs">Minimum absolute effect<input className={`${input} w-full`} type="number" min="0" step=".01" value={minimum} onChange={(e) => setMinimum(Math.max(0, Number(e.target.value)))} /></label>
          <label className="block text-xs">Keep above percentile<input className={`${input} w-full`} type="number" min="0" max="100" value={percentile} onChange={(e) => setPercentile(Number(e.target.value))} /></label>
          {metric === 'separation' && <><label className="block text-xs"><input type="checkbox" checked={reliable} onChange={(e) => setReliable(e.target.checked)} /> Require individual-term CI excluding zero</label><label className="block text-xs">Separation sign<select className={`${input} w-full`} value={direction} onChange={(e) => setDirection(e.target.value)}><option value="all">Both signs</option><option value="positive">Supports endpoint separation</option><option value="negative">Counter-discriminative</option></select></label></>}
          <label className="block text-xs">Saved residue selection<select className={`${input} w-full`} value={selectionId} onChange={(e) => setSelectionId(e.target.value)}><option value="">Whole protein</option>{selections.map((s) => <option key={s.setup_id} value={s.setup_id}>{s.name}</option>)}</select></label>
          {selectionId && <select className={`${input} w-full`} value={selectionMode} onChange={(e) => setSelectionMode(e.target.value)}><option value="exclude">Exclude selection (e.g. pocket)</option><option value="include">Include selection + incident edges</option></select>}
          <label className="block text-xs"><input type="checkbox" checked={blend} onChange={(e) => setBlend(e.target.checked)} /> Edge-weighted residue color</label>
          {blend && <><label className="block text-xs">Edge weight {alpha.toFixed(2)}<input className="w-full" type="range" min="0" max="1" step=".05" value={alpha} onChange={(e) => setAlpha(Number(e.target.value))} /></label><select className={`${input} w-full`} value={aggregation} onChange={(e) => setAggregation(e.target.value)}><option value="mean">Mean incident edges</option><option value="sum">Sum incident edges (degree-sensitive)</option></select><p className="text-xs text-gray-400">Only displayed colors are mixed, using all incident edges before filtering. Ranking and CIs remain individual-term values, not CIs of the blend.</p></>}
        </>}
      </aside>
      <main className="min-w-0 space-y-4">{loading && <p>Loading compact term summaries...</p>}{!payload && !loading && <section className={panel}>Select an analysis or run a new one.</section>}
        {payload && <><section className={panel}><h2 className="font-semibold">{payload.metadata.name}</h2><p className="mt-1 text-sm text-gray-400">{metric === 'parameter_rms' ? 'RMS parameter differences in a common zero-sum gauge, normalized by alphabet size. Strength has no endpoint sign.' : 'ΔE = E(A) − E(B). Separation = mean(B, ΔE) − mean(A, ΔE); positive terms support endpoint discrimination. Reference trajectories have equal weight.'}</p>
          <p className="mt-2 text-xs">{view === 'preference' ? 'Red: more A-like; blue: more B-like. Colors use minus the ensemble term mean relative to the balanced endpoint midpoint.' : 'Pale to orange: increasing absolute importance; no endpoint direction implied.'} Saturation: {maximum.toPrecision(3)}.</p>
          <div className="mt-3 flex flex-wrap gap-2"><button className={input} onClick={() => setTab('table')}>Ranked terms</button><button className={input} onClick={() => setTab('structure')}>3D network</button>{tab === 'structure' && <select className={input} value={network} onChange={(e) => setNetwork(e.target.value)}><option value="combined">Nodes and edges</option><option value="nodes">Nodes only</option><option value="edges">Edges only</option></select>}</div></section>
          {tab === 'table' ? <>{table(rows.nodes, 'Residue terms')}{table(rows.edges, 'Interaction terms')}</> : <section className={panel}><HotspotStructure projectId={projectId} systemId={systemId} clusterId={clusterId} system={system} data={payload.data} nodes={network === 'edges' ? [] : rows.nodes} edges={network === 'nodes' ? [] : rows.edges} signed={view === 'preference'} maximum={maximum} focusKey={focusKey} /></section>}
        </>}
      </main>
    </div>
    <HelpDrawer open={help} onClose={() => setHelp(false)} title="Model-pair hotspot analysis" docPath="/docs/model_pair_hotspots.md" />
    {info && <div className="fixed inset-0 z-50 flex items-center justify-center bg-black/60 p-4"><section className={`${panel} max-h-[85vh] w-full max-w-3xl overflow-auto`}><button className="float-right" aria-label="Close metadata" onClick={() => setInfo(null)}><X /></button><h2 className="font-semibold">{info.name}</h2><pre className="mt-3 whitespace-pre-wrap break-words text-xs">{JSON.stringify(info, null, 2)}</pre></section></div>}
    {modal && <div className="fixed inset-0 z-50 flex items-center justify-center bg-black/60 p-4"><form className={`${panel} max-h-[90vh] w-full max-w-3xl space-y-4 overflow-y-auto`} onSubmit={submit}>
      <div className="flex justify-between"><h2 className="text-lg font-semibold">New discriminative hotspot analysis</h2><button type="button" aria-label="Close new analysis" onClick={() => setModal(false)}><X /></button></div>
      {error && <p role="alert" className="text-red-400">{error}</p>}
      <p className="text-sm text-gray-400">Full or combined models only; not raw delta patches. Endpoint suggestions must be checked. Select independent held-out MD references where possible.</p>
      <label className="block">Method<select className={`${input} w-full`} value={form.mode} onChange={(e) => update('mode', e.target.value)}><option value="ensemble">Ensemble-supported discrimination + parameter contrast</option><option value="parameter">Parameter contrast only (no references)</option></select></label>
      <div className="grid gap-4 md:grid-cols-2">{['a', 'b'].map((side) => <div className="space-y-3" key={side}><label className="block">Model {side.toUpperCase()}<select className={`${input} w-full`} required value={form[`model_${side}_id`]} onChange={(e) => changeModel(side, e.target.value)}><option value="">Choose model</option>{context.models.map((m) => <option key={m.model_id} value={m.model_id}>{m.name}</option>)}</select></label>
        {form.mode === 'ensemble' && <MultiSelect title={`Endpoint ${side.toUpperCase()} reference MDs`} options={referenceOptions} value={form[`reference_${side}`]} onChange={(value) => update(`reference_${side}`, value)} />}</div>)}</div>
      {form.mode === 'ensemble' && <><MultiSelect title="Additional ensembles for preference coloring (optional)" options={context.samples} value={form.sample_ids} onChange={(value) => update('sample_ids', value)} />
        <div className="grid gap-3 md:grid-cols-3">{[['block_length', 'Block length (frames)', 1], ['bootstrap', 'Bootstrap replicates (0 or ≥50)', 0], ['seed', 'Random seed', 0]].map(([key, label, min]) => <label className="text-sm" key={key}>{label}<input required className={`${input} w-full`} type="number" min={min} max={key === 'bootstrap' ? 5000 : undefined} value={form[key]} onChange={(e) => update(key, Number(e.target.value))} /></label>)}</div>
        <p className="text-xs text-gray-400">Use blocks longer than the correlation times of interest; repeat with different block lengths. Confidence intervals describe fixed-model sampling uncertainty only. Short references may have unavailable intervals.</p></>}
      <button disabled={submitting || form.model_a_id === form.model_b_id || (form.mode === 'ensemble' && (!form.reference_a.length || !form.reference_b.length))} className="rounded bg-blue-700 px-4 py-2 font-semibold text-white disabled:opacity-50">{submitting ? 'Submitting...' : 'Run analysis'}</button>
    </form></div>}
  </div>;
}
