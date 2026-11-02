import { useCallback, useEffect, useMemo, useState } from 'react';
import { useLocation, useNavigate, useParams } from 'react-router-dom';
import { CircleHelp, Play, RefreshCw, Trash2, X } from 'lucide-react';
import Plot from 'react-plotly.js';

import HelpDrawer from '../components/common/HelpDrawer';
import DistanceJsResults from '../components/visualization/DistanceJsResults';
import { fetchSystem, fetchClusterAnalyses, fetchClusterAnalysisData, deleteClusterAnalysis } from '../api/projects';
import { fetchJobStatus } from '../api/jobs';
import {
  fetchDistanceStatus, buildDistanceCache, submitDistanceJs, fetchDistancePairs,
  fetchDistancePair, fetchDistancePresets, saveDistancePreset, deleteDistancePreset,
} from '../api/distances';

const panel = 'rounded-xl border border-gray-700 bg-gray-800 p-4';
const field = 'rounded border border-gray-600 bg-gray-900 px-2 py-2 text-sm text-gray-100';
const button = 'inline-flex items-center justify-center gap-2 rounded border border-gray-600 px-3 py-2 text-xs text-gray-200 hover:border-cyan-500 disabled:opacity-40';
const plotLayout = { paper_bgcolor: '#ffffff', plot_bgcolor: '#ffffff', font: { color: '#1f2937' }, margin: { l: 60, r: 25, t: 30, b: 60 }, hoverlabel: { bgcolor: '#ffffff', font: { color: '#111827' } } };
const fmt = (value, digits = 3) => value == null || !Number.isFinite(Number(value)) ? 'n/a' : Number(value).toFixed(digits);
const split = (text) => text.split(',').map((item) => item.trim()).filter(Boolean);

function thresholdsFrom(text) {
  const values = split(text).map(Number);
  if (values.some((value, i) => !Number.isFinite(value) || value <= 0 || (i && value <= values[i - 1]))) {
    throw new Error('Thresholds must be positive and strictly increasing, for example 6,8,10.');
  }
  return values;
}

function StateSelect({ states, value, onChange, label }) {
  return <label className="block space-y-1 text-xs text-gray-300">
    <span>{label}</span>
    <select multiple value={value} onChange={(event) => onChange(Array.from(event.target.selectedOptions, (option) => option.value))} className={`${field} w-full`} size={Math.min(7, Math.max(2, states.length))}>
      {states.map((state) => <option key={state.state_id} value={state.state_id}>{state.name} · {state.frame_count} frames</option>)}
    </select>
  </label>;
}

export default function DistanceAnalysisPage() {
  const { projectId, systemId } = useParams();
  const location = useLocation();
  const navigate = useNavigate();
  const [system, setSystem] = useState(null);
  const [clusterId, setClusterId] = useState('');
  const [status, setStatus] = useState(null);
  const [error, setError] = useState('');
  const [loading, setLoading] = useState(false);
  const [revision, setRevision] = useState(0);
  const [tab, setTab] = useState('inspect');
  const [modal, setModal] = useState('');
  const [help, setHelp] = useState(false);
  const [states, setStates] = useState([]);
  const [metastable, setMetastable] = useState('');
  const [contactCutoff, setContactCutoff] = useState(6);
  const [thresholdText, setThresholdText] = useState('6,8,10');
  const [filters, setFilters] = useState({});
  const [thresholds, setThresholds] = useState([]);
  const [sort, setSort] = useState('range');
  const [descending, setDescending] = useState(true);
  const [search, setSearch] = useState('');
  const [page, setPage] = useState(0);
  const [table, setTable] = useState({ pairs: [], total: 0 });
  const [tableLoading, setTableLoading] = useState(false);
  const [pair, setPair] = useState(null);
  const [detail, setDetail] = useState(null);
  const [detailLoading, setDetailLoading] = useState(false);
  const [presets, setPresets] = useState([]);
  const [presetId, setPresetId] = useState('');
  const [presetName, setPresetName] = useState('');
  const [cutoff, setCutoff] = useState(12);
  const [periodic, setPeriodic] = useState(true);
  const [overrideText, setOverrideText] = useState('');
  const [conditionA, setConditionA] = useState([]);
  const [conditionB, setConditionB] = useState([]);
  const [metaA, setMetaA] = useState('');
  const [metaB, setMetaB] = useState('');
  const [name, setName] = useState('');
  const [minimum, setMinimum] = useState(1);
  const [job, setJob] = useState(null);
  const [jobStatus, setJobStatus] = useState(null);
  const [analyses, setAnalyses] = useState([]);
  const [analysisId, setAnalysisId] = useState('');
  const [result, setResult] = useState(null);
  const [resultLoading, setResultLoading] = useState(false);
  const [contactOnly, setContactOnly] = useState(false);
  const [minContactFraction, setMinContactFraction] = useState(.1);
  const [submitting, setSubmitting] = useState(false);
  const clusters = useMemo(() => (system?.metastable_clusters || []).filter((item) => item.path), [system]);
  const availableStates = status?.states || [];
  const availableMetas = (status?.metastable_ids || []).filter(Boolean);

  useEffect(() => {
    let cancelled = false;
    fetchSystem(projectId, systemId).then((data) => { if (!cancelled) setSystem(data); }).catch((err) => { if (!cancelled) setError(err.message); });
    return () => { cancelled = true; };
  }, [projectId, systemId]);

  useEffect(() => {
    if (!clusters.length) return;
    const requested = new URLSearchParams(location.search).get('cluster_id');
    setClusterId((current) => clusters.some((item) => item.cluster_id === requested) ? requested
      : clusters.some((item) => item.cluster_id === current) ? current : clusters[0].cluster_id);
  }, [clusters, location.search]);

  useEffect(() => {
    setStatus(null); setPair(null); setDetail(null); setPage(0); setAnalysisId(''); setResult(null); setFilters({});
    setJob(null); setJobStatus(null); setTable({ pairs: [], total: 0 });
    setStates([]); setConditionA([]); setConditionB([]); setMetastable(''); setMetaA(''); setMetaB('');
    setPresetId(''); setPresets([]); setAnalyses([]); setError('');
  }, [clusterId]);

  useEffect(() => {
    if (!clusterId) return;
    let cancelled = false;
    setLoading(true);
    Promise.all([
      fetchDistanceStatus(projectId, systemId, clusterId),
      fetchDistancePresets(projectId, systemId, clusterId),
      fetchClusterAnalyses(projectId, systemId, clusterId, { analysisType: 'distance_js' }),
    ]).then(([data, saved, existing]) => {
      if (cancelled) return;
      setStatus(data); setPresets(saved.presets || []); setAnalyses(existing.analyses || []);
      const ids = (data.states || []).map((state) => state.state_id);
      setFilters((current) => current.state_ids ? current : { state_ids: ids.slice(0, 1).join(','), contact_cutoff: 6 });
      setStates((current) => current.length && current.every((id) => ids.includes(id)) ? current : ids.slice(0, 1));
      setConditionA((current) => current.length && current.every((id) => ids.includes(id)) ? current : ids.slice(0, 1));
      setConditionB((current) => current.length && current.every((id) => ids.includes(id)) ? current : ids.slice(1, 2));
    }).catch((err) => { if (!cancelled) setError(err.message); }).finally(() => { if (!cancelled) setLoading(false); });
    return () => { cancelled = true; };
  }, [projectId, systemId, clusterId, revision]);

  useEffect(() => {
    if (!clusterId || !status) return;
    let cancelled = false;
    setTableLoading(true);
    fetchDistancePairs(projectId, systemId, clusterId, {
      ...filters, sort, descending, search, offset: page * 30, limit: 30,
    }).then((data) => {
      if (cancelled) return;
      setTable(data);
      setPair((current) => data.pairs.some((row) => row.r === current?.[0] && row.s === current?.[1])
        ? current : (data.pairs[0] ? [data.pairs[0].r, data.pairs[0].s] : null));
    }).catch((err) => { if (!cancelled) setError(err.message); }).finally(() => { if (!cancelled) setTableLoading(false); });
    return () => { cancelled = true; };
  }, [projectId, systemId, clusterId, status, filters, sort, descending, search, page]);

  useEffect(() => {
    if (!clusterId || !pair || !status) return;
    let cancelled = false;
    setDetailLoading(true); setDetail(null);
    fetchDistancePair(projectId, systemId, clusterId, pair[0], pair[1], {
      condition: { ...(filters.state_ids ? { state_ids: split(filters.state_ids) } : {}), ...(filters.metastable_ids ? { metastable_ids: split(filters.metastable_ids) } : {}) }, thresholds,
    }).then((data) => { if (!cancelled) setDetail(data); }).catch((err) => { if (!cancelled) setError(err.message); }).finally(() => { if (!cancelled) setDetailLoading(false); });
    return () => { cancelled = true; };
  }, [projectId, systemId, clusterId, pair, filters, thresholds, status]);

  useEffect(() => {
    if (!job?.job_id) return undefined;
    let cancelled = false;
    let timer;
    const poll = async () => {
      try {
        const current = await fetchJobStatus(job.job_id);
        if (cancelled) return;
        setJobStatus(current);
        if (current.status === 'finished') {
          setRevision((value) => value + 1);
          const id = current.result?.results?.metadata?.analysis_id;
          if (id) { setAnalysisId(id); setTab('results'); }
          return;
        }
        if (['failed', 'stopped', 'canceled'].includes(current.status)) {
          setError(current.meta?.error || current.result?.error || current.meta?.status || 'Distance job failed.'); return;
        }
        timer = setTimeout(poll, 1800);
      } catch (err) { if (!cancelled) setError(err.message); }
    };
    poll();
    return () => { cancelled = true; clearTimeout(timer); };
  }, [job]);

  useEffect(() => {
    setResult(null);
    if (!analysisId || !clusterId) return undefined;
    let cancelled = false;
    setResultLoading(true);
    fetchClusterAnalysisData(projectId, systemId, clusterId, 'distance_js', analysisId)
      .then((data) => { if (!cancelled) setResult(data); }).catch((err) => { if (!cancelled) setError(err.message); })
      .finally(() => { if (!cancelled) setResultLoading(false); });
    return () => { cancelled = true; };
  }, [projectId, systemId, clusterId, analysisId]);

  const applyFilters = useCallback(() => {
    try {
      if (!states.length) throw new Error('Select at least one state, or use Select all.');
      if (contactOnly && !status?.distance_available) throw new Error('Build a distance cache before filtering by physical distance.');
      const next = thresholdsFrom(thresholdText);
      if (status?.distance_available && next.some((value) => value > status.metadata.max_distance_cutoff)) throw new Error('Threshold exceeds the candidate cutoff. Extend the distance cache first.');
      setThresholds(next); setPage(0); setError('');
      setFilters({ state_ids: states.join(','), ...(metastable ? { metastable_ids: metastable } : {}), contact_cutoff: contactCutoff, ...(contactOnly ? { min_contact_fraction: Number(minContactFraction) } : {}) });
    } catch (err) { setError(err.message); }
  }, [states, thresholdText, metastable, contactCutoff, contactOnly, minContactFraction, status]);

  const startJob = async () => {
    setError(''); setSubmitting(true);
    try {
      let queued;
      if (modal === 'cache') {
        const overrides = {};
        overrideText.split('\n').map((line) => line.trim()).filter(Boolean).forEach((line) => {
          const at = line.indexOf('=');
          if (at < 1) throw new Error('Trajectory overrides must be state_id=/path/to/trajectory, one per line.');
          overrides[line.slice(0, at).trim()] = line.slice(at + 1).trim();
        });
        queued = await buildDistanceCache(projectId, systemId, clusterId, { max_distance_cutoff: Number(cutoff), periodic, trajectory_overrides: overrides });
      } else {
        if (!conditionA.length || !conditionB.length) throw new Error('Select at least one state for each comparison condition.');
        queued = await submitDistanceJs(projectId, systemId, clusterId, {
          condition_a: { state_ids: conditionA, ...(metaA ? { metastable_ids: [metaA] } : {}) },
          condition_b: { state_ids: conditionB, ...(metaB ? { metastable_ids: [metaB] } : {}) },
          thresholds: thresholdsFrom(thresholdText), min_bin_frames: Number(minimum), name: name || undefined,
        });
      }
      setJobStatus(null); setJob(queued); setModal('');
    } catch (err) { setError(err.message); }
    finally { setSubmitting(false); }
  };

  const presetAction = async (action) => {
    setError('');
    try {
      if (action === 'save') {
        if (!presetName.trim()) throw new Error('Give the preset a name.');
        const saved = await saveDistancePreset(projectId, systemId, clusterId, { name: presetName, thresholds: thresholdsFrom(thresholdText), distance_type: 'ca_ca' });
        setPresetId(saved.setup_id);
      } else if (action === 'load') {
        const saved = presets.find((item) => item.setup_id === presetId);
        if (!saved?.compatible) throw new Error(saved?.compatibility_message || 'Select a compatible preset.');
        setThresholdText(saved.payload.thresholds.join(',')); setThresholds(saved.payload.thresholds);
      } else {
        await deleteDistancePreset(projectId, systemId, clusterId, presetId); setPresetId('');
      }
      const data = await fetchDistancePresets(projectId, systemId, clusterId); setPresets(data.presets || []);
    } catch (err) { setError(err.message); }
  };

  const distribution = useMemo(() => {
    const histogram = detail?.histogram;
    if (!histogram?.edges?.length) return null;
    const x = histogram.edges.slice(0, -1).map((edge, i) => (edge + histogram.edges[i + 1]) / 2);
    const widths = histogram.edges.slice(0, -1).map((edge, i) => histogram.edges[i + 1] - edge);
    return { x, widths, series: histogram.by_state?.length ? histogram.by_state : [{ state_id: 'Selected frames', counts: histogram.counts }] };
  }, [detail]);
  const busy = submitting || (job && !['finished', 'failed', 'stopped', 'canceled'].includes(jobStatus?.status));
  const noPairs = status?.pair_count === 0;
  const metaSelect = (value, change, label) => <label className="block space-y-1 text-xs text-gray-300"><span>{label}</span><select value={value} onChange={(event) => change(event.target.value)} className={`${field} w-full`}><option value="">All metastable memberships</option>{availableMetas.map((id) => <option key={id} value={id}>{id}</option>)}</select></label>;

  return <div className="space-y-4">
    <HelpDrawer open={help} onClose={() => setHelp(false)} title="Distance Analysis: purpose and interpretation" docPath="/docs/distance_analysis.md" />
    <div className="flex flex-wrap items-center justify-between gap-3">
      <div><h1 className="text-2xl font-semibold text-white">Distance Analysis</h1><p className="text-sm text-gray-400">Inspect Cα geometry and compare residue-pair states within the same distance regime.</p></div>
      <div className="flex flex-wrap gap-2">
        <button className={button} disabled={!clusterId || busy} onClick={() => setModal('cache')}><RefreshCw size={15} /> Build / refresh cache</button>
        <button className="inline-flex items-center gap-2 rounded bg-cyan-600 px-3 py-2 text-xs font-semibold text-white disabled:opacity-40" disabled={!status || busy || noPairs} title={noPairs ? 'Build the distance cache first: no residue pairs are available.' : undefined} onClick={() => setModal('analysis')}><Play size={15} /> New conditional-JS analysis</button>
        <button className={button} onClick={() => setHelp(true)}><CircleHelp size={15} /> Help</button>
      </div>
    </div>
    {noPairs && <div role="status" className={`${panel} text-sm`}>No residue pairs are available. Click <strong>Build / refresh cache</strong> to compute candidate pairs from the stored trajectories. If the cache is already built, increase its candidate cutoff and rebuild. Without a cache, the single-bin fallback requires existing contact edges.</div>}
    {error && <div role="alert" className="rounded border border-red-700 bg-red-950 p-3 text-sm text-red-200">{error}<button className="float-right" onClick={() => setError('')} aria-label="Dismiss error"><X size={16} /></button></div>}
    {job && <div className={panel}><p className="mb-2 text-xs text-gray-300">{jobStatus?.meta?.status || jobStatus?.status || 'Queued'} · {jobStatus?.meta?.progress || 0}%</p><div className="h-2 rounded bg-gray-700"><div className="h-2 rounded bg-cyan-500" style={{ width: `${jobStatus?.meta?.progress || 0}%` }} /></div></div>}
    <div className="grid gap-4 xl:grid-cols-[300px_minmax(0,1fr)]">
      <aside className="min-w-0 space-y-4">
        <section className={`${panel} space-y-3`}>
          <label className="block text-xs text-gray-300">Cluster<select className={`${field} mt-1 w-full`} value={clusterId} onChange={(event) => {
            const id = event.target.value;
            setClusterId(id);
            navigate(`${location.pathname}${id ? `?cluster_id=${encodeURIComponent(id)}` : ''}`, { replace: true });
          }}><option value="">Select cluster</option>{clusters.map((item) => <option key={item.cluster_id} value={item.cluster_id}>{item.name || item.cluster_id}</option>)}</select></label>
          {loading ? <p className="text-xs text-gray-400">Loading distance metadata…</p> : <p className="text-xs text-gray-400">{status?.distance_available ? `${status.metadata.frame_count} frames · ${status.metadata.edge_count} pairs · candidate cutoff ${status.metadata.max_distance_cutoff} Å` : 'No distance cache: each pair uses one All distances bin. Build a cache to enable distance conditioning.'}</p>}
          <h2 className="text-sm font-semibold">Geometry inspection filters</h2>
          <p className="text-[11px] text-gray-400">These settings affect Pair distances &amp; states only. Saved JS results have their own bin/support filters.</p>
          <StateSelect states={availableStates} value={states} onChange={setStates} label="Inspection trajectories (Ctrl/Cmd to select several)" />
          <button className={button} onClick={() => setStates(availableStates.map((state) => state.state_id))}>Select all trajectories</button>
          {metaSelect(metastable, setMetastable, 'Inspection metastable filter')}
          <label className="block text-xs text-gray-300">Contact distance cutoff (Å)<input type="number" min="0.1" step="0.5" className={`${field} mt-1 w-full`} value={contactCutoff} onChange={(event) => setContactCutoff(event.target.value)} /></label>
          <label className="flex items-center gap-2 text-xs text-gray-300"><input type="checkbox" checked={contactOnly} disabled={!status?.distance_available} onChange={(event) => setContactOnly(event.target.checked)} /> Only pairs contacting below the cutoff</label>
          {contactOnly && <label className="block text-xs text-gray-300">Minimum fraction of selected frames<input className={`${field} ml-2 w-20`} type="number" min="0" max="1" step=".05" value={minContactFraction} onChange={(event) => setMinContactFraction(event.target.value)} /></label>}
          <p className="text-[11px] text-gray-400">Contact means Cα distance strictly below the cutoff. 0 = at least one frame; 0.1 = at least 10%; 1 = every frame. Apply with Refresh inspection. This filters pairs, not their frames.</p>
          <label className="block text-xs text-gray-300">Distance thresholds (Å)<input className={`${field} mt-1 w-full`} placeholder="6,8,10 (blank = global)" value={thresholdText} onChange={(event) => setThresholdText(event.target.value)} /></label>
          <button className={`${button} w-full`} disabled={!status} onClick={applyFilters}><RefreshCw size={14} /> Refresh inspection</button>
          <p className="text-[11px] text-gray-400">Display thresholds split saved frames; they do not change residue cluster definitions or the candidate cutoff.</p>
        </section>
        <section className={`${panel} space-y-2`}>
          <h2 className="text-sm font-semibold text-white">Saved distance presets</h2>
          <select value={presetId} onChange={(event) => setPresetId(event.target.value)} className={`${field} w-full`}><option value="">Select preset</option>{presets.map((item) => <option key={item.setup_id} value={item.setup_id}>{item.name}{item.compatible ? '' : ' (extend cache)'}</option>)}</select>
          <div className="flex gap-2"><button className={button} disabled={!presetId} onClick={() => presetAction('load')}>Load</button><button className={button} disabled={!presetId} onClick={() => presetAction('delete')}><Trash2 size={13} /> Delete</button></div>
          <input value={presetName} onChange={(event) => setPresetName(event.target.value)} placeholder="New preset name" className={`${field} w-full`} />
          <button className={button} disabled={!status} onClick={() => presetAction('save')}>Save current thresholds</button>
        </section>
        <section className={`${panel} space-y-2`}>
          <h2 className="text-sm font-semibold text-white">Conditional-JS analyses</h2>
          {!analyses.length && <p className="text-xs text-gray-400">No analyses yet.</p>}
          {analyses.map((item) => <div key={item.analysis_id} className={`flex gap-2 rounded border p-2 ${analysisId === item.analysis_id ? 'border-cyan-500' : 'border-gray-700'}`}><button className="min-w-0 flex-1 text-left text-xs text-gray-200" onClick={() => { setAnalysisId(item.analysis_id); setTab('results'); }}><span className="block break-words">{item.name || item.analysis_id}</span><span className="text-gray-400">{item.summary?.n_bins} bins · {item.summary?.n_edges} pairs</span></button><button className="text-gray-400 hover:text-red-400" title="Delete analysis" onClick={async () => {
            if (!window.confirm(`Delete ${item.name || 'this analysis'}?`)) return;
            try { await deleteClusterAnalysis(projectId, systemId, clusterId, 'distance_js', item.analysis_id); if (analysisId === item.analysis_id) setAnalysisId(''); setRevision((value) => value + 1); } catch (err) { setError(err.message); }
          }}><Trash2 size={14} /></button></div>)}
        </section>
      </aside>
      <main className="min-w-0 space-y-4">
        <div className="flex gap-2">{[['inspect', 'Pair distances & states'], ['results', 'Conditional-JS results']].map(([id, label]) => <button key={id} onClick={() => setTab(id)} className={`${button} ${tab === id ? 'bg-cyan-700 text-white' : ''}`}>{label}</button>)}</div>
        {tab === 'inspect' ? <>
          <section className={`${panel} space-y-3`}>
            <div className="flex flex-wrap items-center justify-between gap-2"><h2 className="text-base font-semibold text-white">Candidate pair geometry</h2><span className="text-xs text-gray-400">{table.total} pairs · {tableLoading ? 'Loading…' : 'click a row to inspect'}</span></div>
            <div className="flex flex-wrap gap-2"><input className={field} placeholder="Find residue key" value={search} onChange={(event) => { setSearch(event.target.value); setPage(0); }} /><select className={field} value={sort} onChange={(event) => { setSort(event.target.value); setPage(0); }}>{[['range', 'Distance range'], ['contact_occupancy', 'Contact occupancy'], ['min', 'Minimum distance'], ['median', 'Median distance'], ['q10', 'Q10'], ['q90', 'Q90'], ['frame_count', 'Frame count'], ['state_count', 'State count'], ['residue_a', 'Residue A'], ['residue_b', 'Residue B']].map(([id, label]) => <option key={id} value={id}>{label}</option>)}</select><button className={button} onClick={() => setDescending((value) => !value)}>{descending ? 'Highest first ↓' : 'Lowest first ↑'}</button></div>
            <div className="overflow-x-auto"><table className="w-full text-left text-xs text-gray-300"><thead className="border-b border-gray-600 text-gray-400"><tr>{['Residue A', 'Residue B', 'Contact fraction', 'Min Å', 'Median Å', 'Q10 Å', 'Q90 Å', 'Range Å', 'Frames', 'States'].map((label) => <th key={label} className="whitespace-nowrap p-2">{label}</th>)}</tr></thead><tbody>{table.pairs.map((row) => <tr key={`${row.r}-${row.s}`} tabIndex={0} onKeyDown={(event) => { if (event.key === 'Enter') setPair([row.r, row.s]); }} onClick={() => setPair([row.r, row.s])} className={`cursor-pointer border-b border-gray-700 hover:bg-gray-700 ${pair?.[0] === row.r && pair?.[1] === row.s ? 'bg-gray-700' : ''}`}><td className="p-2">{row.residue_a}</td><td className="p-2">{row.residue_b}</td>{['contact_occupancy', 'min', 'median', 'q10', 'q90', 'range'].map((key) => <td className="p-2" key={key}>{fmt(row[key])}</td>)}<td className="p-2">{row.frame_count}</td><td className="p-2">{row.state_count}</td></tr>)}</tbody></table></div>
            <div className="flex items-center gap-3 text-xs text-gray-400"><button className={button} disabled={!page} onClick={() => setPage((value) => value - 1)}>Previous</button>Page {page + 1}<button className={button} disabled={(page + 1) * 30 >= table.total} onClick={() => setPage((value) => value + 1)}>Next</button></div>
          </section>
          <section className={`${panel} space-y-3`}><h2 className="text-base font-semibold text-white">{detail?.residue_keys?.join(' — ') || 'Select a pair'} · distance distribution</h2>
            {detailLoading ? <p className="text-sm text-gray-400">Loading aligned pair data…</p> : distribution ? <Plot data={distribution.series.map((series) => ({ x: distribution.x, y: series.counts, width: distribution.widths, type: 'bar', name: availableStates.find((state) => state.state_id === series.state_id)?.name || series.state_id, opacity: .55, hovertemplate: '%{x:.3f} Å<br>%{y} frames<extra>%{fullData.name}</extra>' }))} layout={{ ...plotLayout, height: 330, barmode: 'overlay', xaxis: { title: 'Cα distance (Å)' }, yaxis: { title: 'Frame count' }, legend: { orientation: 'h', y: -.25 }, shapes: thresholds.map((value) => ({ type: 'line', x0: value, x1: value, y0: 0, y1: 1, yref: 'paper', line: { color: '#475569', dash: 'dash' } })) }} useResizeHandler config={{ responsive: true, displayModeBar: false }} style={{ width: '100%' }} /> : <p className="text-sm text-gray-400">Physical distance distribution unavailable without a cache, or no frames match this filter.</p>}
            <div className="flex flex-wrap gap-3 text-xs text-gray-300">{detail?.bin_labels?.map((label, i) => <span key={label}>{label}: {detail.bin_counts[i]} frames ({fmt(100 * detail.bin_counts[i] / Math.max(1, detail.frame_count), 1)}%)</span>)}</div>
          </section>
          <section className={`${panel} space-y-3`}><h2 className="text-base font-semibold text-white">Pair-state distributions conditioned on distance</h2><p className="text-xs text-gray-400">Each heatmap is normalized within its own frame subset. Residue cluster IDs stay fixed. Invalid labels are excluded for this pair only; empty subsets are shown explicitly.</p>
            <div className="grid gap-3 2xl:grid-cols-2">{detail?.panels?.map((item) => <div key={item.label} className="min-w-0 rounded border border-gray-600 bg-white p-2"><p className="text-center text-xs font-semibold text-gray-800">{item.label} · {item.count} pair-valid frames</p>{item.count ? <Plot data={[{ type: 'heatmap', z: item.probabilities, x: detail.cluster_labels_s.map((value) => `c${value}`), y: detail.cluster_labels_r.map((value) => `c${value}`), colorscale: 'Blues', zmin: 0, zmax: 1, colorbar: { title: 'p' }, hovertemplate: '%{y}, %{x}<br>p = %{z:.4f}<extra></extra>' }]} layout={{ ...plotLayout, height: 300, xaxis: { title: `${detail.residue_keys[1]} cluster` }, yaxis: { title: `${detail.residue_keys[0]} cluster` } }} useResizeHandler config={{ displayModeBar: false, responsive: true }} style={{ width: '100%' }} /> : <p className="py-20 text-center text-xs text-gray-600">No valid frames in this subset</p>}</div>)}</div>
          </section>
        </> : <section className={`${panel} space-y-4`}>
          <h2 className="text-base font-semibold text-white">{result?.metadata?.name || 'Select a saved conditional-JS analysis'}</h2>
          {resultLoading && <p className="text-sm text-gray-400">Loading results…</p>}
          {result && <DistanceJsResults key={analysisId} result={result} projectId={projectId} systemId={systemId} clusterId={clusterId} system={system} residueNames={status?.residue_names} />}
        </section>}
      </main>
    </div>
    {modal && <div className="fixed inset-0 z-40 flex items-center justify-center bg-black/70 p-4"><section role="dialog" aria-modal="true" aria-label={modal === 'cache' ? 'Build distance cache' : 'New conditional-JS analysis'} className={`${panel} max-h-[90vh] w-full max-w-2xl space-y-4 overflow-y-auto shadow-xl`}>
      <div className="flex items-center justify-between"><h2 className="font-semibold text-white">{modal === 'cache' ? 'Build / refresh Cα distance cache' : 'New conditional-JS analysis'}</h2><button className={button} onClick={() => setModal('')}><X size={15} /> Close</button></div>
      {error && <p role="alert" className="rounded border border-red-700 bg-red-950 p-2 text-sm text-red-200">{error}</p>}
      {modal === 'cache' ? <>
        <p className="text-sm text-gray-400">Uses every merged cluster frame with saved state slices, strides and residue mappings. This background job can run with one worker.</p>
        <label className="block text-xs text-gray-300">Candidate-edge cutoff (Å)<input type="number" min=".1" step=".5" value={cutoff} onChange={(event) => setCutoff(event.target.value)} className={`${field} ml-3`} /></label>
        <label className="flex items-center gap-2 text-xs text-gray-300"><input type="checkbox" checked={periodic} onChange={(event) => setPeriodic(event.target.checked)} /> Minimum-image distance when a valid periodic box is available</label>
        <label className="block space-y-1 text-xs text-gray-300"><span>Optional raw trajectory overrides (state_id=/path, one per line)</span><textarea className={`${field} w-full`} rows={3} value={overrideText} onChange={(event) => setOverrideText(event.target.value)} /></label>
      </> : <>
        <input className={`${field} w-full`} placeholder="Analysis name (optional)" value={name} onChange={(event) => setName(event.target.value)} />
        <div className="grid gap-4 sm:grid-cols-2"><div className="space-y-2"><StateSelect states={availableStates} value={conditionA} onChange={setConditionA} label="Condition A trajectories" />{metaSelect(metaA, setMetaA, 'Condition A metastable subset')}</div><div className="space-y-2"><StateSelect states={availableStates} value={conditionB} onChange={setConditionB} label="Condition B trajectories" />{metaSelect(metaB, setMetaB, 'Condition B metastable subset')}</div></div>
        <label className="block text-xs text-gray-300">Distance thresholds (Å)<input value={thresholdText} placeholder="6,8,10 (blank = global)" onChange={(event) => setThresholdText(event.target.value)} className={`${field} ml-2`} /></label>
        <label className="block text-xs text-gray-300">Minimum valid frames per condition per bin<input type="number" min="1" value={minimum} onChange={(event) => setMinimum(event.target.value)} className={`${field} ml-2 w-24`} /></label>
        <p className="text-xs text-gray-400">Uses assigned labels. Frames with invalid labels for this pair are excluded. JS is saved independently in nats; occupancy and valid counts accompany every bin. {!status?.distance_available && 'No cache exists, so requested bins collapse to All distances.'}</p>
      </>}
      <button className="rounded bg-cyan-600 px-4 py-2 text-sm font-semibold text-white disabled:opacity-40" disabled={submitting} onClick={startJob}>{submitting ? 'Submitting…' : 'Run background job'}</button>
    </section></div>}
  </div>;
}
