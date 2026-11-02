import { useCallback, useMemo, useState } from 'react';
import Plot from 'react-plotly.js';
import DistanceJsStructure from './DistanceJsStructure';
import { binsBelowCutoff, distanceJsRows, distanceJsResidues, JS_MAX } from '../../utils/distanceJs';

const input = 'rounded border border-gray-600 bg-gray-900 p-2 text-xs text-gray-100';
const layout = { paper_bgcolor: '#fff', plot_bgcolor: '#fff', font: { color: '#1f2937' }, hoverlabel: { bgcolor: '#fff', font: { color: '#111827' } }, margin: { l: 65, r: 30, t: 35, b: 90 } };
const fmt = (value) => value == null ? 'n/a' : Number(value).toFixed(3);
const emptyNames = {};

export default function DistanceJsResults({ result, projectId, systemId, clusterId, system, residueNames = emptyNames }) {
  const { data, metadata } = result;
  const [normalized, setNormalized] = useState(true);
  const [sort, setSort] = useState(data.distance_bin_labels.length > 1 ? 'change' : 'global');
  const [minimum, setMinimum] = useState(20);
  const [minChange, setMinChange] = useState(0);
  const [cutoff, setCutoff] = useState('');
  const [search, setSearch] = useState('');
  const [page, setPage] = useState(0);
  const [selected, setSelected] = useState(null);
  const [view, setView] = useState('ranking');
  const residueLabel = useCallback((key) => residueNames[key] ? `${residueNames[key]}${key.match(/-?\d+/)?.[0] || key}` : key, [residueNames]);
  const pairLabel = useCallback((i) => [data.residue_keys[data.edge_index[0][i]], data.residue_keys[data.edge_index[1][i]]].map(residueLabel).join(' — '), [data, residueLabel]);
  const bins = useMemo(() => binsBelowCutoff(metadata, data.distance_bin_labels.length, cutoff), [metadata, data, cutoff]);
  const rows = useMemo(() => distanceJsRows(data, { bins, minFrames: Math.max(1, Number(minimum) || 1), minChange: Number(minChange), sort })
    .filter((row) => `${data.edge_labels[row.index]} ${pairLabel(row.index)}`.toLowerCase().includes(search.toLowerCase())), [data, bins, minimum, minChange, sort, search, pairLabel]);
  const residues = useMemo(() => distanceJsResidues(rows, data), [rows, data]);
  const visible = rows.slice(page * 30, (page + 1) * 30);
  const current = rows.some((r) => r.index === selected) ? selected : rows[0]?.index;
  const display = (raw) => raw == null ? null : raw / (normalized ? JS_MAX : 1);
  const unit = normalized ? 'normalized JS (0–1)' : 'JS (nats)';
  const set = (setter) => (event) => { setter(event.target.value); setPage(0); };
  const conditionName = (condition) => condition?.label || (condition?.state_ids || []).map((id) => system?.states?.[id]?.name || id).join(', ') || 'All frames';

  return <div className="space-y-4 text-gray-200">
    <p className="text-sm"><strong>A:</strong> {conditionName(metadata.condition_a)} · <strong>B:</strong> {conditionName(metadata.condition_b)}</p>
    <div className="rounded border border-cyan-700 p-3 text-sm space-y-2">
      <p><strong>What changes with distance?</strong> Each value compares the joint residue-cluster distributions of A and B, using only frames in the same distance bin. High JS means the ensembles differ even at comparable distances.</p>
      <p><strong>JS range = maximum − minimum bin JS.</strong> A large range means the A/B difference depends on distance. Global JS alone cannot tell you this. Only bins with enough valid frames in both ensembles are compared.</p>
      <p className="text-xs text-gray-400">This is an effect-size screen, not a significance test or evidence of causation/allostery. Sparse bins and correlated MD frames require caution; validate candidates with independent trajectories or block-resampled uncertainty.</p>
    </div>
    <div className="flex flex-wrap gap-3 text-xs">
      <label>Rank pairs by <select className={input} value={sort} onChange={set(setSort)}><option value="change">JS range across bins</option><option value="global">Global JS</option>{data.distance_bin_labels.map((label, i) => <option key={label} value={i}>{label} JS</option>)}</select></label>
      <label>Minimum nA and nB / bin <input className={`${input} w-20`} type="number" min="1" value={minimum} onChange={set(setMinimum)} /></label>
      <label>Minimum JS range (0–1) <input className={`${input} w-20`} type="number" min="0" max="1" step=".05" value={minChange} onChange={set(setMinChange)} /></label>
      <label>Use distances below <select className={input} value={cutoff} onChange={set(setCutoff)}><option value="">All saved bins</option>{(metadata.thresholds || []).map((value) => <option key={value} value={value}>{value} Å</option>)}</select></label>
      <label>Find pair <input className={input} value={search} onChange={set(setSearch)} placeholder="residue key" /></label>
      <label className="flex items-center gap-2"><input type="checkbox" checked={normalized} onChange={(event) => setNormalized(event.target.checked)} /> Show normalized JS</label>
    </div>
    <p className="text-xs text-gray-400">{rows.length} supported pairs · {bins.length} selected distance bins. Filters use whole saved bins, never approximated new boundaries. Range filtering always uses 0–1 units; raw display only changes JS labels. A JS range requires at least two supported bins.</p>
    {!metadata.distance_available && <p className="text-sm">No geometry was available for this analysis: only global/unconditional JS can be shown. Build a cache and rerun with at least one threshold to measure distance dependence.</p>}
    {metadata.distance_available && data.distance_bin_labels.length < 2 && <p className="text-sm text-amber-300">This analysis used no distance thresholds and contains only global JS. Run a new analysis with thresholds such as 6,8,10 Å to measure changes with distance; existing results cannot reconstruct missing bins.</p>}
    <div className="flex gap-2">{[['ranking', 'Pair & residue rankings'], ['structure', '3D structure']].map(([id, label]) => <button key={id} className={`${input} ${view === id ? 'border-cyan-500' : ''}`} onClick={() => setView(id)}>{label}</button>)}</div>
    {!rows.length && <p role="status" className="text-sm">No pairs pass these settings. Lower the minimum counts/range, include more distance bins, or rank by global JS. Bins saved as n/a require rerunning the analysis with suitable thresholds/counts.</p>}
    {view === 'structure' ? <DistanceJsStructure projectId={projectId} systemId={systemId} clusterId={clusterId} system={system} data={data} rows={rows} /> : <>
      <h3 className="font-semibold">Distance-bin JS heatmap · click a pair below for its profile</h3>
      {!!visible.length && <Plot data={[{ type: 'heatmap', x: visible.map((row) => pairLabel(row.index)), y: bins.map((b) => data.distance_bin_labels[b]),
        z: bins.map((b) => visible.map(({ index, supported }) => supported.includes(b) ? display(data.edge_bin_js[b][index]) : null)),
        zmin: 0, zmax: normalized ? 1 : JS_MAX, colorscale: 'YlOrRd', colorbar: { title: unit }, hoverongaps: false,
        hovertemplate: '%{x}<br>%{y}<br>JS = %{z:.3f}<extra></extra>' }]} layout={{ ...layout, height: 380, margin: { l: 130, r: 30, t: 30, b: 140 }, xaxis: { title: 'Residue pair', tickangle: -60 }, yaxis: { title: 'Distance bin' } }} config={{ responsive: true }} useResizeHandler style={{ width: '100%' }} />}
      <div className="overflow-x-auto"><table className="w-full text-left text-xs"><thead><tr><th className="p-2">Pair (select)</th><th className="p-2">JS range {normalized ? '(0–1)' : '(nats)'}</th><th className="p-2">Global JS</th>{bins.map((b) => <th key={b} className="p-2">{data.distance_bin_labels[b]}: JS · nA/nB · occupancy A/B</th>)}</tr></thead><tbody>{visible.map((row) => <tr className={`border-t border-gray-700 ${current === row.index ? 'bg-gray-700' : ''}`} key={row.index}>
        <td className="p-2"><button className="text-cyan-400 underline" onClick={() => setSelected(row.index)}>{pairLabel(row.index)}</button></td>
        <td className="p-2">{fmt(row.change == null ? null : row.change * (normalized ? 1 : JS_MAX))}</td><td className="p-2">{fmt(display(row.global))}</td>{bins.map((b) => <td key={b} className="whitespace-nowrap p-2">{row.supported.includes(b) ? fmt(display(data.edge_bin_js[b][row.index])) : 'n/a'} · {data.edge_bin_counts_a[b][row.index]}/{data.edge_bin_counts_b[b][row.index]} · {fmt(data.edge_bin_occupancy_a[b][row.index])}/{fmt(data.edge_bin_occupancy_b[b][row.index])}</td>)}</tr>)}</tbody></table></div>
      <div className="flex items-center gap-3 text-xs"><button className={input} disabled={!page} onClick={() => setPage((v) => v - 1)}>Previous</button>Page {page + 1}<button className={input} disabled={(page + 1) * 30 >= rows.length} onClick={() => setPage((v) => v + 1)}>Next</button></div>
      {current != null && <>
        <h3 className="font-semibold">{pairLabel(current)} · JS versus distance</h3>
        <Plot data={[{ type: 'scatter', mode: 'lines+markers', x: bins.map((b) => data.distance_bin_labels[b]),
          y: bins.map((b) => rows.find((r) => r.index === current)?.supported.includes(b) ? display(data.edge_bin_js[b][current]) : null), connectgaps: false,
          customdata: bins.map((b) => [data.edge_bin_counts_a[b][current], data.edge_bin_counts_b[b][current]]), name: 'Conditional JS',
          hovertemplate: '%{x}<br>JS %{y:.3f}<br>nA %{customdata[0]}, nB %{customdata[1]}<extra></extra>', line: { color: '#0891b2' } }]} layout={{ ...layout, height: 310, xaxis: { title: 'Distance bin (categorical)' }, yaxis: { title: unit, range: [0, normalized ? 1 : JS_MAX] }, shapes: rowGlobalLine(display(data.edge_global_js[current])) }} config={{ responsive: true }} useResizeHandler style={{ width: '100%' }} />
        <p className="text-xs text-gray-400">Dashed horizontal line: unconditional/global JS. Missing bins remain gaps, not zero. Conditional versus global JS also reflects changing bin populations; prioritize the between-bin range.</p>
      </>}
      <h3 className="font-semibold">Residues participating in distance-dependent pairs</h3>
      <p className="text-xs text-gray-400">Ranked by mean JS range across incident pairs passing the filters; this is an edge-derived residue score, not single-residue JS. Partner count and maximum are shown to distinguish broad effects from one strong partner.</p>
      <div className="max-h-72 overflow-auto"><table className="w-full text-left text-xs"><thead><tr>{['Residue', 'Mean range (0–1)', 'Maximum range (0–1)', 'Partners'].map((label) => <th key={label} className="p-2">{label}</th>)}</tr></thead><tbody>{residues.map((row) => <tr key={row.key}><td className="p-2">{residueLabel(row.key)}</td><td className="p-2">{fmt(row.mean)}</td><td className="p-2">{fmt(row.max)}</td><td className="p-2">{row.count}</td></tr>)}</tbody></table></div>
    </>}
  </div>;
}

function rowGlobalLine(value) {
  return value == null ? [] : [{ type: 'line', xref: 'paper', x0: 0, x1: 1, y0: value, y1: value, line: { dash: 'dash', color: '#64748b' } }];
}
