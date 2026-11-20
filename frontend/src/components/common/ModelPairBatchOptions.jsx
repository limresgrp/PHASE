import { modelPairs } from '../../utils/modelPairs';

export default function ModelPairBatchOptions({ models, enabled, onEnabled, ids, onIds, references, onReferences, samples = [] }) {
  const selected = ids.filter((id) => models.some((m) => m.model_id === id));
  const count = selected.length * (selected.length - 1) / 2;
  return <details className="space-y-3 rounded border border-gray-700 p-3">
    <summary className="cursor-pointer text-sm font-semibold">Advanced: compare multiple models</summary>
    <label className="block text-sm"><input type="checkbox" checked={enabled} onChange={(e) => onEnabled(e.target.checked)} /> Run all unique model pairs</label>
    {enabled && <>
      <p className="text-xs text-gray-400">Each pair creates an independent analysis/job, also with one server worker. The same trajectory selection and settings are reused. Earlier selected model = A; later = B (ΔE = E(A) − E(B)).</p>
      <div className="max-h-48 space-y-1 overflow-auto">{models.map((m) => <label key={m.model_id} className="block text-sm"><input type="checkbox" checked={selected.includes(m.model_id)} onChange={(e) => onIds(e.target.checked ? [...selected, m.model_id] : selected.filter((id) => id !== m.model_id))} /> {m.name || m.model_id}</label>)}</div>
      <p className="text-sm">{count} independent analyses</p>
      {selected.length >= 2 && <p className="text-xs">{modelPairs(selected).map(([a, b]) => `${models.find((m) => m.model_id === a)?.name || a} vs ${models.find((m) => m.model_id === b)?.name || b}`).join('; ')}</p>}
      {references && selected.map((id) => <label className="block text-sm" key={id}>Reference MDs for {models.find((m) => m.model_id === id)?.name || id}
        <select multiple className="block w-full rounded border border-gray-600 bg-gray-900 p-2" size={Math.min(5, Math.max(2, samples.length))} value={references[id] || []} onChange={(e) => onReferences(id, Array.from(e.target.selectedOptions, (o) => o.value))}>
          {samples.map((s) => <option key={s.sample_id} value={s.sample_id}>{s.name || s.sample_id}</option>)}
        </select><span className="text-xs text-gray-400">Check suggestions; select disjoint endpoint reference groups. Ctrl/Cmd-click for multiple MDs.</span>
      </label>)}
    </>}
  </details>;
}
