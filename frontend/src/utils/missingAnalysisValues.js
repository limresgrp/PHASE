// JSON carries unavailable scientific values as null. Number(null) is zero,
// so restore NaN before plotting, filtering, centering or combining metrics.
const metricKey = /^(?:node_js$|edge_js$|js_node_|js_edge_|q_residue|q_edge|frustration_|D_residue$|D_edge$|p_node|p_edge|global_node_|global_edge_|frame_score_)/;
const numericArray = (value) => Array.isArray(value)
  ? value.map(numericArray)
  : value === null ? NaN : value;

export function preserveMissingAnalysisValues(payload) {
  if (!payload || typeof payload !== 'object') return payload;
  const convert = (data, metricContainer = false) => Object.fromEntries(Object.entries(data).map(([key, value]) => [
    key, (metricContainer || metricKey.test(key)) && Array.isArray(value) ? numericArray(value)
      : value && typeof value === 'object' && !Array.isArray(value) ? convert(value, key === 'node' || key === 'edge') : value,
  ]));
  return payload.data && typeof payload.data === 'object'
    ? { ...payload, data: convert(payload.data) } : convert(payload);
}
