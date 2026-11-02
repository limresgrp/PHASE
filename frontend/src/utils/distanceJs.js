// Saved JS is in nats. Display filters and effect sizes use normalized units.
export const JS_MAX = Math.log(2);
const finite = (value) => typeof value === 'number' && Number.isFinite(value);

export function distanceJsRows(data, { bins, minFrames = 20, minChange = 0, sort = 'change' } = {}) {
  if (!data?.edge_global_js) return [];
  const selected = bins || (data.distance_bin_labels || []).map((_, index) => index);
  return data.edge_global_js.map((global, index) => {
    const supported = selected.filter((b) => finite(data.edge_bin_js?.[b]?.[index])
      && data.edge_bin_counts_a?.[b]?.[index] >= minFrames && data.edge_bin_counts_b?.[b]?.[index] >= minFrames);
    const values = supported.map((b) => data.edge_bin_js[b][index] / JS_MAX);
    const change = values.length >= 2 ? Math.max(...values) - Math.min(...values) : null;
    const value = sort === 'change' ? change : sort === 'global' ? (finite(global) ? global / JS_MAX : null)
      : supported.includes(Number(sort)) ? data.edge_bin_js[Number(sort)][index] / JS_MAX : null;
    return { index, global, supported, change, value };
  }).filter((row) => row.supported.length > 0 && finite(row.value) && (minChange === 0 || (row.change !== null && row.change >= minChange)))
    .sort((a, b) => b.value - a.value || a.index - b.index);
}

export function distanceJsResidues(rows, data) {
  const residues = new Map();
  rows.filter((row) => row.change !== null).forEach((row) => {
    for (const i of [data.edge_index[0][row.index], data.edge_index[1][row.index]]) {
      if (!residues.has(i)) residues.set(i, { index: i, key: data.residue_keys[i], sum: 0, count: 0, max: 0 });
      const residue = residues.get(i);
      residue.sum += row.change; residue.count += 1; residue.max = Math.max(residue.max, row.change);
    }
  });
  return [...residues.values()].map((row) => ({ ...row, mean: row.sum / row.count })).sort((a, b) => b.mean - a.mean || a.index - b.index);
}

export function binsBelowCutoff(metadata, count, cutoff) {
  if (!cutoff) return Array.from({ length: count }, (_, i) => i);
  if (!metadata?.distance_available) return [];
  // Retain complete saved bins only; never imply an arbitrary new boundary.
  return (metadata.thresholds || []).map((upper, i) => upper <= Number(cutoff) ? i : -1).filter((i) => i >= 0);
}
