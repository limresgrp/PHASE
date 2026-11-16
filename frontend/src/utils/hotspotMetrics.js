export const numberOrNaN = (value) => value === null || value === undefined ? NaN : Number(value);

export function hotspotRows(data, { metric = 'separation', preferenceSample = 0, view = 'importance',
  blend = false, alpha = .75, aggregation = 'mean', minimum = 0, top = 30, percentile = 0,
  reliable = false, direction = 'all', selection = [], selectionMode = 'exclude', labels = {} } = {}) {
  if (!data) return { nodes: [], edges: [] };
  const keys = data.residue_keys || [];
  const get = (kind, name) => (data[`${kind}_${name}`] || []).map(numberOrNaN);
  const nodeEffect = get('node', metric), edgeEffect = get('edge', metric);
  const preference = (kind) => (data[`${kind}_preference`]?.[preferenceSample] || []).map(numberOrNaN);
  let nodeColor = view === 'preference' ? preference('node') : nodeEffect.map(Math.abs);
  const edgeColor = view === 'preference' ? preference('edge') : edgeEffect.map(Math.abs);
  if (blend) {
    const incident = keys.map(() => []);
    (data.edges || []).forEach(([a, b], i) => {
      if (Number.isFinite(edgeColor[i])) { incident[a]?.push(edgeColor[i]); incident[b]?.push(edgeColor[i]); }
    });
    nodeColor = nodeColor.map((value, i) => {
      const terms = incident[i];
      if (!terms.length) return value;
      const sum = terms.reduce((a, b) => a + b, 0);
      const edge = aggregation === 'sum' ? sum : sum / terms.length;
      if (alpha === 1) return edge;
      return Number.isFinite(value) ? (1 - alpha) * value + alpha * edge : NaN;
    });
  }
  const selectionSet = new Set(selection.map(Number));
  const selected = (indices) => !selectionSet.size || (selectionMode === 'include'
    ? indices.some((i) => selectionSet.has(i)) : indices.every((i) => !selectionSet.has(i)));
  const build = (kind, effects, colors, indices) => {
    const low = get(kind, 'ci_low'), high = get(kind, 'ci_high');
    const countA = get(kind, 'counts_a'), countB = get(kind, 'counts_b');
    const rows = effects.map((effect, index) => ({ index, indices: indices[index], effect, color: colors[index],
      label: indices[index].map((i) => labels[keys[i]] || keys[i]).join(' - '),
      low: low[index], high: high[index], countA: countA[index], countB: countB[index],
    })).filter((row) => Number.isFinite(row.effect) && Number.isFinite(row.color) && Math.abs(row.effect) >= minimum
      && selected(row.indices) && (!reliable || (row.low > 0 || row.high < 0))
      && (direction === 'all' || (direction === 'positive' ? row.effect > 0 : row.effect < 0)));
    rows.sort((a, b) => Math.abs(b.effect) - Math.abs(a.effect));
    const keep = Math.max(1, Math.ceil(rows.length * (1 - Math.min(100, Math.max(0, percentile)) / 100)));
    return rows.slice(0, Math.min(top || rows.length, keep));
  };
  return { nodes: build('node', nodeEffect, nodeColor, keys.map((_, i) => [i])),
    edges: build('edge', edgeEffect, edgeColor, data.edges || []) };
}

export function hotspotColor(value, maximum, signed = false) {
  const t = Math.min(1, Math.abs(value) / Math.max(maximum, 1e-12));
  const end = signed ? (value >= 0 ? [210, 45, 55] : [45, 105, 210]) : [220, 95, 25];
  return end.reduce((color, channel) => (color << 8) | Math.round(238 + (channel - 238) * t), 0);
}
