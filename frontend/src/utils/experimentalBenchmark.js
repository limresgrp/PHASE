export const benchmarkNumber = (value) => value === null || value === undefined ? NaN : Number(value);

export function axisPairs(count) {
  const result = [];
  for (let i = 0; i < count; i += 1) for (let j = i + 1; j < count; j += 1) result.push([i, j]);
  return result;
}

const palette = ['#2563eb', '#e76f51', '#16a085', '#a855a1', '#d4a017', '#6b7280', '#d63384', '#4f46e5'];
const symbols = ['circle', 'square', 'diamond', 'triangle-up', 'cross', 'star'];

export function benchmarkTraces(result, xAxis, yAxis, variant = 'common') {
  if (!result) return [];
  const categories = [...result.label_scheme.categories];
  if (result.scores.some((s) => !s.label)) categories.push(null);
  return categories.map((category, i) => {
    const points = result.states.map((state) => {
      const x = result.scores.find((s) => s.state_id === state.state_id && s.axis === xAxis);
      const y = result.scores.find((s) => s.state_id === state.state_id && s.axis === yAxis);
      return { state, x, y };
    }).filter(({ x, y }) => x && y && x.label === category && Number.isFinite(benchmarkNumber(x[`${variant}_rank`])) && Number.isFinite(benchmarkNumber(y[`${variant}_rank`])));
    return { type: 'scatter', mode: 'markers', name: category || 'Unlabelled',
      x: points.map((p) => p.x[`${variant}_rank`]), y: points.map((p) => p.y[`${variant}_rank`]),
      customdata: points.map(({ state, x, y }) => [state.name, category || 'Unlabelled', x.raw_score, x.common_score, x.raw_rank, x.common_rank, x.observed_nodes, x.observed_edges, y.raw_score, y.common_score, y.raw_rank, y.common_rank, y.observed_nodes, y.observed_edges, state.state_id]),
      marker: { color: category === null ? '#6b7280' : palette[i % palette.length], symbol: symbols[i % symbols.length], size: 12, line: { color: '#ffffff', width: 1 } },
      hovertemplate: '<b>%{customdata[0]}</b> (%{customdata[14]})<br>User label: %{customdata[1]}<br>X raw/common: %{customdata[2]:.4f} / %{customdata[3]:.4f}<br>X ranks raw/common: %{customdata[4]:.3f} / %{customdata[5]:.3f}<br>X observed nodes/edges: %{customdata[6]} / %{customdata[7]}<br>Y raw/common: %{customdata[8]:.4f} / %{customdata[9]:.4f}<br>Y ranks raw/common: %{customdata[10]:.3f} / %{customdata[11]:.3f}<br>Y observed nodes/edges: %{customdata[12]} / %{customdata[13]}<extra></extra>',
    };
  }).filter((trace) => trace.x.length);
}

export function discordancePairs(result, globalThreshold, localThreshold) {
  return (result?.discordance || []).filter((p) => Number.isFinite(benchmarkNumber(p.global_rank_difference))
    && Number.isFinite(benchmarkNumber(p.local_distance)) && p.global_rank_difference <= globalThreshold && p.local_distance >= localThreshold)
    .sort((a, b) => b.local_distance - a.local_distance);
}
