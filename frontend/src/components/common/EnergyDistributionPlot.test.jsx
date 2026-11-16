import { buildEnergyDistributionPlot } from './EnergyDistributionPlot';

jest.mock('react-plotly.js', () => () => null);

const distribution = { label: 'MD', bins: [0, 1, 2], density: [2, 4], color: '#2563eb' };
const pose = { label: 'Crystal', values: [0.35], bins: [0, 1, 2], density: [1000, 0], color: '#dc2626', lineDash: 'dash' };

test.each(['histogram', 'curves'])('single-frame states are lines at the tallest distribution in %s mode', (mode) => {
  const reference = buildEnergyDistributionPlot({ series: [distribution], mode });
  const plot = buildEnergyDistributionPlot({ series: [pose, distribution], mode });
  const line = plot.data.find((trace) => trace.name === 'Crystal');
  const tallest = Math.max(...reference.data.flatMap((trace) => trace.y || []));
  expect(line.x).toEqual([0.35, 0.35]);
  expect(line.y).toEqual([0, tallest]);
  expect(line.line).toMatchObject({ color: '#dc2626', dash: 'dash' });
  expect(line.hovertemplate).toContain('single frame');
  expect(plot.data.filter((trace) => trace.type === 'bar')).toHaveLength(mode === 'histogram' ? 1 : 0);
  expect(plot.data.filter((trace) => trace.name?.startsWith('MD'))).toEqual(reference.data);
});

test('raw-value histogram heights participate in pose-line height without adding a pose histogram', () => {
  const plot = buildEnergyDistributionPlot({ series: [{ label: 'MD', values: [0, 0, 1, 1] }, { label: 'Pose', values: [0.5] }] });
  const line = plot.data.find((trace) => trace.name === 'Pose');
  const histogram = plot.data.find((trace) => trace.type === 'histogram');
  const curve = plot.data.find((trace) => trace.name === 'MD');
  const tallestBar = 2 / (4 * histogram.xbins.size);
  expect(line.y[1]).toBe(Math.max(tallestBar, ...curve.y));
  expect(plot.data.filter((trace) => trace.type === 'histogram')).toHaveLength(1);
});

test('multiple pose-only states use a finite fallback height and remain in the legend', () => {
  const plot = buildEnergyDistributionPlot({ series: [{ label: 'A', values: [-1] }, { label: 'B', values: [2] }] });
  expect(plot.data.map((trace) => trace.y)).toEqual([[0, 1], [0, 1]]);
  expect(plot.data.map((trace) => trace.name)).toEqual(['A', 'B']);
});

test('decorative peak-marker heights do not stretch single-frame state lines', () => {
  const plot = buildEnergyDistributionPlot({ series: [{ ...distribution, showPeakLine: true, peakLineHeight: 100 }, pose], mode: 'curves' });
  expect(plot.data.find((trace) => trace.name === 'Crystal').y[1]).toBe(3.5);
  expect(plot.layout.shapes[0].y1).toBe(100);
});
