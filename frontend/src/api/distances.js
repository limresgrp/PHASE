import { requestJSON } from './client';

function base(projectId, systemId, clusterId) {
  return `/projects/${encodeURIComponent(projectId)}/systems/${encodeURIComponent(systemId)}/metastable/clusters/${encodeURIComponent(clusterId)}/distance`;
}

export function fetchDistanceStatus(p, s, c) { return requestJSON(`${base(p, s, c)}/cache`); }
export function buildDistanceCache(p, s, c, payload) { return requestJSON(`${base(p, s, c)}/cache`, { method: 'POST', body: payload }); }
export function submitDistanceJs(p, s, c, payload) { return requestJSON(`${base(p, s, c)}/analyses`, { method: 'POST', body: payload }); }
export function fetchDistancePairs(p, s, c, options) {
  return requestJSON(`${base(p, s, c)}/pairs?${new URLSearchParams(options)}`);
}
export function fetchDistancePair(p, s, c, r, t, payload) {
  return requestJSON(`${base(p, s, c)}/pairs/${r}/${t}`, { method: 'POST', body: payload });
}
export function fetchDistancePresets(p, s, c) { return requestJSON(`${base(p, s, c)}/presets`); }
export function saveDistancePreset(p, s, c, payload) { return requestJSON(`${base(p, s, c)}/presets`, { method: 'POST', body: payload }); }
export function deleteDistancePreset(p, s, c, id) { return requestJSON(`${base(p, s, c)}/presets/${encodeURIComponent(id)}`, { method: 'DELETE' }); }
