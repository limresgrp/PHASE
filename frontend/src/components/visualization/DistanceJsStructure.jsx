import { useCallback, useEffect, useRef, useState } from 'react';
import { Shape } from 'molstar/lib/mol-model/shape';
import { Mesh } from 'molstar/lib/mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from 'molstar/lib/mol-geo/geometry/mesh/mesh-builder';
import { addCylinder } from 'molstar/lib/mol-geo/geometry/mesh/builder/cylinder';
import { ShapeRepresentation } from 'molstar/lib/mol-repr/shape/representation';
import { MolScriptBuilder as MS } from 'molstar/lib/mol-script/language/builder';
import { Script } from 'molstar/lib/mol-script/script';
import { StructureSelection } from 'molstar/lib/mol-model/structure';
import { setStructureOverpaint, clearStructureOverpaint } from 'molstar/lib/mol-plugin-state/helpers/structure-overpaint';
import MolstarTrajectoryViewer from './MolstarTrajectoryViewer';
import { requestJSON } from '../../api/client';

function color(value, maximum) {
  const t = Math.max(0, Math.min(1, value / Math.max(.001, maximum)));
  return (Math.round(34 + 205 * t) << 16) | (Math.round(197 - 129 * t) << 8) | Math.round(94 - 26 * t);
}

export default function DistanceJsStructure({ projectId, systemId, clusterId, system, data, rows }) {
  const viewer = useRef(null);
  const queue = useRef(Promise.resolve());
  const loaded = useRef('');
  const mapping = useRef([]);
  const overlay = useRef(null);
  const [ready, setReady] = useState(false);
  const [stateId, setStateId] = useState('');
  const [limit, setLimit] = useState(40);
  const [saturation, setSaturation] = useState(.3);
  const [error, setError] = useState('');
  const [message, setMessage] = useState('');
  const states = Object.values(system?.states || {}).filter((state) => state.pdb_file);
  const selected = stateId || states[0]?.state_id;
  const onStatus = useCallback((status) => { if (status === 'ready') setReady(true); }, []);

  useEffect(() => {
    if (!ready || !selected) return undefined;
    let cancelled = false;
    queue.current = queue.current.catch(() => {}).then(async () => {
      if (cancelled) return;
      const plugin = viewer.current?.getPlugin();
      if (!plugin) return;
      setError(''); setMessage('Updating structure…');
      if (overlay.current) { plugin.canvas3d?.remove(overlay.current); overlay.current.destroy(); overlay.current = null; }
      const signature = `${projectId}/${systemId}/${clusterId}/${selected}`;
      if (loaded.current !== signature) {
        loaded.current = '';
        const base = `/projects/${encodeURIComponent(projectId)}/systems/${encodeURIComponent(systemId)}`;
        const map = await requestJSON(`${base}/metastable/clusters/${encodeURIComponent(clusterId)}/distance/structure-map/${encodeURIComponent(selected)}`);
        if (cancelled) return;
        await viewer.current.loadStructure({ structureUrl: `/api/v1${base}/structures/${encodeURIComponent(selected)}`, structureName: 'reference.pdb',
          presetParams: { model: { modelIndex: 0 }, structure: { name: 'model', params: {} } } });
        // Use an opaque neutral cartoon; mapping coordinates are for the first
        // deposited model, not a symmetry-expanded biological assembly.
        const root = plugin.managers.structure.hierarchy.current.structures[0];
        if (!root) throw new Error('No structure was loaded.');
        await plugin.managers.structure.component.clear([root]);
        const component = await plugin.builders.structure.tryCreateComponentFromExpression(root.cell, MS.struct.generator.all(), 'distance-js-base');
        if (!component) throw new Error('Could not create the reference cartoon.');
        await plugin.builders.structure.representation.addRepresentation(component, { type: 'cartoon', typeParams: { alpha: 1 }, color: 'uniform', colorParams: { value: 0x9ca3af } });
        mapping.current = map.residues;
        loaded.current = signature;
      }
      if (cancelled) return;
      const structure = plugin.managers.structure.hierarchy.current.structures[0];
      if (!structure) throw new Error('No structure was loaded.');
      const components = structure.components;
      await clearStructureOverpaint(plugin, components);
      const shown = rows.filter((row) => row.change !== null).slice(0, Math.max(1, Number(limit) || 1));
      const residues = new Map();
      const mesh = MeshBuilder.createState();
      const labels = [];
      const colors = [];
      for (const row of shown) {
        const r = data.edge_index[0][row.index], s = data.edge_index[1][row.index];
        const a = mapping.current[r], b = mapping.current[s];
        if (!a || !b || a.key !== data.residue_keys[r] || b.key !== data.residue_keys[s]) throw new Error('Structure residue mapping differs from this analysis.');
        mesh.currentGroup = labels.length;
        addCylinder(mesh, a.position, b.position, 1, { radiusTop: .18, radiusBottom: .18, radialSegments: 8 });
        labels.push(`${a.resname}${a.resid} — ${b.resname}${b.resid}: normalized JS range ${row.change.toFixed(3)}`);
        colors.push(color(row.change, saturation));
        for (const i of [r, s]) {
          const values = residues.get(i) || []; values.push(row.change); residues.set(i, values);
        }
      }
      for (const [i, values] of residues) {
        if (cancelled) return;
        const atom = mapping.current[i];
        const tests = { 'residue-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_seq_id(), atom.resid]) };
        if (atom.chain) tests['chain-test'] = MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_asym_id(), atom.chain]);
        const expression = MS.struct.generator.atomGroups(tests);
        await setStructureOverpaint(plugin, components, color(values.reduce((sum, v) => sum + v, 0) / values.length, saturation), (s) =>
          StructureSelection.toLociWithSourceUnits(Script.getStructureSelection(expression, s)));
      }
      if (shown.length && !cancelled) {
        const shape = Shape.create('Distance-dependent JS', {}, MeshBuilder.getMesh(mesh), (i) => colors[i], () => 1, (i) => labels[i]);
        const repr = ShapeRepresentation(() => shape, Mesh.Utils);
        await plugin.runTask(repr.createOrUpdate({ alpha: 1 }, {}));
        if (cancelled) repr.destroy();
        else { plugin.canvas3d.add(repr); overlay.current = repr; }
      }
      if (!cancelled) setMessage(`${shown.length} ranked edges; residue colors average incident displayed-edge changes. Unselected residues retain their base color.`);
    }).catch((err) => { if (!cancelled) { setError(err.message); setMessage(''); } });
    return () => { cancelled = true; };
  }, [ready, selected, projectId, systemId, clusterId, data, rows, limit, saturation]);

  useEffect(() => () => { overlay.current?.destroy(); overlay.current = null; }, []);

  return <div className="space-y-3">
    <h3 className="font-semibold">Distance-dependent JS on the reference structure</h3>
    <div className="flex flex-wrap gap-3 text-xs">
      <label>Reference PDB <select className="rounded bg-gray-900 p-2" value={selected || ''} onChange={(e) => setStateId(e.target.value)}>{states.map((state) => <option key={state.state_id} value={state.state_id}>{state.name || state.state_id}</option>)}</select></label>
      <label>Top edges <input className="w-20 rounded bg-gray-900 p-2" type="number" min="1" max="200" value={limit} onChange={(e) => setLimit(Math.min(200, Number(e.target.value)))} /></label>
      <label>Red at JS range ≥ <input className="w-20 rounded bg-gray-900 p-2" type="number" min=".001" max="1" step=".05" value={saturation} onChange={(e) => setSaturation(Number(e.target.value))} /></label>
    </div>
    <p className="text-xs text-gray-400">Green = small change; red = large change (normalized 0–1). Lines connect Cα atoms in this reference PDB, not a distance-bin average structure. Hover lines for pair identities. Colors describe effect size, not statistical significance. This view needs two supported bins; global-only analyses have no JS-range edges.</p>
    {error && <p role="alert" className="text-sm text-red-400">{error}</p>}
    <p className="text-xs text-gray-400">{message}</p>
    <MolstarTrajectoryViewer ref={viewer} height={580} onStatusChange={onStatus} onError={setError} />
  </div>;
}
