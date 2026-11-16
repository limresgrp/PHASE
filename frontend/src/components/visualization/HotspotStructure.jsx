import { useCallback, useEffect, useRef, useState } from 'react';
import { Shape } from 'molstar/lib/mol-model/shape';
import { Mesh } from 'molstar/lib/mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from 'molstar/lib/mol-geo/geometry/mesh/mesh-builder';
import { addCylinder } from 'molstar/lib/mol-geo/geometry/mesh/builder/cylinder';
import { ShapeRepresentation } from 'molstar/lib/mol-repr/shape/representation';
import { MolScriptBuilder as MS } from 'molstar/lib/mol-script/language/builder';
import { Script } from 'molstar/lib/mol-script/script';
import { StructureSelection, StructureElement, StructureProperties } from 'molstar/lib/mol-model/structure';
import { setStructureOverpaint, clearStructureOverpaint } from 'molstar/lib/mol-plugin-state/helpers/structure-overpaint';
import MolstarTrajectoryViewer from './MolstarTrajectoryViewer';
import { requestJSON } from '../../api/client';

import { hotspotColor } from '../../utils/hotspotMetrics';

export default function HotspotStructure({ projectId, systemId, clusterId, system, data, nodes, edges, signed, maximum, focusKey }) {
  const viewer = useRef(null);
  const queue = useRef(Promise.resolve());
  const loaded = useRef('');
  const mapping = useRef([]);
  const overlay = useRef(null);
  const hover = useRef(new Map());
  const [ready, setReady] = useState(false);
  const [stateId, setStateId] = useState('');
  const [error, setError] = useState('');
  const [message, setMessage] = useState('');
  const states = Object.values(system?.states || {}).filter((state) => state.pdb_file);
  const selected = stateId || (system?.reference_pdb_file ? '_reference' : states[0]?.state_id);
  const onStatus = useCallback((status) => { if (status === 'ready') setReady(true); }, []);

  useEffect(() => {
    if (!ready) return undefined;
    const plugin = viewer.current?.getPlugin();
    if (!plugin) return undefined;
    const provider = {
      priority: 200,
      label: (loci) => {
        if (!StructureElement.Loci.is(loci)) return undefined;
        const location = StructureElement.Loci.getFirstLocation(loci);
        if (!location) return undefined;
        const identity = JSON.stringify([
          StructureProperties.chain.auth_asym_id(location),
          StructureProperties.residue.auth_seq_id(location),
          StructureProperties.residue.pdbx_PDB_ins_code(location) || '',
        ]);
        const row = hover.current.get(identity);
        if (!row) return undefined;
        const interval = Number.isFinite(row.low) && Number.isFinite(row.high)
          ? `${row.low.toPrecision(3)} / ${row.high.toPrecision(3)}` : 'unavailable';
        return `Hotspot effect: ${row.effect.toPrecision(4)}; displayed: ${row.color.toPrecision(4)}; 95% CI: ${interval}; observed A/B: ${row.countA ?? '-'} / ${row.countB ?? '-'}`;
      },
      group: (label) => `phase-hotspots:${label}`,
    };
    plugin.managers.lociLabels.addProvider(provider);
    return () => {
      try { plugin.managers.lociLabels.removeProvider(provider); } catch { /* Viewer may already be disposed. */ }
    };
  }, [ready]);

  useEffect(() => {
    if (!ready || !selected) return undefined;
    let cancelled = false;
    queue.current = queue.current.catch(() => {}).then(async () => {
      if (cancelled) return;
      const plugin = viewer.current?.getPlugin();
      if (!plugin) return;
      setError(''); setMessage('Updating structure...');
      hover.current = new Map();
      if (overlay.current) { plugin.canvas3d?.remove(overlay.current); overlay.current.destroy(); overlay.current = null; }
      const signature = `${projectId}/${systemId}/${clusterId}/${selected}`;
      if (loaded.current !== signature) {
        loaded.current = '';
        const base = `/projects/${encodeURIComponent(projectId)}/systems/${encodeURIComponent(systemId)}`;
        const map = await requestJSON(`${base}/metastable/clusters/${encodeURIComponent(clusterId)}/hotspots/structure-map/${encodeURIComponent(selected)}`);
        if (cancelled) return;
        await viewer.current.loadStructure({ structureUrl: selected === '_reference' ? `/api/v1${base}/reference` : `/api/v1${base}/structures/${encodeURIComponent(selected)}`, structureName: 'reference.pdb',
          presetParams: { model: { modelIndex: 0 }, structure: { name: 'model', params: {} } } });
        // Use an opaque neutral cartoon; mapping coordinates are for the first
        // deposited model, not a symmetry-expanded biological assembly.
        const root = plugin.managers.structure.hierarchy.current.structures[0];
        if (!root) throw new Error('No structure was loaded.');
        await plugin.managers.structure.component.clear([root]);
        const component = await plugin.builders.structure.tryCreateComponentFromExpression(root.cell, MS.struct.generator.all(), 'hotspot-base');
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
      const residues = new Map(nodes.map((row) => [row.index, row.color]));
      for (const row of nodes) {
        const atom = mapping.current[row.index];
        if (atom) hover.current.set(JSON.stringify([atom.chain || '', atom.resid, atom.insertion_code || '']), row);
      }
      const mesh = MeshBuilder.createState();
      const labels = [];
      const colors = [];
      let skipped = 0;
      for (const row of edges) {
        const [r, s] = data.edges[row.index];
        const a = mapping.current[r], b = mapping.current[s];
        if (!a || !b) { skipped += 1; continue; }
        if (a.key !== data.residue_keys[r] || b.key !== data.residue_keys[s]) throw new Error('Structure mapping differs from this analysis.');
        mesh.currentGroup = labels.length;
        addCylinder(mesh, a.position, b.position, 1, { radiusTop: .18, radiusBottom: .18, radialSegments: 8 });
        labels.push(`${row.label}: effect ${row.effect.toPrecision(4)}; displayed ${row.color.toPrecision(4)}`);
        colors.push(hotspotColor(row.color, maximum, signed));
      }
      for (const [i, value] of residues) {
        if (cancelled) return;
        const atom = mapping.current[i];
        if (!atom) { skipped += 1; continue; }
        const tests = { 'residue-test': MS.core.logic.and([
          MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_seq_id(), atom.resid]),
          MS.core.rel.eq([MS.struct.atomProperty.macromolecular.pdbx_PDB_ins_code(), atom.insertion_code || '']),
        ]), 'chain-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_asym_id(), atom.chain || '']) };
        const expression = MS.struct.generator.atomGroups(tests);
        await setStructureOverpaint(plugin, components, hotspotColor(value, maximum, signed), (s) =>
          StructureSelection.toLociWithSourceUnits(Script.getStructureSelection(expression, s)));
      }
      if (labels.length && !cancelled) {
        const shape = Shape.create('Endpoint discriminative edges', {}, MeshBuilder.getMesh(mesh), (i) => colors[i], () => 1, (i) => labels[i]);
        const repr = ShapeRepresentation(() => shape, Mesh.Utils);
        await plugin.runTask(repr.createOrUpdate({ alpha: 1 }, {}));
        if (cancelled) repr.destroy();
        else { plugin.canvas3d.add(repr); overlay.current = repr; }
      }
      if (!cancelled) {
        setMessage(`${residues.size} selected residues, ${labels.length} edges; ${skipped} missing-coordinate items skipped. Unselected residues remain opaque gray.`);
        const atom = mapping.current.find((item) => item?.key === focusKey);
        const structureData = structure.cell.obj?.data;
        if (atom && structureData) {
          const expression = MS.struct.generator.atomGroups({
            'chain-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_asym_id(), atom.chain || '']),
            'residue-test': MS.core.logic.and([
              MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_seq_id(), atom.resid]),
              MS.core.rel.eq([MS.struct.atomProperty.macromolecular.pdbx_PDB_ins_code(), atom.insertion_code || '']),
            ]),
          });
          plugin.managers.camera.focusLoci(StructureSelection.toLociWithSourceUnits(Script.getStructureSelection(expression, structureData)));
        }
      }
    }).catch((err) => { if (!cancelled) { setError(err.message); setMessage(''); } });
    return () => { cancelled = true; };
  }, [ready, selected, projectId, systemId, clusterId, data, nodes, edges, signed, maximum, focusKey]);

  useEffect(() => () => { overlay.current?.destroy(); overlay.current = null; }, []);

  return <div className="space-y-3">
    <h3 className="font-semibold">Discriminative hotspot network</h3>
    <label className="text-sm">Reference PDB <select className="rounded bg-gray-900 p-2" value={selected || ''} onChange={(e) => setStateId(e.target.value)}>
      {system?.reference_pdb_file && <option value="_reference">System reference</option>}
      {states.map((state) => <option key={state.state_id} value={state.state_id}>{state.name || state.state_id}</option>)}
    </select></label>
    <p className="text-xs text-gray-400">{signed ? 'Red = A-like; blue = B-like relative to the balanced endpoint midpoint.' : 'Pale = low importance; orange = high importance. Unsigned strength has no A/B direction.'} Gray = not selected. Click a table row to focus a residue; hover a line for its pair and scores.</p>
    {!selected && <p className="text-sm">No reference PDB is available. Add a structure on the System page.</p>}
    {error && <p role="alert" className="text-sm text-red-400">{error}</p>}
    <p className="text-xs text-gray-400">{message}</p>
    <MolstarTrajectoryViewer ref={viewer} height={580} onStatusChange={onStatus} onError={setError} />
  </div>;
}
