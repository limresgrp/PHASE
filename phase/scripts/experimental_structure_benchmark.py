"""Run reproducible experimental structure benchmarks and manage analysis labels."""
import argparse
import json
import os
from pathlib import Path


def interactive_config(store, project, system, cluster):
    from phase.analysis.experimental_benchmark import benchmark_context
    from phase.services.analysis_annotations import save_scheme
    context = benchmark_context(store, project, system, cluster)
    def choose(title, rows, identity, multiple=False):
        print(title)
        for i, row in enumerate(rows, 1):
            print(f"  [{i}] {row.get('name') or row[identity]}")
        text = input('Select numbers (comma separated): ' if multiple else 'Select number: ').strip()
        indices = [int(n.strip()) - 1 for n in text.split(',') if n.strip()]
        if not indices or any(i < 0 or i >= len(rows) for i in indices) or len(set(indices)) != len(indices) or (not multiple and len(indices) != 1):
            raise ValueError('Select valid, distinct menu numbers.')
        return [rows[i] for i in indices] if multiple else rows[indices[0]]
    states = choose('Select PDB-only/static states (assignments must contain exactly one frame):',
                    [s for s in context['states'] if s['n_frames'] <= 1], 'state_id', True)
    state_ids = [s['state_id'] for s in states]
    if len(states) < 2:
        raise ValueError('Select at least two states.')
    saved = context['label_schemes']
    use_saved = saved and input('Use a saved label scheme? [Y/n]: ').lower() != 'n'
    config = {'state_ids': state_ids, 'contrasts': [], 'state_samples': {}}
    if use_saved:
        scheme = choose('Select label scheme:', saved, 'scheme_id')
        config['scheme_id'] = scheme['scheme_id']
    else:
        name = input('Label scheme name: ').strip()
        kind = input('Scheme kind (categorical/ordinal) [categorical]: ').strip() or 'categorical'
        categories = [v.strip() for v in input('Categories, comma separated (for ordinal, explicit low-to-high order): ').split(',')]
        labels = {}
        for state in states:
            label = input(f"Label for {state['name']} (blank = unlabelled): ").strip()
            if label:
                labels[state['state_id']] = label
        scheme = save_scheme(store, project, system, {'name': name, 'kind': kind, 'categories': categories, 'labels': labels})
        config['scheme_id'] = scheme['scheme_id']
    print(f"Label scheme: {scheme['name']} ({scheme['kind']}); supplied category list/order: {', '.join(scheme['categories'])}")
    for state in states:
        if len(state['assigned_samples']) > 1:
            config['state_samples'][state['state_id']] = choose(f"Choose single-frame assignment for {state['name']}:", state['assigned_samples'], 'sample_id')['sample_id']
    while True:
        a = choose('Select model A:', context['models'], 'model_id')
        b = choose('Select model B:', context['models'], 'model_id')
        config['contrasts'].append({'model_a_id': a['model_id'], 'model_b_id': b['model_id']})
        if input('Add another contrast? [y/N]: ').lower() != 'y':
            break
    config['hotspots_mode'] = input('Hotspots mode (parameter/ensemble) [parameter]: ').strip() or 'parameter'
    if config['hotspots_mode'] == 'ensemble':
        for i, contrast in enumerate(config['contrasts'], 1):
            for side in ('a', 'b'):
                refs = choose(f'Axis {i}, endpoint {side.upper()} reference MD trajectories:', context['md_samples'], 'sample_id', True)
                contrast[f'reference_{side}'] = [r['sample_id'] for r in refs]
        config['block_length'] = int(input('Bootstrap block length [50]: ') or 50)
        config['bootstrap'] = int(input('Bootstrap replicates (0 or >=50) [200]: ') or 200)
    if input('Compute AUC for two supplied label categories? [y/N]: ').lower() == 'y':
        config['positive_label'] = input('Positive category (higher score is positive): ').strip()
        config['negative_label'] = input('Negative category: ').strip()
        config['best_oriented'] = input('Also report best-oriented AUC separately? [y/N]: ').lower() == 'y'
    config['fingerprint_mode'] = input('Fingerprint mode (standardized/raw) [standardized]: ').strip() or 'standardized'
    config['distance'] = input('PCoA/local distance (euclidean/correlation) [euclidean]: ').strip() or 'euclidean'
    config['top_k'] = int(input('Hotspot comparison top-k (ties included) [20]: ') or 20)
    config['global_threshold'] = float(input('Global rank RMS difference threshold [0.15]: ') or .15)
    local = input('Local distance threshold (blank = 75th percentile): ').strip()
    config['local_threshold'] = float(local) if local else None
    config['seed'] = int(input('Random seed [0]: ') or 0)
    config['name'] = input('Benchmark name [Experimental structure benchmark]: ').strip() or 'Experimental structure benchmark'
    return config


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', default=os.getenv('PHASE_DATA_ROOT', 'data'))
    parser.add_argument('--project-id', required=True)
    parser.add_argument('--system-id', required=True)
    parser.add_argument('--cluster-id')
    parser.add_argument('--config', help='JSON containing state_ids, contrasts and scheme_id or label_scheme; all benchmark options supported')
    parser.add_argument('--interactive', action='store_true')
    parser.add_argument('--progress', action='store_true')
    parser.add_argument('--list-label-schemes', action='store_true')
    parser.add_argument('--save-label-scheme', help='JSON scheme file; include scheme_id to update')
    parser.add_argument('--delete-label-scheme')
    args = parser.parse_args(argv)
    root = Path(args.root).expanduser().resolve()
    os.environ['PHASE_DATA_ROOT'] = str(root)
    from phase.services.project_store import ProjectStore
    from phase.services.analysis_annotations import list_schemes, save_scheme, delete_scheme
    from phase.analysis.experimental_benchmark import run_experimental_benchmark
    store = ProjectStore(base_dir=root / 'projects')
    if args.list_label_schemes:
        print(json.dumps(list_schemes(store, args.project_id, args.system_id), indent=2)); return 0
    if args.save_label_scheme:
        value = json.loads(Path(args.save_label_scheme).read_text())
        print(json.dumps(save_scheme(store, args.project_id, args.system_id, value, value.get('scheme_id')), indent=2)); return 0
    if args.delete_label_scheme:
        delete_scheme(store, args.project_id, args.system_id, args.delete_label_scheme); return 0
    if not args.cluster_id or bool(args.config) == args.interactive:
        parser.error('Supply --cluster-id and one of --config or --interactive.')
    config = interactive_config(store, args.project_id, args.system_id, args.cluster_id) if args.interactive else json.loads(Path(args.config).read_text())
    bar = None
    def progress(message, current, total):
        nonlocal bar
        if args.progress:
            from tqdm import tqdm
            if bar is None:
                bar = tqdm(total=total, unit='stage')
            bar.total = total; bar.set_description_str(message); bar.update(current - bar.n)
    try:
        result = run_experimental_benchmark(args.project_id, args.system_id, args.cluster_id, store=store, progress_callback=progress, **config)
    except (ValueError, FileNotFoundError) as exc:
        parser.exit(2, f'Benchmark error: {exc}\n')
    finally:
        if bar is not None:
            bar.close()
    print(json.dumps(result, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
