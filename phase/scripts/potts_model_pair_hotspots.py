"""Command-line interface for independent endpoint hotspot prioritization."""
import argparse
import json
from pathlib import Path
import os


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', default=os.getenv('PHASE_DATA_ROOT', 'data'))
    for field in ('project-id', 'system-id', 'cluster-id'):
        parser.add_argument('--' + field, required=True)
    parser.add_argument('--model-a-id')
    parser.add_argument('--model-b-id')
    parser.add_argument('--model-ids', help='Advanced: comma-separated models; run all unique pairs in order')
    parser.add_argument('--model-reference', action='append', default=[], help='Per-model MD references as model_id:sample_id,sample_id (batch mode)')
    parser.add_argument('--mode', choices=['parameter', 'ensemble'], default='ensemble')
    parser.add_argument('--reference-a', help='Comma-separated assigned MD sample IDs; omitted = infer from model states')
    parser.add_argument('--reference-b', help='Comma-separated assigned MD sample IDs')
    parser.add_argument('--sample-ids', default='', help='Additional trajectories for endpoint preference visualization')
    parser.add_argument('--block-length', type=int, default=50, help='Contiguous assigned frames per bootstrap block')
    parser.add_argument('--bootstrap', type=int, default=200, help='Replicates; 0 disables uncertainty, otherwise >=50')
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--progress', action='store_true')
    args = vars(parser.parse_args(argv))
    from phase.potts.model_pairs import model_pairs, run_model_pairs
    try:
        ids = args.pop('model_ids')
        pairs = model_pairs(ids.split(',') if ids else None, args.pop('model_a_id'), args.pop('model_b_id'))
        refs = {}
        for item in args.pop('model_reference'):
            model, samples = item.split(':', 1)
            refs[model] = [v.strip() for v in samples.split(',') if v.strip()]
        if len(pairs) > 1 and (args['reference_a'] is not None or args['reference_b'] is not None):
            raise ValueError('Use --model-reference for batch references, not --reference-a/b.')
        if set(refs) - {m for pair in pairs for m in pair}:
            raise ValueError('References specified for an unselected model.')
    except ValueError as exc:
        parser.error(str(exc))
    os.environ['PHASE_DATA_ROOT'] = str(Path(args.pop('root')).expanduser().resolve())
    # Import after resolving the root; ProjectStore defaults are evaluated on import.
    from phase.services.project_store import ProjectStore
    from phase.potts.hotspots import run_hotspot_analysis
    show = args.pop('progress')
    bar = None
    def progress(message, current, total):
        nonlocal bar
        if show:
            from tqdm import tqdm
            if bar is None:
                bar = tqdm(total=total, unit='stage')
            bar.total = total
            bar.set_description_str(message)
            bar.update(current - bar.n)
    for field in ('reference_a', 'reference_b', 'sample_ids'):
        if args[field] is not None:
            args[field] = [v.strip() for v in args[field].split(',') if v.strip()]
    store = ProjectStore(base_dir=Path(os.environ['PHASE_DATA_ROOT']) / 'projects')
    def run(a, b):
        nonlocal bar
        bar = None
        options = dict(args, model_a_id=a, model_b_id=b)
        if a in refs:
            options['reference_a'] = refs[a]
        if b in refs:
            options['reference_b'] = refs[b]
        try:
            result = run_hotspot_analysis(store=store, progress_callback=progress, **options)
        finally:
            if bar is not None:
                bar.close()
        print(json.dumps(result, indent=2))
    return run_model_pairs(pairs, run)


if __name__ == '__main__':
    raise SystemExit(main())
