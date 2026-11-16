"""Command-line interface for independent endpoint hotspot prioritization."""
import argparse
import json
from pathlib import Path
import os


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', default=os.getenv('PHASE_DATA_ROOT', 'data'))
    for field in ('project-id', 'system-id', 'cluster-id', 'model-a-id', 'model-b-id'):
        parser.add_argument('--' + field, required=True)
    parser.add_argument('--mode', choices=['parameter', 'ensemble'], default='ensemble')
    parser.add_argument('--reference-a', help='Comma-separated assigned MD sample IDs; omitted = infer from model states')
    parser.add_argument('--reference-b', help='Comma-separated assigned MD sample IDs')
    parser.add_argument('--sample-ids', default='', help='Additional trajectories for endpoint preference visualization')
    parser.add_argument('--block-length', type=int, default=50, help='Contiguous assigned frames per bootstrap block')
    parser.add_argument('--bootstrap', type=int, default=200, help='Replicates; 0 disables uncertainty, otherwise >=50')
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--progress', action='store_true')
    args = vars(parser.parse_args(argv))
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
    try:
        result = run_hotspot_analysis(store=ProjectStore(base_dir=Path(os.environ['PHASE_DATA_ROOT']) / 'projects'), progress_callback=progress, **args)
    finally:
        if bar is not None:
            bar.close()
    print(json.dumps(result, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
