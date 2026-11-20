"""Small batch driver; each pair retains the ordinary independent analysis format."""
from itertools import combinations


def model_pairs(model_ids=None, model_a_id=None, model_b_id=None):
    if model_ids:
        if model_a_id or model_b_id:
            raise ValueError("Use --model-ids OR --model-a-id/--model-b-id, not both.")
        ids = list(dict.fromkeys(v.strip() for v in model_ids if v.strip()))
    else:
        ids = list(dict.fromkeys(v for v in (model_a_id, model_b_id) if v))
    if len(ids) < 2:
        raise ValueError("Select at least two different models.")
    return list(combinations(ids, 2))


def run_model_pairs(pairs, run, report=print):
    """Attempt every pair, reporting failures without discarding successful results."""
    failures = []
    for index, (a, b) in enumerate(pairs, 1):
        if len(pairs) > 1:
            report(f"[model_pairs] Running pair {index}/{len(pairs)}: {a} vs {b}")
        try:
            run(a, b)
        except Exception as exc:
            if len(pairs) == 1:
                raise
            failures.append((a, b, str(exc)))
            report(f"[model_pairs] FAILED {a} vs {b}: {exc}")
    if len(pairs) > 1:
        report(f"[model_pairs] {len(pairs) - len(failures)}/{len(pairs)} pairs succeeded; {len(failures)} failed.")
    return 1 if failures else 0
