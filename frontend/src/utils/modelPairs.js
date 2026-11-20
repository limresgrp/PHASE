// Selection order fixes the sign: earlier model is A, later model is B.
export function modelPairs(ids) {
  const unique = [...new Set(ids.filter(Boolean))];
  if (unique.length < 2) throw new Error('Select at least two different models.');
  return unique.flatMap((a, i) => unique.slice(i + 1).map((b) => [a, b]));
}

export async function submitModelPairs(pairs, submit) {
  const jobs = [], errors = [];
  for (const [a, b] of pairs) {
    try { jobs.push({ ...await submit(a, b), model_a_id: a, model_b_id: b }); }
    catch (error) { errors.push(`${a} vs ${b}: ${error.message}`); }
  }
  return { jobs, errors };
}

export async function fetchPairJobStatus(jobs, fetchStatus) {
  const statuses = await Promise.all(jobs.map((job) => fetchStatus(job.job_id)));
  if (statuses.length === 1) return statuses[0];
  const terminal = (s) => ['finished', 'failed', 'stopped', 'canceled'].includes(s.status);
  const done = statuses.filter(terminal).length;
  const failed = statuses.filter((s) => terminal(s) && s.status !== 'finished').length;
  const allDone = done === statuses.length;
  const progress = statuses.reduce((sum, s) => sum + (terminal(s) ? 100 : Math.min(100, Math.max(0, Number(s.meta?.progress) || 0))), 0) / statuses.length;
  return { status: allDone ? (failed ? 'failed' : 'finished') : 'started',
    meta: { progress: Math.round(progress), status: `${done}/${statuses.length} pairs completed; ${failed} failed`,
      error: failed ? statuses.filter((s) => terminal(s) && s.status !== 'finished').map((s, i) => s.meta?.error || s.result?.error || `Pair job failed (${i + 1})`).join('; ') : undefined },
    jobs: jobs.map((job, i) => ({ ...job, status: statuses[i].status })) };
}
