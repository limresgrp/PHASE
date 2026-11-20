import { modelPairs, submitModelPairs, fetchPairJobStatus } from './modelPairs';

test('three models generate three oriented pairs without duplicates', () => {
  expect(modelPairs(['a', 'b', 'a', 'c'])).toEqual([['a', 'b'], ['a', 'c'], ['b', 'c']]);
  expect(() => modelPairs(['a', 'a'])).toThrow(/at least two/);
});

test('submits independent jobs and continues after submission failures', async () => {
  const submit = jest.fn(async (a, b) => {
    if (a === 'a' && b === 'c') throw new Error('invalid references');
    return { job_id: `${a}-${b}` };
  });
  const result = await submitModelPairs(modelPairs(['a', 'b', 'c']), submit);
  expect(submit).toHaveBeenCalledTimes(3);
  expect(result.jobs.map((j) => j.job_id)).toEqual(['a-b', 'b-c']);
  expect(result.errors).toEqual(['a vs c: invalid references']);
});

test('one failed job does not stop polling queued or running pairs', async () => {
  const jobs = [{ job_id: 'ab' }, { job_id: 'ac' }, { job_id: 'bc' }];
  const fetch = jest.fn().mockResolvedValueOnce({ status: 'failed', meta: { error: 'bad pair' } })
    .mockResolvedValueOnce({ status: 'finished' }).mockResolvedValueOnce({ status: 'queued' });
  const result = await fetchPairJobStatus(jobs, fetch);
  expect(result.status).toBe('started');
  expect(result.meta.progress).toBe(67);
  expect(result.meta.status).toBe('2/3 pairs completed; 1 failed');
  fetch.mockResolvedValue({ status: 'finished' });
  expect((await fetchPairJobStatus(jobs, fetch)).status).toBe('finished');
});
