import { useEffect, useMemo, useState } from 'react';
import { useNavigate, useParams } from 'react-router-dom';
import Loader from '../components/common/Loader';
import ErrorMessage from '../components/common/ErrorMessage';
import { fetchResult } from '../api/jobs';

export default function SimulationResultPage() {
  const { jobId } = useParams();
  const navigate = useNavigate();
  const [result, setResult] = useState(null);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState(null);

  useEffect(() => {
    const load = async () => {
      setIsLoading(true);
      setError(null);
      try {
        const data = await fetchResult(jobId);
        if (data.analysis_type !== 'simulation') {
          throw new Error('This result is not a sampling job.');
        }
        setResult(data);
      } catch (err) {
        setError(err.message || 'Failed to load sampling result.');
      } finally {
        setIsLoading(false);
      }
    };
    load();
  }, [jobId]);

  const systemRef = result?.system_reference;
  const samplingUrl = useMemo(() => {
    if (!systemRef?.project_id || !systemRef.system_id || !systemRef.cluster_id) return null;
    const params = new URLSearchParams({ cluster_id: systemRef.cluster_id });
    if (systemRef.sample_id) params.set('sample_id', systemRef.sample_id);
    return `/projects/${systemRef.project_id}/systems/${systemRef.system_id}/sampling/visualize?${params}`;
  }, [systemRef]);

  useEffect(() => {
    if (samplingUrl) {
      navigate(samplingUrl, { replace: true });
    }
  }, [samplingUrl, navigate]);

  if (isLoading) return <Loader message="Loading sampling result..." />;
  if (error) return <ErrorMessage message={error} />;
  if (!result) return null;

  return (
    <div className="space-y-4">
      <h1 className="text-2xl font-bold text-white">Sampling results</h1>
      <p className="text-sm text-gray-400">
        Sampling reports now open in the Sampling Explorer. Use the button below if the redirect does not happen.
      </p>
      <button
        type="button"
        onClick={() => samplingUrl && navigate(samplingUrl)}
        className="px-3 py-2 rounded-md bg-cyan-600 text-white text-sm disabled:opacity-60"
        disabled={!samplingUrl}
      >
        Open Sampling Explorer
      </button>
    </div>
  );
}
