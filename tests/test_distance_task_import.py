import os
from pathlib import Path
import subprocess
import sys


def test_rq_distance_task_prefers_source_over_stale_installed_package(tmp_path):
    repository = Path(__file__).resolve().parents[1]
    installed = tmp_path / "site-packages"
    (installed / "phase").mkdir(parents=True)
    (installed / "phase/__init__.py").write_text("OLD_INSTALLED_PACKAGE = True\n")
    code = """
import sys
from pathlib import Path
repository, installed = map(Path, sys.argv[1:])
sys.path = [str(installed)] + [
    value for value in sys.path if Path(value or '.').resolve() != repository
] + [str(repository)]
from rq.utils import import_attribute
task = import_attribute('backend.distance_tasks.run_distance_job')
import phase
assert Path(phase.__file__).resolve() == repository / 'phase/__init__.py'
assert callable(task)
print('RQ distance task imported from current sources')
"""
    result = subprocess.run(
        [sys.executable, "-c", code, str(repository), str(installed)],
        cwd=tmp_path, env={**os.environ, "PYTHONPATH": ""}, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "current sources" in result.stdout
