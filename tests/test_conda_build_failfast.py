import os
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


class CondaBuildFailFastTests(unittest.TestCase):
    def _run_build(self, bfb_project, ac_project):
        with tempfile.TemporaryDirectory() as tmpdir:
            workdir = Path(tmpdir)
            (workdir / "bfbarchitectlib").mkdir()
            (workdir / "ampliconclassifierlib").mkdir()
            if bfb_project:
                (workdir / "bfbarchitectlib" / "pyproject.toml").touch()
            if ac_project:
                (workdir / "ampliconclassifierlib" / "pyproject.toml").touch()
            prefix = workdir / "prefix"
            prefix.mkdir()
            marker = workdir / "python_was_called"
            fake_python = workdir / "fake_python"
            fake_python.write_text(
                "#!/bin/sh\n"
                "touch \"$FAKE_PYTHON_MARKER\"\n"
                "exit 0\n"
            )
            fake_python.chmod(0o755)
            env = os.environ.copy()
            env.update({
                "PREFIX": str(prefix),
                "PYTHON": str(fake_python),
                "FAKE_PYTHON_MARKER": str(marker),
            })

            result = subprocess.run(
                ["bash", str(REPO_ROOT / "conda-recipe" / "build.sh")],
                cwd=workdir,
                env=env,
                capture_output=True,
                text=True,
            )

            return result, marker.exists()

    def test_missing_bfba_source_fails_before_any_package_install(self):
        result, python_was_called = self._run_build(bfb_project=False, ac_project=True)

        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn("BFBArchitect source is not independently installable", result.stderr)
        self.assertFalse(python_was_called, "The build attempted a partial package installation")

    def test_legacy_ac_source_fails_before_any_package_install(self):
        result, python_was_called = self._run_build(bfb_project=True, ac_project=False)

        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn("does not support the four-tool package layout", result.stderr)
        self.assertFalse(python_was_called, "The build attempted a partial package installation")


if __name__ == "__main__":
    unittest.main()
