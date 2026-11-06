import os
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


class InstallSolverBindingTests(unittest.TestCase):
    def setUp(self):
        self._tmpdir = tempfile.TemporaryDirectory()
        self.tmpdir = Path(self._tmpdir.name)
        self.fake_bin = self.tmpdir / "bin"
        self.fake_bin.mkdir()
        self.pip_log = self.tmpdir / "pip.log"
        self.mosek_marker = self.tmpdir / "mosek_installed"
        fake_python = self.fake_bin / "python3"
        fake_python.write_text(
            "#!/bin/sh\n"
            "case \"$*\" in\n"
            "  *\"find_spec('bfbarchitect')\"*) exit \"${FAKE_BFBA_EXIT:-1}\" ;;\n"
            "  *\"import gurobipy\"*) exit \"${FAKE_GUROBI_EXIT:-1}\" ;;\n"
            "  *\"import mosek\"*) test -f \"$MOSEK_MARKER\"; exit $? ;;\n"
            "  *\"-m pip install\"*)\n"
            "    echo \"$*\" >> \"$PIP_LOG\"\n"
            "    if [ \"${FAKE_PIP_FIXES_MOSEK:-1}\" = 1 ]; then touch \"$MOSEK_MARKER\"; fi\n"
            "    exit 0 ;;\n"
            "esac\n"
            "exit 0\n"
        )
        fake_python.chmod(0o755)

        self.home = self.tmpdir / "home"
        self.home.mkdir()
        self.data_repo = self.tmpdir / "data_repo"
        self.data_repo.mkdir()

    def tearDown(self):
        self._tmpdir.cleanup()

    def _pip_calls(self):
        return self.pip_log.read_text().splitlines() if self.pip_log.exists() else []

    def _run_finalizer(self, bfb_exit, gurobi_exit, mosek_installed=True, **extra_env):
        if mosek_installed:
            self.mosek_marker.touch()
        env = os.environ.copy()
        env.update({
            "PATH": str(self.fake_bin) + os.pathsep + env.get("PATH", ""),
            "HOME": str(self.home),
            "AA_DATA_REPO": str(self.data_repo),
            "MOSEKLM_LICENSE_FILE": str(self.tmpdir / "missing_mosek"),
            "FAKE_BFBA_EXIT": str(bfb_exit),
            "FAKE_GUROBI_EXIT": str(gurobi_exit),
            "PIP_LOG": str(self.pip_log),
            "MOSEK_MARKER": str(self.mosek_marker),
        })
        env.update(extra_env)
        return subprocess.run(
            ["bash", str(REPO_ROOT / "install.sh"), "--finalize_only"],
            cwd=self.tmpdir,
            env=env,
            capture_output=True,
            text=True,
        )

    def test_finalizer_accepts_installed_gurobi_binding(self):
        result = self._run_finalizer(bfb_exit=0, gurobi_exit=0)

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("Gurobi Python binding is installed", result.stdout)

    def test_finalizer_rejects_missing_gurobi_for_native_install(self):
        result = self._run_finalizer(bfb_exit=0, gurobi_exit=1)

        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn("mamba install -c gurobi gurobi", result.stderr)

    def test_container_host_without_local_bfba_skips_binding_check(self):
        result = self._run_finalizer(bfb_exit=1, gurobi_exit=1, mosek_installed=False)

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("expected when only preparing a host for container use", result.stdout)
        self.assertEqual(self._pip_calls(), [])

    def test_finalizer_installs_missing_mosek(self):
        result = self._run_finalizer(bfb_exit=0, gurobi_exit=0, mosek_installed=False)

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self._pip_calls(), ["-m pip install mosek"])

    def test_finalizer_leaves_working_mosek_alone(self):
        result = self._run_finalizer(bfb_exit=0, gurobi_exit=0)

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self._pip_calls(), [])
        self.assertIn("Mosek Python binding is installed", result.stdout)

    def test_finalizer_warns_when_mosek_cannot_load(self):
        result = self._run_finalizer(bfb_exit=0, gurobi_exit=0, mosek_installed=False, FAKE_PIP_FIXES_MOSEK="0")

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("could not be installed or loaded", result.stderr)
        self.assertIn("will use Clarabel", result.stderr)
        self.assertNotIn("Mosek Python binding is installed", result.stdout)
        self.assertIn("Finished configuring", result.stdout)


if __name__ == "__main__":
    unittest.main()
