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
        fake_python = self.fake_bin / "python3"
        fake_python.write_text(
            "#!/bin/sh\n"
            "case \"$*\" in\n"
            "  *\"find_spec('bfbarchitect')\"*) exit \"${FAKE_BFBA_EXIT:-1}\" ;;\n"
            "  *\"import gurobipy\"*) exit \"${FAKE_GUROBI_EXIT:-1}\" ;;\n"
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

    def _run_finalizer(self, bfb_exit, gurobi_exit):
        env = os.environ.copy()
        env.update({
            "PATH": str(self.fake_bin) + os.pathsep + env.get("PATH", ""),
            "HOME": str(self.home),
            "AA_DATA_REPO": str(self.data_repo),
            "MOSEKLM_LICENSE_FILE": str(self.tmpdir / "missing_mosek"),
            "FAKE_BFBA_EXIT": str(bfb_exit),
            "FAKE_GUROBI_EXIT": str(gurobi_exit),
        })
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
        result = self._run_finalizer(bfb_exit=1, gurobi_exit=1)

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("expected when only preparing a host for container use", result.stdout)


if __name__ == "__main__":
    unittest.main()
