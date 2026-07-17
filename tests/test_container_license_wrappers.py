import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


class ContainerLicenseWrapperTests(unittest.TestCase):
    def setUp(self):
        self._tmpdir = tempfile.TemporaryDirectory()
        self.tmpdir = Path(self._tmpdir.name)

        self.fake_bin = self.tmpdir / "bin"
        self.fake_bin.mkdir()
        self._make_executable("docker", "#!/bin/sh\nexit 0\n")
        self._make_executable(
            "singularity",
            "#!/bin/sh\n"
            "if [ \"$1\" = \"--version\" ]; then\n"
            "    echo \"singularity version 3.11.0\"\n"
            "    exit 0\n"
            "fi\n"
            "exit \"${FAKE_SINGULARITY_EXIT:-0}\"\n",
        )

        self.home = self.tmpdir / "home"
        self.home.mkdir()
        self.data_repo = self.tmpdir / "data_repo"
        (self.data_repo / "GRCh38").mkdir(parents=True)
        (self.data_repo / "coverage.stats").touch()
        self.bam = self.tmpdir / "sample.bam"
        self.bam.write_bytes(b"test bam placeholder")
        self.sif = self.tmpdir / "ampliconsuite.sif"
        self.sif.write_bytes(b"test sif placeholder")
        self.group_input = self.tmpdir / "grouped.input"
        self.group_input.write_text("sample {} tumor\n".format(self.bam))

    def tearDown(self):
        self._tmpdir.cleanup()

    def _make_executable(self, name, contents):
        path = self.fake_bin / name
        path.write_text(contents)
        path.chmod(0o755)

    def _environment(self, **updates):
        env = os.environ.copy()
        env.update({
            "PATH": str(self.fake_bin) + os.pathsep + env.get("PATH", ""),
            "HOME": str(self.home),
            "AA_DATA_REPO": str(self.data_repo),
        })
        env.pop("MOSEKLM_LICENSE_FILE", None)
        env.pop("GRB_LICENSE_FILE", None)
        env.update(updates)
        return env

    def _run(self, relative_script, arguments, env=None):
        result = subprocess.run(
            [sys.executable, str(REPO_ROOT / relative_script)] + arguments,
            cwd=self.tmpdir,
            env=env or self._environment(),
            capture_output=True,
            text=True,
        )
        combined_output = result.stdout + result.stderr
        self.assertEqual(result.returncode, 0, combined_output)
        return combined_output

    def _standalone_arguments(self, output_name):
        return [
            "-o", str(self.tmpdir / output_name),
            "-s", "sample",
            "-t", "1",
            "--ref", "GRCh38",
            "--bam", str(self.bam),
            "--no_cstats",
        ]

    def _singularity_arguments(self, output_name):
        return ["--sif", str(self.sif)] + self._standalone_arguments(output_name)

    def _grouped_arguments(self, output_name):
        return [
            "--sif", str(self.sif),
            "-i", str(self.group_input),
            "-o", str(self.tmpdir / output_name),
            "-t", "1",
            "--ref", "GRCh38",
            "--no_cstats",
        ]

    def _licensed_environment(self):
        mosek_dir = self.tmpdir / "custom_mosek"
        mosek_dir.mkdir()
        (mosek_dir / "mosek.lic").touch()
        gurobi_file = self.tmpdir / "custom_gurobi.lic"
        gurobi_file.touch()
        return self._environment(
            MOSEKLM_LICENSE_FILE=str(mosek_dir),
            GRB_LICENSE_FILE=str(gurobi_file),
        ), mosek_dir, gurobi_file

    def test_docker_explicit_clarabel_has_no_license_fallback_warning(self):
        arguments = self._standalone_arguments("docker_clarabel")
        arguments.extend(["--run_AA", "--run_AC", "--AA_solver", "clarabel"])
        output = self._run("docker/run_paa_docker.py", arguments)

        self.assertNotIn("Mosek is unavailable", output)
        self.assertNotIn("No Gurobi license", output)

    def test_docker_requested_mosek_warns(self):
        arguments = self._standalone_arguments("docker_mosek")
        arguments.extend(["--run_AA", "--AA_solver", "mosek"])
        output = self._run("docker/run_paa_docker.py", arguments)

        self.assertIn("Mosek is unavailable", output)
        self.assertIn("AA also retries with Clarabel", output)

    def test_docker_custom_license_mounts_are_read_only(self):
        env, mosek_dir, gurobi_file = self._licensed_environment()
        arguments = self._standalone_arguments("docker_licensed")
        arguments.extend(["--run_AA", "--run_AC"])
        output = self._run("docker/run_paa_docker.py", arguments, env)

        self.assertIn(str(mosek_dir) + ":/home/mosek/:ro", output)
        self.assertIn(str(gurobi_file) + ":/home/gurobi/gurobi.lic:ro", output)

    def test_singularity_explicit_clarabel_has_no_license_fallback_warning(self):
        arguments = self._singularity_arguments("singularity_clarabel")
        arguments.extend(["--run_AA", "--run_AC", "--AA_solver", "clarabel"])
        output = self._run("singularity/run_paa_singularity.py", arguments)

        self.assertNotIn("Mosek is unavailable", output)
        self.assertNotIn("No Gurobi license", output)

    def test_singularity_custom_license_mounts_are_read_only(self):
        env, mosek_dir, gurobi_file = self._licensed_environment()
        arguments = self._singularity_arguments("singularity_licensed")
        arguments.extend(["--run_AA", "--run_AC"])
        output = self._run("singularity/run_paa_singularity.py", arguments, env)

        self.assertIn(str(mosek_dir) + ":/home/mosek/:ro", output)
        self.assertIn(str(gurobi_file) + ":/home/gurobi/gurobi.lic:ro", output)

    def test_grouped_singularity_requested_mosek_warns(self):
        arguments = self._grouped_arguments("grouped_mosek")
        output = self._run("singularity/run_ga_singularity.py", arguments)

        self.assertIn("Mosek is unavailable", output)
        self.assertNotIn("No Gurobi license", output)

    def test_grouped_singularity_custom_license_mounts_are_read_only(self):
        env, mosek_dir, gurobi_file = self._licensed_environment()
        arguments = self._grouped_arguments("grouped_licensed")
        output = self._run("singularity/run_ga_singularity.py", arguments, env)

        self.assertIn(str(mosek_dir) + ":/home/mosek/:ro", output)
        self.assertIn(str(gurobi_file) + ":/home/gurobi/gurobi.lic:ro", output)

    def test_grouped_no_aa_skips_license_discovery_and_mounts(self):
        env, mosek_dir, gurobi_file = self._licensed_environment()
        arguments = self._grouped_arguments("grouped_no_aa")
        arguments.append("--no_AA")
        output = self._run("singularity/run_ga_singularity.py", arguments, env)

        self.assertNotIn("Mosek is unavailable", output)
        self.assertNotIn(str(mosek_dir), output)
        self.assertNotIn(str(gurobi_file), output)

    def test_grouped_singularity_propagates_container_failure(self):
        arguments = self._grouped_arguments("grouped_failure")
        env = self._environment(FAKE_SINGULARITY_EXIT="9")
        result = subprocess.run(
            [sys.executable, str(REPO_ROOT / "singularity/run_ga_singularity.py")] + arguments,
            cwd=self.tmpdir,
            env=env,
            capture_output=True,
            text=True,
        )

        self.assertEqual(result.returncode, 9, result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
