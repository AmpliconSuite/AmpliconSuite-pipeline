import logging
import os
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from paalib import config_validator


class AASolverLicenseValidationTests(unittest.TestCase):
    @staticmethod
    def _args(solver="mosek", run_aa=True):
        return SimpleNamespace(run_AA=run_aa, AA_solver=solver, AA_src=None)

    def _validate_environment(self, args):
        with patch.dict(os.environ, {"AA_SRC": "/aa", "AC_SRC": "/ac"}, clear=True):
            return config_validator.validate_aa_environment(args)

    def test_explicit_clarabel_skips_mosek_validation(self):
        args = self._args(solver="clarabel")
        with patch.object(config_validator, "_validate_mosek_license") as validate:
            self.assertEqual(self._validate_environment(args), ("/aa", "/ac"))

        validate.assert_not_called()
        self.assertEqual(args.AA_solver, "clarabel")

    def test_available_mosek_license_preserves_mosek(self):
        args = self._args()
        with patch.object(config_validator, "_validate_mosek_license", return_value=True) as validate:
            self._validate_environment(args)

        validate.assert_called_once_with("python3")
        self.assertEqual(args.AA_solver, "mosek")

    def test_missing_mosek_license_selects_clarabel(self):
        args = self._args()
        with patch.object(config_validator, "_validate_mosek_license", return_value=False), \
                self.assertLogs(level=logging.INFO) as logs:
            self._validate_environment(args)

        self.assertEqual(args.AA_solver, "clarabel")
        self.assertIn("AA will use the Clarabel solver", "\n".join(logs.output))
        self.assertFalse([record for record in logs.records if record.levelno >= logging.WARNING])

    def test_missing_license_detection_selects_clarabel(self):
        args = self._args()
        with tempfile.TemporaryDirectory() as tmpdir, patch.dict(os.environ, {
                "AA_SRC": "/aa",
                "AC_SRC": "/ac",
                "HOME": tmpdir,
        }, clear=True), self.assertLogs(level=logging.INFO):
            config_validator.validate_aa_environment(args)

        self.assertEqual(args.AA_solver, "clarabel")

    def test_no_aa_run_skips_mosek_validation(self):
        args = self._args(run_aa=False)
        with patch.object(config_validator, "_validate_mosek_license") as validate:
            self._validate_environment(args)

        validate.assert_not_called()

    def test_missing_license_returns_false(self):
        with tempfile.TemporaryDirectory() as tmpdir, \
                patch.dict(os.environ, {"HOME": tmpdir}, clear=True):
            self.assertFalse(config_validator._validate_mosek_license())

    def test_present_license_returns_true(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            license_dir = Path(tmpdir) / "mosek"
            license_dir.mkdir()
            (license_dir / "mosek.lic").touch()
            with patch.dict(os.environ, {
                    "HOME": str(Path(tmpdir) / "other-home"),
                    "MOSEKLM_LICENSE_FILE": str(license_dir),
            }, clear=True):
                self.assertTrue(config_validator._validate_mosek_license())

    def _validate_with_license(self, license_text, installed_major):
        with tempfile.TemporaryDirectory() as tmpdir:
            license_dir = Path(tmpdir) / "mosek"
            license_dir.mkdir()
            (license_dir / "mosek.lic").write_text(license_text)
            with patch.dict(os.environ, {"HOME": tmpdir}, clear=True), \
                    patch.object(config_validator, "_installed_mosek_major", return_value=installed_major):
                return config_validator._validate_mosek_license()

    LICENSE_V11 = (
        "FEATURE PTS MOSEKLM 11 20-dec-2026 uncounted VENDOR_STRING=x \\\n"
        "    HOSTID=ANY SIGN=\"abc\"\n"
        "FEATURE PTON MOSEKLM 11 20-dec-2026 uncounted SIGN=\"abc\"\n"
    )

    def test_newer_mosek_than_license_returns_false(self):
        self.assertFalse(self._validate_with_license(self.LICENSE_V11, installed_major=12))

    def test_same_major_as_license_returns_true(self):
        self.assertTrue(self._validate_with_license(self.LICENSE_V11, installed_major=11))

    def test_older_mosek_than_license_returns_true(self):
        self.assertTrue(self._validate_with_license(self.LICENSE_V11, installed_major=10))

    def test_undetermined_mosek_version_returns_true(self):
        self.assertTrue(self._validate_with_license(self.LICENSE_V11, installed_major=None))

    def test_license_major_uses_lowest_required_feature(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            license_path = Path(tmpdir) / "mosek.lic"
            license_path.write_text("FEATURE PTS MOSEKLM 12 x\nFEATURE PTON MOSEKLM 11 x\nFEATURE OTHER MOSEKLM 9 x\n")
            self.assertEqual(config_validator._mosek_license_major(license_path), 11)

    def test_license_major_unstated_returns_none(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            license_path = Path(tmpdir) / "mosek.lic"
            license_path.write_text("USE_SERVER\n")
            self.assertIsNone(config_validator._mosek_license_major(license_path))
            self.assertIsNone(config_validator._mosek_license_major(Path(tmpdir) / "missing.lic"))


if __name__ == "__main__":
    unittest.main()
