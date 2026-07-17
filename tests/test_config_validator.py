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

        validate.assert_called_once_with()
        self.assertEqual(args.AA_solver, "mosek")

    def test_missing_mosek_license_selects_clarabel(self):
        args = self._args()
        with patch.object(config_validator, "_validate_mosek_license", return_value=False), \
                self.assertLogs(level=logging.WARNING) as logs:
            self._validate_environment(args)

        self.assertEqual(args.AA_solver, "clarabel")
        self.assertIn("falling back to the Clarabel solver", "\n".join(logs.output))

    def test_missing_license_detection_selects_clarabel(self):
        args = self._args()
        with tempfile.TemporaryDirectory() as tmpdir, patch.dict(os.environ, {
                "AA_SRC": "/aa",
                "AC_SRC": "/ac",
                "HOME": tmpdir,
        }, clear=True), self.assertLogs(level=logging.WARNING):
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


if __name__ == "__main__":
    unittest.main()
