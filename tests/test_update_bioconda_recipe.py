import importlib.util
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
UPDATER_PATH = REPO_ROOT / "conda-recipe" / "update_bioconda_recipe.py"
SPEC = importlib.util.spec_from_file_location("update_bioconda_recipe", UPDATER_PATH)
UPDATER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(UPDATER)


class UpdateBiocondaRecipeTests(unittest.TestCase):
    def test_updates_all_four_versions_hashes_and_build_number(self):
        original = (REPO_ROOT / "conda-recipe" / "meta.yaml").read_text()
        versions = {
            "AS": "9.1.0",
            "AA": "9.2.r0",
            "AC": "9.3.0",
            "BFBA": "9.4.0",
        }
        hashes = {
            "AS": "a" * 64,
            "AA": "b" * 64,
            "AC": "c" * 64,
            "BFBA": "d" * 64,
        }

        updated = UPDATER.update_recipe(original, versions, hashes, build_number=7)

        for short_name, variable in UPDATER.COMPONENTS:
            self.assertIn(
                '{{% set {}="{}" %}}'.format(variable, versions[short_name]),
                updated,
            )
        for digest in hashes.values():
            self.assertIn("sha256: " + digest, updated)
        self.assertIn("    number: 7", updated)
        self.assertIn("folder: bfbarchitectlib", updated)
        self.assertIn("bfba-conda-solvers.patch", updated)


if __name__ == "__main__":
    unittest.main()
