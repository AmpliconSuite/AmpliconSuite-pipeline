import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


REPO_ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("AC_SRC", str(REPO_ROOT))

import GroupedAnalysisAmpSuite as grouped


class GroupedAnalysisACTests(unittest.TestCase):
    def setUp(self):
        self._tmpdir = tempfile.TemporaryDirectory()
        self.tmpdir = Path(self._tmpdir.name)
        self.output_dir = self.tmpdir / "output"
        self.output_dir.mkdir()
        self.ac_src = self.tmpdir / "ac"
        self.ac_src.mkdir()
        self.group_input = self.tmpdir / "cohort.input"
        self.group_input.touch()
        self.sample_lines = []
        self.grouped_seeds = {}

        for sample_name in ("sample_a", "sample_b"):
            sample_dir = self.output_dir / sample_name
            aa_dir = sample_dir / (sample_name + "_AA_results")
            aa_dir.mkdir(parents=True)
            (aa_dir / (sample_name + "_amplicon1_cycles.txt")).touch()
            (aa_dir / (sample_name + "_amplicon1_graph.txt")).touch()
            (aa_dir / (sample_name + "_summary.txt")).touch()
            (sample_dir / (sample_name + "_sample_metadata.json")).write_text(
                json.dumps({"sample_name": sample_name}))
            (sample_dir / (sample_name + "_run_metadata.json")).write_text(
                json.dumps({"AA_version": "test"}))
            seed_path = sample_dir / (sample_name + "_AA_CNV_SEEDS.bed")
            seed_path.touch()
            self.grouped_seeds[sample_name] = str(seed_path)
            self.sample_lines.append([sample_name, str(self.tmpdir / (sample_name + ".bam")), "tumor",
                                      None, None, None])

        self._make_executable("make_input.sh", self._make_input_script())
        self._make_executable("amplicon_classifier.py", self._classifier_script())
        self._make_executable("make_results_table.py", self._results_script())

    def tearDown(self):
        self._tmpdir.cleanup()

    def _make_executable(self, name, contents):
        path = self.ac_src / name
        path.write_text(contents)
        path.chmod(0o755)

    @staticmethod
    def _make_input_script():
        return """#!/usr/bin/env python3
import os
import sys

prefix = sys.argv[-1]
with open(prefix + '.input', 'w') as input_file, open(prefix + '_summary_map.txt', 'w') as summary_file:
    for aa_dir in sys.argv[1:-1]:
        sample = os.path.basename(aa_dir).rsplit('_AA_results', 1)[0]
        cycles = os.path.join(aa_dir, sample + '_amplicon1_cycles.txt')
        graph = os.path.join(aa_dir, sample + '_amplicon1_graph.txt')
        summary = os.path.join(aa_dir, sample + '_summary.txt')
        input_file.write('{}\\t{}\\t{}\\n'.format(sample, cycles, graph))
        summary_file.write('{}\\t{}\\n'.format(sample, summary))
"""

    @staticmethod
    def _classifier_script():
        return """#!/usr/bin/env python3
import os
import sys

if '--help' in sys.argv:
    print('--jobs --bfb_threads --no_results_table')
    sys.exit(0)
if '--version' in sys.argv:
    print('mock-ac-version')
    sys.exit(0)
if os.environ.get('FAIL_GROUPED_AC'):
    sys.exit(7)

input_path = sys.argv[sys.argv.index('-i') + 1]
prefix = sys.argv[sys.argv.index('-o') + 1]
bed_dir = prefix + '_classification_bed_files'
os.makedirs(bed_dir, exist_ok=True)
with open(prefix + '_invocation.txt', 'w') as outfile:
    outfile.write(' '.join(sys.argv[1:]))
with open(prefix + '_amplicon_classification_profiles.tsv', 'w') as profiles, \\
        open(prefix + '_features_to_graph.txt', 'w') as features:
    profiles.write('sample_name\\tamplicon_number\\n')
    with open(input_path) as input_file:
        for line in input_file:
            sample, _, graph = line.rstrip().split('\\t')
            profiles.write('{}\\tamplicon1\\n'.format(sample))
            bed = os.path.join(bed_dir, sample + '_amplicon1_Linear_1_intervals.bed')
            open(bed, 'w').close()
            features.write('{}\\t{}\\n'.format(bed, graph))
with open(prefix + '_feature_similarity_scores.tsv', 'w') as outfile:
    outfile.write('feature1\\tfeature2\\tscore\\n')
"""

    @staticmethod
    def _results_script():
        return """#!/usr/bin/env python3
import sys

classification_file = sys.argv[sys.argv.index('--classification_file') + 1]
prefix = classification_file.rsplit('_amplicon_classification_profiles.tsv', 1)[0]
with open(prefix + '_results_invocation.txt', 'w') as outfile:
    outfile.write(' '.join(sys.argv[1:]))
with open(prefix + '_result_table.tsv', 'w') as outfile:
    outfile.write('Sample name\\n')
"""

    def test_aa_commands_do_not_run_ac_per_sample(self):
        commands = grouped.create_AA_cmds(
            self.sample_lines, " --ref GRCh38", self.grouped_seeds, str(self.output_dir) + os.sep, 2)

        self.assertEqual(set(commands), {"sample_a", "sample_b"})
        for command in commands.values():
            self.assertIn("--run_AA", command)
            self.assertNotIn("--run_AC", command)

    def test_combined_ac_uses_shared_thread_policy_and_updates_metadata(self):
        class_prefix = grouped.run_grouped_ac(
            self.sample_lines,
            self.grouped_seeds,
            str(self.output_dir),
            str(self.group_input),
            "GRCh38",
            8,
            os.environ.get("PYTHON", "python3"),
            str(self.ac_src),
        )

        invocation = Path(class_prefix + "_invocation.txt").read_text()
        self.assertIn("--jobs 2", invocation)
        self.assertIn("--bfb_threads 3", invocation)
        self.assertTrue(Path(class_prefix + "_feature_similarity_scores.tsv").exists())
        self.assertTrue(Path(class_prefix + "_result_table.tsv").exists())

        results_invocation = Path(class_prefix + "_results_invocation.txt").read_text()
        self.assertIn("--sample_metadata_list", results_invocation)
        self.assertIn("--run_metadata_list", results_invocation)
        self.assertIn("--sample_cnv_bed_list", results_invocation)

        for sample_name in ("sample_a", "sample_b"):
            sample_dir = self.output_dir / sample_name
            sample_metadata = json.loads(
                (sample_dir / (sample_name + "_sample_metadata.json")).read_text())
            run_metadata = json.loads(
                (sample_dir / (sample_name + "_run_metadata.json")).read_text())
            self.assertEqual(sample_metadata["number_of_AA_amplicons"], 1)
            self.assertEqual(sample_metadata["number_of_AA_features"], 1)
            self.assertEqual(run_metadata["AC_version"], "mock-ac-version")
            self.assertIn("--jobs 2", run_metadata["AC_cmd"])

    def test_combined_ac_failure_is_reported(self):
        with patch.dict(os.environ, {"FAIL_GROUPED_AC": "1"}):
            with self.assertRaisesRegex(RuntimeError, "Combined AmpliconClassifier stage failed"):
                grouped.run_grouped_ac(
                    self.sample_lines,
                    self.grouped_seeds,
                    str(self.output_dir),
                    str(self.group_input),
                    "GRCh38",
                    4,
                    os.environ.get("PYTHON", "python3"),
                    str(self.ac_src),
                )


if __name__ == "__main__":
    unittest.main()
