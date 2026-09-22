import hashlib
import io
import os
import sys
import tempfile
import subprocess
import unittest
from pathlib import Path
from unittest.mock import patch

import pysam

from paalib import cram_compat


class CramCompatibilityTests(unittest.TestCase):
    def test_version_detection_tolerates_non_utf8_build_details(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            executable = Path(tmpdir) / "samtools"
            executable.write_text("#!/bin/sh\nprintf 'samtools 1.13\\nBuild: \\253\\n'\n")
            executable.chmod(0o755)
            self.assertEqual(cram_compat._samtools_version(str(executable)), "samtools 1.13")

    def test_cram31_rejects_old_or_unknown_samtools_before_fetching(self):
        for version in ("samtools 1.9", "samtools 1.14", "samtools 1.20", "version unknown"):
            with self.subTest(version=version), \
                    patch.object(cram_compat, "read_cram_version", return_value=(3, 1)), \
                    patch.object(cram_compat, "_samtools_version", return_value=version), \
                    patch.object(cram_compat, "fetch_missing_references") as fetch, \
                    patch.object(cram_compat, "_pysam_probe") as probe:
                with self.assertRaisesRegex(cram_compat.CramCompatibilityError, "requires samtools >=1.21"):
                    cram_compat.preflight_cram("sample.cram", "samtools", fetch_refs=True)
                fetch.assert_not_called()
                probe.assert_not_called()

    def test_version_gate_handles_newer_versions_and_leaves_cram30_unchanged(self):
        for cram_version, samtools_version in (((3, 0), "samtools 1.14"),
                                               ((3, 1), "samtools 1.21 (htslib 1.21)"),
                                               ((3, 1), "samtools 1.22.1"),
                                               ((4, 0), "samtools 2.0")):
            with self.subTest(cram=cram_version, samtools=samtools_version), \
                    patch.object(cram_compat, "read_cram_version", return_value=cram_version), \
                    patch.object(cram_compat, "_samtools_version", return_value=samtools_version), \
                    patch.object(cram_compat, "_pysam_probe", return_value="decoder"):
                self.assertEqual(cram_compat.preflight_cram("sample.cram", "samtools"), cram_version)

    def test_future_cram_still_requires_modern_samtools(self):
        with patch.object(cram_compat, "read_cram_version", return_value=(4, 0)), \
                patch.object(cram_compat, "_samtools_version", return_value="samtools 1.20"):
            with self.assertRaisesRegex(cram_compat.CramCompatibilityError, "CRAM 4.0 requires"):
                cram_compat.preflight_cram("sample.cram", "samtools")

    def test_missing_reference_error_gives_fetch_command_without_unconditional_upgrade(self):
        result = subprocess.CompletedProcess([], 1, "pysam=0.23.3 htslib=1.21\n",
                                             "Unable to fetch reference #0\nOSError: truncated file")
        with patch.object(cram_compat, "read_cram_version", return_value=(3, 1)), \
                patch.object(cram_compat, "_samtools_version", return_value="samtools 1.21"), \
                patch.object(cram_compat.subprocess, "run", return_value=result):
            with self.assertRaises(cram_compat.CramCompatibilityError) as caught:
                cram_compat.preflight_cram("a sample.cram", "samtools")
        message = str(caught.exception)
        self.assertIn("Reference lookup failed", message)
        self.assertIn("-m paalib.cram_compat", message)
        self.assertIn("'a sample.cram'", message)
        self.assertIn("If the error indicates an unsupported", message)
        self.assertNotIn("Upgrade the AmpliconSuite environment", message)
        self.assertIn("OSError: truncated file", message)

    def test_other_decoder_errors_do_not_claim_reference_lookup_failed(self):
        result = subprocess.CompletedProcess([], 1, "pysam=0.20.0\n", "Unsupported codec")
        with patch.object(cram_compat, "read_cram_version", return_value=(3, 1)), \
                patch.object(cram_compat, "_samtools_version", return_value="samtools 1.21"), \
                patch.object(cram_compat.subprocess, "run", return_value=result):
            with self.assertRaises(cram_compat.CramCompatibilityError) as caught:
                cram_compat.preflight_cram("sample.cram", "samtools")
        self.assertNotIn("Reference lookup failed", str(caught.exception))
        self.assertIn("Unsupported codec", str(caught.exception))

    def test_reads_cram_version_from_file_definition(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cram = Path(tmpdir) / "sample.cram"
            cram.write_bytes(b"CRAM\x03\x01more-data")
            self.assertEqual(cram_compat.read_cram_version(str(cram)), (3, 1))

    def test_rejects_invalid_cram_file_definition(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cram = Path(tmpdir) / "sample.cram"
            cram.write_bytes(b"not-a-cram")
            with self.assertRaisesRegex(cram_compat.CramCompatibilityError, "valid CRAM"):
                cram_compat.read_cram_version(str(cram))

    def test_bam_skips_cram_preflight(self):
        with patch.object(cram_compat, "_pysam_probe") as probe:
            self.assertIsNone(cram_compat.preflight_cram("sample.bam", "samtools"))
        probe.assert_not_called()

    def test_pysam_probe_decodes_a_real_cram_alignment(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            reference = Path(tmpdir) / "reference.fa"
            reference.write_text(">chr1\n" + "A" * 100 + "\n")
            pysam.faidx(str(reference))

            cram = Path(tmpdir) / "sample.cram"
            header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 100}]}
            read = pysam.AlignedSegment()
            read.query_name = "read1"
            read.query_sequence = "A"
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 0
            read.mapping_quality = 60
            read.cigar = ((0, 1),)
            read.query_qualities = pysam.qualitystring_to_array("I")
            with pysam.AlignmentFile(
                    str(cram), "wc", header=header, reference_filename=str(reference)) as outfile:
                outfile.write(read)

            decoder = cram_compat._pysam_probe(sys.executable, str(cram))
            self.assertIn("pysam=", decoder)
            self.assertIn("htslib=", decoder)

    def test_cram_preflight_decodes_one_alignment(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cram = Path(tmpdir) / "sample.cram"
            cram.write_bytes(b"CRAM\x03\x01")
            with patch.object(cram_compat, "_samtools_version", return_value="samtools 1.21") as samtools_version, \
                    patch.object(cram_compat, "_pysam_probe", return_value="pysam=0.23.3 htslib=1.21") as probe:
                self.assertEqual(
                    cram_compat.preflight_cram(str(cram), "/tools/samtools", "/env/python3"),
                    (3, 1),
                )
        samtools_version.assert_called_once_with("/tools/samtools")
        probe.assert_called_once_with("/env/python3", str(cram))

    def test_cram_preflight_reports_decoder_context(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cram = Path(tmpdir) / "sample.cram"
            cram.write_bytes(b"CRAM\x03\x01")
            failure = cram_compat.CramCompatibilityError("decoder failed")
            with patch.object(cram_compat, "_samtools_version", return_value="samtools 1.21"), \
                    patch.object(cram_compat, "_pysam_probe", side_effect=failure):
                with self.assertRaisesRegex(cram_compat.CramCompatibilityError, "CRAM 3.1 preflight failed"):
                    cram_compat.preflight_cram(str(cram), "samtools")

    def test_expands_flat_and_nested_ref_cache_templates(self):
        digest = "0123456789abcdef0123456789abcdef"
        self.assertEqual(cram_compat._expand_ref_cache("/cache/%s", digest), "/cache/" + digest)
        self.assertEqual(
            cram_compat._expand_ref_cache("/cache/%2s/%2s/%s", digest),
            "/cache/01/23/456789abcdef0123456789abcdef",
        )

    def test_fetch_requires_explicit_ref_cache(self):
        with patch.dict(os.environ, {}, clear=True):
            with self.assertRaisesRegex(cram_compat.CramCompatibilityError, "requires REF_CACHE"):
                cram_compat.fetch_missing_references("sample.cram", "samtools")

    def test_fetch_rejects_cache_without_checksum_placeholder(self):
        with self.assertRaisesRegex(cram_compat.CramCompatibilityError, "checksum placeholder"):
            cram_compat.fetch_missing_references("sample.cram", "samtools", "/cache/reference")

    def test_fetches_and_md5_validates_missing_reference(self):
        sequence = b"ACGT\n"
        digest = hashlib.md5(sequence).hexdigest()
        with tempfile.TemporaryDirectory() as tmpdir, \
                patch.object(cram_compat, "_cram_md5s", return_value=[digest]), \
                patch.object(cram_compat.urllib.request, "urlopen", return_value=io.BytesIO(sequence)):
            template = str(Path(tmpdir) / "%2s" / "%2s" / "%s")
            self.assertEqual(
                cram_compat.fetch_missing_references("sample.cram", "samtools", template),
                1,
            )
            destination = Path(cram_compat._expand_ref_cache(template, digest))
            self.assertEqual(destination.read_bytes(), sequence)


if __name__ == "__main__":
    unittest.main()
