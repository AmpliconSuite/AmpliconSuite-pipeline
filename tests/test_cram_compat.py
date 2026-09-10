import hashlib
import io
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pysam

from paalib import cram_compat


class CramCompatibilityTests(unittest.TestCase):
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
