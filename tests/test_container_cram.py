import io
import os
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from paalib import container_cram, cram_compat


class ContainerCramTests(unittest.TestCase):
    def test_bam_does_not_launch_preparation_or_create_cache(self):
        with patch.object(container_cram.subprocess, "run") as run, \
                patch.object(container_cram.os, "makedirs") as mkdir:
            self.assertEqual(container_cram.prepare_cram_container("docker", "image", ["sample.bam", None]), [])
        run.assert_not_called()
        mkdir.assert_not_called()

    def test_cache_template_mount_preserves_nested_and_prefixed_layouts(self):
        for template, root, target in (("/cache/%2s/%2s/%s", "/cache", "/mounted/%2s/%2s/%s"),
                                       ("/cache/ref-%s", "/cache", "/mounted/ref-%s")):
            self.assertEqual(container_cram.cache_mount(template, "/mounted"), (root, target))
        with self.assertRaises(cram_compat.CramCompatibilityError):
            container_cram.cache_mount("/cache/no-placeholder", "/mounted")
        with self.assertRaises(cram_compat.CramCompatibilityError):
            container_cram.cache_mount("/%s", "/mounted")

    def test_prepare_uses_container_tools_host_uid_and_readonly_analysis_cache(self):
        with tempfile.TemporaryDirectory(prefix="cram cache '") as directory:
            root = Path(directory)
            cram = root / "normal.cram"
            cram.touch()
            refs = root / "other refs"
            refs.mkdir()
            cache = root / "cache"
            with patch.dict(os.environ, {"REF_CACHE": str(cache / "%2s/%2s/%s"),
                                         "REF_PATH": str(refs / "%s")}), \
                    patch.object(container_cram.subprocess, "run") as run:
                result = container_cram.prepare_cram_container("docker", "test:dev", ["tumor.bam", str(cram), str(cram)])
            command = run.call_args[0][0]
            self.assertIn(str(cache) + ":/home/cram_cache", command)
            self.assertIn(str(refs) + ":/home/cram_ref_1:ro", command)
            self.assertIn("--user", command)
            self.assertIn("--fetch_cram_refs", command)
            self.assertIn("--use_ref_path", command)
            self.assertEqual(command.count("/home/cram_input_0.cram"), 1)
            self.assertNotIn("/home/cram_input_1.cram", command)
            self.assertIn(str(cache) + ":/home/cram_cache:ro", result)
            self.assertIn("REF_PATH=/home/cram_cache/%2s/%2s/%s:/home/cram_ref_1/%s", result)
            self.assertTrue(run.call_args.kwargs["check"])

    def test_singularity_default_cache_survives_no_home_cleanenv(self):
        with tempfile.TemporaryDirectory() as directory:
            cram = Path(directory) / "input.cram"
            cram.touch()
            with patch.dict(os.environ, {"HOME": directory}, clear=True), \
                    patch.object(container_cram.subprocess, "run") as run:
                result = container_cram.prepare_cram_container("singularity", "image.sif", [str(cram)])
            self.assertTrue((Path(directory) / ".cache/hts-ref").is_dir())
            self.assertIn("--cleanenv", run.call_args[0][0])
            self.assertIn("--no-home", run.call_args[0][0])
            self.assertIn("--env", result)
            self.assertIn("REF_CACHE=/home/cram_cache/%2s/%2s/%s", result)

    def test_preparation_failure_stops_launcher(self):
        with tempfile.TemporaryDirectory() as directory:
            cram = Path(directory) / "input.cram"
            cram.touch()
            with patch.dict(os.environ, {"REF_CACHE": directory + "/cache/%s", "REF_PATH": ""}), \
                    patch.object(container_cram.subprocess, "run",
                                 side_effect=subprocess.CalledProcessError(1, "container")):
                with self.assertRaisesRegex(cram_compat.CramCompatibilityError, "analysis was not started"):
                    container_cram.prepare_cram_container("singularity", "image.sif", [str(cram)])

    def test_fetch_reuses_secondary_checksum_store_and_downloads_only_missing(self):
        import hashlib
        sequence = b"ACGT"
        missing = hashlib.md5(sequence).hexdigest()
        present = hashlib.md5(b"AAAA").hexdigest()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "refs").mkdir()
            (root / "refs" / present).write_bytes(b"AAAA")
            with patch.object(cram_compat, "_cram_md5s", return_value=[present, missing]), \
                    patch.object(cram_compat.urllib.request, "urlopen", return_value=io.BytesIO(sequence)) as fetch:
                count = cram_compat.fetch_missing_references("input.cram", "samtools", str(root / "cache/%s"),
                                                            [str(root / "refs/%s")])
            self.assertEqual(count, 1)
            fetch.assert_called_once()
            self.assertFalse((root / "cache" / present).exists())
            self.assertEqual((root / "cache" / missing).read_bytes(), sequence)
