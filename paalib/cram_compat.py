"""Small, CRAM-only compatibility checks used before alignment processing."""

import hashlib
import logging
import os
import re
import shlex
import subprocess
import sys
import tempfile
import urllib.error
import urllib.request


CRAM_MAGIC = b"CRAM"
ENA_CRAM_MD5_URL = "https://www.ebi.ac.uk/ena/cram/md5/{}"


class CramCompatibilityError(RuntimeError):
    """Raised when a CRAM cannot be decoded by a tool used by the pipeline."""


def read_cram_version(path):
    """Return the ``(major, minor)`` version from the CRAM file definition."""
    with open(path, "rb") as infile:
        file_definition = infile.read(6)

    if len(file_definition) != 6 or file_definition[:4] != CRAM_MAGIC:
        raise CramCompatibilityError(
            "{} has a .cram suffix but does not begin with a valid CRAM file definition.".format(path)
        )
    return file_definition[4], file_definition[5]


def _samtools_version(samtools_path):
    try:
        result = subprocess.run(
            [samtools_path, "--version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            errors="replace",
            timeout=30,
        )
    except (OSError, subprocess.TimeoutExpired):
        return "version unknown"
    if result.returncode == 0 and result.stdout.splitlines():
        return result.stdout.splitlines()[0]
    return "version unknown"


def _pysam_probe(python_path, cram_path):
    """Decode one alignment using the same Python interpreter used by CNVKit."""
    probe = r'''
import pysam
import sys

print("pysam={} htslib={}".format(
    getattr(pysam, "__version__", "unknown"),
    getattr(pysam, "__htslib_version__", getattr(pysam, "__samtools_version__", "unknown"))))
with pysam.AlignmentFile(sys.argv[1], "rc") as alignment_file:
    next(alignment_file.fetch(until_eof=True), None)
'''
    try:
        result = subprocess.run(
            [python_path, "-c", probe, cram_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            timeout=120,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise CramCompatibilityError("Could not run the pysam CRAM preflight: {}".format(exc))

    version_line = result.stdout.splitlines()[0] if result.stdout.splitlines() else "pysam/HTSlib version unknown"
    if result.returncode != 0:
        detail = result.stderr.strip() or "pysam exited with status {}".format(result.returncode)
        raise CramCompatibilityError(
            "CNVKit's CRAM decoder failed while reading the first alignment ({}).\n"
            "CNVKit uses pysam from '{}'; --samtools_path does not replace this embedded HTSlib.\n"
            "Check the decoder error below: unavailable references, unsupported codecs, or damaged input can cause this failure.\n"
            "If the error indicates an unsupported format or codec, update the AmpliconSuite environment "
            "(recommended: pysam>=0.23.3).\n"
            "Decoder error: {}".format(version_line, python_path, detail)
        )
    return version_line


def _cram_md5s(samtools_path, cram_path):
    """Read unique reference MD5s from @SQ records without decoding alignments."""
    try:
        result = subprocess.run(
            [samtools_path, "view", "-H", cram_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            timeout=120,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise CramCompatibilityError("Could not read the CRAM header with samtools: {}".format(exc))

    if result.returncode != 0:
        raise CramCompatibilityError(
            "samtools could not read the CRAM header: {}".format(result.stderr.strip())
        )

    md5s = []
    seen = set()
    for line in result.stdout.splitlines():
        if not line.startswith("@SQ\t"):
            continue
        match = re.search(r"(?:^|\t)M5:([0-9a-fA-F]{32})(?:\t|$)", line)
        if match:
            digest = match.group(1).lower()
            if digest not in seen:
                seen.add(digest)
                md5s.append(digest)
    return md5s


def _expand_ref_cache(template, digest):
    """Expand HTSlib's sequential %s / %Ns REF_CACHE placeholders."""
    offset = 0

    def replace(match):
        nonlocal offset
        width_text = match.group(1)
        if width_text:
            width = int(width_text)
            value = digest[offset:offset + width]
            offset += width
            return value
        value = digest[offset:]
        offset = len(digest)
        return value

    return re.sub(r"%(\d*)s", replace, os.path.expanduser(template))


def fetch_missing_references(cram_path, samtools_path, ref_cache=None, ref_paths=()):
    """Fetch absent CRAM references from ENA into REF_CACHE, without prompting."""
    cache_template = ref_cache or os.environ.get("REF_CACHE")
    if not cache_template:
        raise CramCompatibilityError(
            "--fetch_cram_refs requires REF_CACHE to be set, for example "
            "REF_CACHE=$HOME/.cache/samtools/%s."
        )
    if not re.search(r"%\d*s", cache_template):
        raise CramCompatibilityError(
            "REF_CACHE must contain an HTSlib checksum placeholder such as %s or %2s/%2s/%s."
        )

    md5s = _cram_md5s(samtools_path, cram_path)
    if not md5s:
        logging.info("No M5 reference checksums were present in the CRAM header; nothing to fetch.")
        return 0

    fetched = 0
    for digest in md5s:
        destination = _expand_ref_cache(cache_template, digest)
        if os.path.isfile(destination):
            continue
        if any(os.path.isfile(_expand_ref_cache(template, digest)) for template in ref_paths):
            continue

        destination_dir = os.path.dirname(destination) or "."
        os.makedirs(destination_dir, exist_ok=True)
        url = ENA_CRAM_MD5_URL.format(digest)
        logging.info("Fetching missing CRAM reference {} from ENA".format(digest))
        tmp_name = None
        try:
            with urllib.request.urlopen(url, timeout=120) as response, \
                    tempfile.NamedTemporaryFile(dir=destination_dir, delete=False) as tmpfile:
                tmp_name = tmpfile.name
                observed = hashlib.md5()
                while True:
                    block = response.read(1024 * 1024)
                    if not block:
                        break
                    observed.update(block)
                    tmpfile.write(block)

            if observed.hexdigest() != digest:
                raise CramCompatibilityError(
                    "Reference {} downloaded from {} had MD5 {}.".format(
                        digest, url, observed.hexdigest())
                )
            os.replace(tmp_name, destination)
            tmp_name = None
            fetched += 1
        except (OSError, urllib.error.URLError) as exc:
            raise CramCompatibilityError(
                "Failed to download CRAM reference {} from {}: {}".format(digest, url, exc)
            )
        finally:
            if tmp_name and os.path.exists(tmp_name):
                os.unlink(tmp_name)
    return fetched


def preflight_cram(cram_path, samtools_path, python_path="python3", fetch_refs=False, ref_paths=()):
    """Run fast CRAM checks; do nothing for BAM and other input formats."""
    if not cram_path or not cram_path.lower().endswith(".cram"):
        return None

    version = read_cram_version(cram_path)
    logging.info("Detected CRAM {}.{} input".format(*version))

    samtools_decoder = _samtools_version(samtools_path)
    match = re.match(r"^samtools\s+(\d+)\.(\d+)(?:\D|$)", samtools_decoder)
    if version >= (3, 1) and (not match or tuple(map(int, match.groups())) < (1, 21)):
        raise CramCompatibilityError(
            "CRAM {}.{} requires samtools >=1.21 in AmpliconSuite; '{}' reports '{}'. "
            "Install a supported samtools and select it with --samtools_path. "
            "A missing or unrecognized version cannot establish compatibility.".format(
                version[0], version[1], samtools_path, samtools_decoder)
        )

    if fetch_refs:
        fetched = fetch_missing_references(cram_path, samtools_path, ref_paths=ref_paths)
        logging.info("Fetched {} missing CRAM reference sequence(s)".format(fetched))

    try:
        decoder = _pysam_probe(python_path, cram_path)
    except CramCompatibilityError as exc:
        message = str(exc)
        if any(marker in message.lower() for marker in (
                "failed to populate reference", "unable to fetch reference",
                "failed to open reference", "failed to load reference")):
            package_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            command = " ".join(shlex.quote(arg) for arg in (
                "env", "PYTHONPATH=" + package_root, python_path, "-m", "paalib.cram_compat",
                "--samtools_path", samtools_path,
                "--fetch_cram_refs", cram_path))
            message += (
                "\nReference lookup failed. Make the matching sequences available through REF_PATH/REF_CACHE, "
                "or download missing sequences and retry:\n"
                'export REF_CACHE="${REF_CACHE:-$HOME/.cache/hts-ref/%2s/%2s/%s}"\n' + command
            )
        raise CramCompatibilityError(
            "CRAM {}.{} preflight failed for '{}'.\nExternal decoder: {}.\n{}".format(
                version[0], version[1], cram_path, samtools_decoder, message
            )
        )

    logging.info(
        "CRAM preflight decoded one alignment successfully with CNVKit's decoder ({}); external {}".format(
            decoder, samtools_decoder
        )
    )
    return version


def main(argv=None):
    """Check CRAM inputs, optionally populating a persistent reference cache."""
    from paalib.argument_parser import setup_cram_argument_parser

    parser = setup_cram_argument_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    try:
        for path in args.crams:
            if not path.lower().endswith(".cram"):
                parser.error("expected a .cram input: {}".format(path))
            ref_paths = os.environ.get("REF_PATH", "").split(os.pathsep) if args.use_ref_path else ()
            preflight_cram(path, args.samtools_path, sys.executable, args.fetch_cram_refs, ref_paths)
    except (CramCompatibilityError, OSError) as exc:
        logging.error(str(exc))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
