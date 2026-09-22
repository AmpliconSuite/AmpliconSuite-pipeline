"""Prepare persistent CRAM references for the container launchers."""

import os
import re
import shlex
import subprocess

from paalib.cram_compat import CramCompatibilityError


def cache_mount(template, destination):
    """Map a local checksum template into a container without mounting all of HOME."""
    template = os.path.abspath(os.path.expanduser(template))
    match = re.search(r"%\d*s", template)
    if not match:
        raise CramCompatibilityError("REF_CACHE must contain a checksum placeholder such as %s or %2s/%2s/%s.")
    root = os.path.dirname(template[:match.start()])
    if root == os.path.sep:
        raise CramCompatibilityError("Place CRAM references in a dedicated cache directory, not the filesystem root.")
    relative_template = os.path.relpath(template, root)
    root = os.path.realpath(root)
    if any(char in root for char in (":", ",", "\n")):
        raise CramCompatibilityError("Container reference cache paths cannot contain colons, commas or newlines.")
    return root, destination + "/" + relative_template


def prepare_cram_container(runtime, image, alignment_paths):
    """Fetch missing references as the host user; return read-only analysis mounts.

    Header inspection and downloads run in the selected image, using its samtools
    and the shared fetcher. BAM-only launches do not create or inspect a cache.
    """
    crams = list(dict.fromkeys(os.path.realpath(path) for path in alignment_paths
                              if path and path.lower().endswith(".cram")))
    if not crams:
        return []
    template = os.environ.get("REF_CACHE") or "~/.cache/hts-ref/%2s/%2s/%s"
    cache_root, container_template = cache_mount(template, "/home/cram_cache")
    os.makedirs(cache_root, exist_ok=True)
    mounts = [(cache_root, "/home/cram_cache")]
    ref_paths = [container_template]
    # Reuse additional local checksum stores advertised through REF_PATH.
    # Remote URLs and FASTA paths are not checksum stores; ENA fills any gaps.
    for entry in os.environ.get("REF_PATH", "").split(os.pathsep):
        if not entry.startswith(("/", "~")) or not re.search(r"%\d*s", entry):
            continue
        if os.path.abspath(os.path.expanduser(entry)) == os.path.abspath(os.path.expanduser(template)):
            continue
        source, target = cache_mount(entry, "/home/cram_ref_{}".format(len(mounts)))
        if os.path.isdir(source):
            # Keep the destination independent of the template's nested layout.
            mount_point = "/home/cram_ref_{}".format(len(mounts))
            mounts.append((source, mount_point))
            ref_paths.append(target)

    def options(read_only):
        result = []
        for i, (source, target) in enumerate(mounts):
            mode = ":ro" if read_only or i else ""
            result.extend(["-v" if runtime == "docker" else "--bind", source + ":" + target + mode])
        for name, value in (("REF_CACHE", container_template), ("REF_PATH", ":".join(ref_paths))):
            result.extend(["-e" if runtime == "docker" else "--env", name + "=" + value])
        return result

    uid_options = ["--user", "{}:{}".format(os.getuid(), os.getgid())] if runtime == "docker" else []
    input_options = []
    container_inputs = []
    for i, path in enumerate(crams):
        if not os.path.isfile(path):
            raise CramCompatibilityError("CRAM input not found: {}".format(path))
        target = "/home/cram_input_{}.cram".format(i)
        input_options.extend(["-v" if runtime == "docker" else "--bind", path + ":" + target + ":ro"])
        container_inputs.append(target)
    package_dir = "/home/programs/AmpliconSuite-pipeline-master"
    command = (["docker", "run", "--rm"] + uid_options if runtime == "docker" else
               ["singularity", "exec", "--cleanenv", "--no-home"])
    command += options(False) + input_options
    command += (["-w", package_dir, "--entrypoint", "python3", image] if runtime == "docker" else
                ["--pwd", package_dir, image, "python3"])
    command += ["-m", "paalib.cram_compat", "--fetch_cram_refs", "--use_ref_path"] + container_inputs
    print("Preparing CRAM references in persistent cache {}".format(cache_root), flush=True)
    print(" ".join(shlex.quote(arg) for arg in command), flush=True)
    try:
        subprocess.run(command, check=True)
    except (OSError, subprocess.CalledProcessError) as exc:
        raise CramCompatibilityError(
            "CRAM reference preparation failed; analysis was not started. "
            "Use an image with the current CRAM support and samtools >=1.21. {}".format(exc)
        ) from exc
    return uid_options + options(True)
