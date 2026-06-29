"""
Data repository download functionality for AmpliconSuite-pipeline
"""
import logging
import os
import sys
import tarfile
import time


DATA_REPO_BASE_URL = "https://refs.ampliconrepository.org/data/module_support_files/AmpliconArchitect/"

# Freshness-check tuning
REPO_CHECK_TTL = 86400        # seconds to trust an "up to date" verdict before re-checking (24h)
REPO_CHECK_NETFAIL_TTL = 3600  # seconds to wait before re-attempting after a network failure (1h)
REPO_CHECK_TIMEOUT = 5         # seconds to wait on the network before giving up
REPO_CHECK_CACHE_NAME = ".repo_check_cache"  # stored at the top of $AA_DATA_REPO


def download_file(url, destination_folder):
    import urllib.request  # here because python2 not work with it
    filename = os.path.join(destination_folder, url.split("/")[-1])
    try:
        response = urllib.request.urlopen(url)
        file_size = int(response.headers.get('Content-Length', 0))
        response.close()
        file_size_gb = round(file_size / (1024**3), 2)
        if file_size_gb > 0.1:
            print("\nDownloading " + url + " ... (" + str(file_size_gb) + "GB)")
        else:
            print("\nDownloading " + url + " ...")

        bar_width = 40

        def reporthook(count, block_size, total_size):
            if total_size <= 0:
                downloaded = count * block_size
                sys.stdout.write("\r  Downloaded: {:.1f} MB".format(downloaded / (1024**2)))
            else:
                downloaded = min(count * block_size, total_size)
                frac = downloaded / total_size
                filled = int(bar_width * frac)
                bar = '#' * filled + '-' * (bar_width - filled)
                pct = frac * 100
                mb_done = downloaded / (1024**2)
                mb_total = total_size / (1024**2)
                sys.stdout.write("\r  [{bar}] {pct:5.1f}%  {done:.1f}/{total:.1f} MB".format(
                    bar=bar, pct=pct, done=mb_done, total=mb_total))
            sys.stdout.flush()

        urllib.request.urlretrieve(url, filename, reporthook=reporthook)
        sys.stdout.write("\n")
        print("File downloaded and saved to: " + str(filename))
    except Exception as e:
        sys.stderr.write("Failed to download file. Error: " + str(e) + "\n")


def extract_tar_gz(file_path, destination_folder):
    if not file_path.endswith('.tar.gz'):
        sys.stderr.write("Cannot extract file " + file_path)
        sys.exit(1)

    with tarfile.open(file_path, 'r:gz') as tar:
        tar.extractall(destination_folder)

    os.remove(file_path)


def handle_repo_download(args, AA_REPO):
    """Handle data repository download and exit"""
    for ref in args.download_repo:
        print("Downloading " + ref + " into " + AA_REPO)
        ref_base_url = DATA_REPO_BASE_URL + ref
        md5file = ref_base_url + "_md5sum.txt"
        ref_file = ref_base_url + ".tar.gz"

        # The extracted directory is the bare reference name; the '_indexed' tarball variant extracts to the same dir.
        ref_dir = ref[:-len("_indexed")] if ref.endswith("_indexed") else ref
        if os.path.exists(AA_REPO + ref_dir):
            # Acknowledge whether the version being pulled differs from what is already installed.
            local_md5 = _read_md5_token(os.path.join(AA_REPO, ref + "_md5sum.txt"))
            remote_md5 = _fetch_remote_md5(ref)
            if local_md5 and remote_md5 and local_md5 == remote_md5:
                print("An AA data repo for " + ref_dir + " already exists and is already up to date (md5 matches "
                      "the online version). It will be re-downloaded and replaced.")
            elif local_md5 and remote_md5:
                print("An AA data repo for " + ref_dir + " already exists but is OUT OF DATE compared to the online "
                      "version (md5 differs). It will be replaced with the newer version.")
            else:
                print("An AA data repo directory already exists for " + ref_dir + " and it will be replaced! "
                      "(Could not compare versions by md5.)")

        download_file(md5file, AA_REPO)
        download_file(ref_file, AA_REPO)
        print("Extracting...\n")
        extract_tar_gz(AA_REPO + ref + ".tar.gz", AA_REPO)
    
    print("Finished")
    sys.exit(0)


def _read_md5_token(path):
    """Read the md5 hash (first whitespace-delimited token) from an md5sum.txt-style file/URL content."""
    try:
        with open(path) as infile:
            return infile.read().split()[0].strip()
    except (IOError, OSError, IndexError):
        return None


def _local_md5_variant(AA_REPO, ref):
    """
    Determine which data-repo variant the user has locally, based on which '<ref>_md5sum.txt' file is present
    at the top of $AA_DATA_REPO. Returns (variant, local_md5) where variant is e.g. 'GRCh38' or 'GRCh38_indexed',
    or (None, None) if no md5sum file is present (a very old or manually-installed repo).
    """
    for variant in (ref + "_indexed", ref):
        local_md5_path = os.path.join(AA_REPO, variant + "_md5sum.txt")
        if os.path.isfile(local_md5_path):
            return variant, _read_md5_token(local_md5_path)

    return None, None


def _fetch_remote_md5(variant):
    """Fetch the remote md5 hash for a data-repo variant. Returns the hash string, or None on any failure."""
    import urllib.request  # local import: python2 compatibility, matches download_file()
    url = DATA_REPO_BASE_URL + variant + "_md5sum.txt"
    try:
        with urllib.request.urlopen(url, timeout=REPO_CHECK_TIMEOUT) as response:
            return response.read().decode("utf-8").split()[0].strip()
    except Exception:
        return None


def _read_check_cache(AA_REPO):
    """Read the freshness-check cache. Returns {key: (status, value, timestamp)}."""
    cache = {}
    cache_path = os.path.join(AA_REPO, REPO_CHECK_CACHE_NAME)
    try:
        with open(cache_path) as infile:
            for line in infile:
                fields = line.split()
                if len(fields) == 4:
                    key, status, value, ts = fields
                    try:
                        cache[key] = (status, value, float(ts))
                    except ValueError:
                        continue
    except (IOError, OSError):
        pass

    return cache


def _write_check_cache(AA_REPO, cache):
    """Persist the freshness-check cache. Best-effort: a read-only repo simply means we re-check next time."""
    cache_path = os.path.join(AA_REPO, REPO_CHECK_CACHE_NAME)
    try:
        with open(cache_path, 'w') as outfile:
            for key, (status, value, ts) in cache.items():
                outfile.write("{} {} {} {}\n".format(key, status, value, ts))
    except (IOError, OSError):
        pass


def _local_last_updated(AA_REPO, ref):
    """Return the human-readable construction date of the local repo for messaging, or None."""
    try:
        with open(os.path.join(AA_REPO, ref, "last_updated.txt")) as infile:
            return infile.read().strip()
    except (IOError, OSError):
        return None


def check_repo_freshness(AA_REPO, ref, strict=False):
    """
    Compare the locally-installed $AA_DATA_REPO for `ref` against the latest version hosted online and warn the
    user (or exit, if strict) when it is out of date.

    The comparison uses the small '<ref>_md5sum.txt' files: the one stored locally at download time vs. a freshly
    fetched remote copy. Results are cached in $AA_DATA_REPO/.repo_check_cache so that an up-to-date repo only
    incurs one network call per REPO_CHECK_TTL window. An out-of-date repo is NOT cached, so the warning fires on
    every run until the user updates. The check always fails open: any network error is a soft note, never a block.
    """
    last_updated = _local_last_updated(AA_REPO, ref)
    if last_updated:
        logging.info(ref + " data repo constructed on " + last_updated)

    variant, local_md5 = _local_md5_variant(AA_REPO, ref)

    # Case: neither an md5sum nor a last_updated.txt -> a repo predating version tracking entirely.
    if variant is None and last_updated is None:
        msg = ("Your $AA_DATA_REPO for " + ref + " predates data-repo version tracking and is very likely out of "
               "date. Please re-download it with '--download_repo " + ref + "' (or '" + ref + "_indexed').")
        _emit_staleness(msg, strict)
        return

    # Case: we have a last_updated.txt but no md5sum.txt -> can't do the comparison; soft, non-blocking note.
    if variant is None:
        logging.warning("Could not verify whether your $AA_DATA_REPO for " + ref + " is up to date (no "
                        "'" + ref + "_md5sum.txt' found). If you did not download it via '--download_repo', consider "
                        "re-downloading to enable update checks. Use --no_repo_check to silence this.\n")
        return

    key = variant
    now = time.time()
    cache = _read_check_cache(AA_REPO)
    cached = cache.get(key)
    if cached:
        status, value, ts = cached
        # Trust a recent "up to date" verdict only if the local content is unchanged since we confirmed it.
        if status == "ok" and value == local_md5 and (now - ts) < REPO_CHECK_TTL:
            return
        # Avoid hammering an unreachable server after a recent network failure.
        if status == "netfail" and (now - ts) < REPO_CHECK_NETFAIL_TTL:
            return

    remote_md5 = _fetch_remote_md5(variant)
    if remote_md5 is None:
        logging.info("Could not reach the data-repo server to check for updates; continuing. "
                     "Use --no_repo_check to disable this check.\n")
        cache[key] = ("netfail", "na", now)
        _write_check_cache(AA_REPO, cache)
        return

    if local_md5 == remote_md5:
        cache[key] = ("ok", local_md5, now)
        _write_check_cache(AA_REPO, cache)
        return

    # Out of date: warn (or exit) every run, and do NOT cache an "ok" verdict.
    date_note = (" (last updated " + last_updated + ")") if last_updated else ""
    msg = ("Your $AA_DATA_REPO for " + ref + date_note + " is OUT OF DATE -- a newer version is available online. "
           "Update it by running: AmpliconSuite-pipeline.py --download_repo " + variant)
    _emit_staleness(msg, strict)


def _emit_staleness(msg, strict):
    """Emit a data-repo staleness message as a hard error (strict) or a warning, with the disable hint appended."""
    hint = " (To skip this check, pass --no_repo_check or set AS_NO_REPO_CHECK=1.)"
    if strict:
        logging.error(msg + hint + "\n")
        sys.exit(1)
    else:
        logging.warning(msg + hint + "\n")
