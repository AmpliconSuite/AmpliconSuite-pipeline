#!/usr/bin/env python3
"""Update the AmpliconSuite Bioconda recipe versions and source hashes."""

import argparse
import difflib
import hashlib
import re
import sys
import urllib.request
from pathlib import Path


COMPONENTS = (
    ("AS", "AS_version"),
    ("AA", "AA_version"),
    ("AC", "AC_version"),
    ("BFBA", "BFBA_version"),
)

SOURCES = (
    (
        "AS",
        "AS_version",
        "https://github.com/AmpliconSuite/AmpliconSuite-pipeline/archive/v{version}.tar.gz",
    ),
    (
        "AA",
        "AA_version",
        "https://github.com/AmpliconSuite/AmpliconArchitect/archive/v{version}.tar.gz",
    ),
    (
        "AC",
        "AC_version",
        "https://github.com/AmpliconSuite/AmpliconClassifier/archive/v{version}.tar.gz",
    ),
    (
        "BFBA",
        "BFBA_version",
        "https://github.com/AmpliconSuite/BFBArchitect/archive/v{version}.tar.gz",
    ),
)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Update recipes/ampliconsuite/meta.yaml with new AmpliconSuite, "
            "AmpliconArchitect, AmpliconClassifier, and BFBArchitect versions "
            "and sha256 sums."
        )
    )
    parser.add_argument(
        "--meta-yaml",
        required=True,
        type=Path,
        help="Path to the Bioconda recipes/ampliconsuite/meta.yaml file.",
    )
    parser.add_argument("--as-version", required=True, help="AmpliconSuite-pipeline version.")
    parser.add_argument("--aa-version", required=True, help="AmpliconArchitect version.")
    parser.add_argument("--ac-version", required=True, help="AmpliconClassifier version.")
    parser.add_argument("--bfba-version", required=True, help="BFBArchitect version.")
    parser.add_argument(
        "--build-number",
        type=int,
        default=0,
        help="Build number to set in the recipe. Defaults to 0 for version bumps.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the recipe diff without writing changes.",
    )
    return parser.parse_args()


def source_sha256(url):
    digest = hashlib.sha256()
    with urllib.request.urlopen(url) as response:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)

    return digest.hexdigest()


def update_jinja_version(text, variable, version):
    pattern = re.compile(r'(\{%\s*set\s+' + re.escape(variable) + r'\s*=\s*")[^"]+("\s*%\})')
    text, count = pattern.subn(r"\g<1>" + version + r"\2", text, count=1)
    if count != 1:
        raise ValueError(f"Could not find exactly one Jinja assignment for {variable}.")

    return text


def update_sha_after_source_url(text, variable, sha256):
    source_pattern = re.compile(
        r"(?P<prefix>- url:\s*https://github\.com/AmpliconSuite/[^/\n]+/archive/v"
        + r"\{\{\s*"
        + re.escape(variable)
        + r"\s*\}\}\.tar\.gz\s*\n\s*sha256:\s*)"
        + r"(?P<sha>[0-9a-fA-F]{64})",
        re.MULTILINE,
    )
    text, count = source_pattern.subn(r"\g<prefix>" + sha256, text, count=1)
    if count != 1:
        raise ValueError(f"Could not find exactly one source sha256 for {variable}.")

    return text


def update_build_number(text, build_number):
    pattern = re.compile(r"(^build:\s*\n(?:^[ \t].*\n)*?^[ \t]*number:\s*)\d+", re.MULTILINE)
    text, count = pattern.subn(r"\g<1>" + str(build_number), text, count=1)
    if count != 1:
        raise ValueError("Could not find exactly one build:number entry.")

    return text


def update_recipe(original_text, versions, hashes, build_number):
    updated_text = original_text
    for short_name, variable in COMPONENTS:
        updated_text = update_jinja_version(updated_text, variable, versions[short_name])

    for short_name, variable, _url_template in SOURCES:
        updated_text = update_sha_after_source_url(updated_text, variable, hashes[short_name])

    return update_build_number(updated_text, build_number)


def print_diff(path, before, after):
    diff = difflib.unified_diff(
        before.splitlines(keepends=True),
        after.splitlines(keepends=True),
        fromfile=str(path),
        tofile=str(path),
    )
    sys.stdout.writelines(diff)


def main():
    args = parse_args()
    versions = {
        "AS": args.as_version,
        "AA": args.aa_version,
        "AC": args.ac_version,
        "BFBA": args.bfba_version,
    }

    if not args.meta_yaml.exists():
        raise SystemExit(f"Recipe file does not exist: {args.meta_yaml}")

    hashes = {}
    for short_name, _variable, url_template in SOURCES:
        url = url_template.format(version=versions[short_name])
        print(f"Fetching {short_name} source: {url}", file=sys.stderr)
        hashes[short_name] = source_sha256(url)
        print(f"{short_name} sha256: {hashes[short_name]}", file=sys.stderr)

    original_text = args.meta_yaml.read_text()
    updated_text = update_recipe(original_text, versions, hashes, args.build_number)

    if original_text == updated_text:
        print("No recipe changes needed.", file=sys.stderr)
        return

    if args.dry_run:
        print_diff(args.meta_yaml, original_text, updated_text)
    else:
        args.meta_yaml.write_text(updated_text)
        print(f"Updated {args.meta_yaml}", file=sys.stderr)


if __name__ == "__main__":
    main()
