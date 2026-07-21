# Release Process

Use this checklist for coordinated AmpliconSuite-pipeline (PAA), AmpliconArchitect (AA), AmpliconClassifier (AC), and BFBArchitect releases. Commands that publish tags, releases, images, or pull requests require explicit developer approval. Keep validation outputs outside every repository.

## 1. Define the release set

Record the intended versions and commit SHAs for all four components. Update version files and dependency constraints before testing:

- PAA: `paalib/_version.py`
- AA: `src/AmpliconArchitect.py`
- AC: `ampclasslib/_version.py`
- BFBArchitect: `bfbarchitect/_version.py`
- PAA containers: `docker/requirements/requirements.txt` and `singularity/requirements/requirements.txt`
- AC packaging: `pyproject.toml` and `requirements.txt`
- standalone installer: `install.sh`

An AmpliconSuite-pipeline release is an exact, tested four-tool bundle: one PAA
version pins one AA version, one AC version, and one BFBArchitect version. Keep
that mapping identical in the containers, standalone installer, and Bioconda
recipe. Libraries such as NumPy, Matplotlib, and PuLP may use compatibility
ranges, but component versions must not float independently within a published
AmpliconSuite release.

Confirm every worktree is clean, on its release branch, and synchronized with its upstream:

```bash
git status --short --branch
git rev-list --left-right --count HEAD...@{upstream}
git log -1 --oneline --decorate
```

Review user-facing changes, compatibility notes, licenses, and installation instructions in each repository. Do not update the checked-in Bioconda recipe to unreleased versions: its hashes must refer to immutable GitHub release archives.

## 2. Validate the release candidates

Push the reviewed release commits, but do not create release tags yet. Test those exact remote commits.

For PAA, run the complete unit suite, compilation, shell syntax, and CLI checks:

```bash
python -m unittest discover -s tests -v
python -m py_compile AmpliconSuite-pipeline.py GroupedAnalysisAmpSuite.py paalib/*.py
bash -n install.sh docker/internal_docker_script.sh singularity/*.sh
python AmpliconSuite-pipeline.py --help
python GroupedAnalysisAmpSuite.py --help
```

Also run each component's native test suite. Verify imported package versions and locations with `importlib.metadata`; editable checkouts can leave stale distribution metadata.

Exercise three independent installation surfaces:

1. Run `install.sh` in a disposable, minimal environment and complete a small end-to-end AA+AC analysis.
2. Build Docker from `docker/Dockerfile` with a local candidate tag, then launch a standalone run through `docker/run_paa_docker.py`.
3. Build Singularity from `singularity/ampliconsuite-pipeline.def`, then launch standalone and grouped runs through the Python wrappers.

Cover Clarabel without a Mosek license, a licensed Mosek run, and
AC/BFBArchitect with Gurobi, Mosek, and CBC available. The Gurobi and Mosek
bindings are required parts of a supported installation even though
unrestricted commercial licenses are optional. Also test the no-license
fallback path and confirm the logs identify the solver actually used. Include
at least two distinct samples in grouped testing and verify that AC runs once
across the cohort, writes similarity scores, and respects the shared thread
budget.

Record commands, component versions, image IDs or SIF checksum, input sample/reference, solver, exit status, and output location. Use an external validation directory such as `/tmp/ampliconsuite-release-<date>`; never add BAMs, FASTQs, licenses, or generated results to Git.

## 3. Publish GitHub releases

After candidate tests pass and all repositories are clean, create annotated version tags and GitHub releases in dependency order: BFBArchitect, AA, AC, then PAA. Confirm that each tag resolves to the tested SHA and that the release archive reports the expected CLI version. Do not move or recreate published tags.

## 4. Build and publish containers

Build the final Docker and Singularity artifacts only after the GitHub releases exist. The current recipes download the PAA `master`, AA `master`, and AC `main` branch archives, so build immediately from verified release branch heads and do not merge later work until both builds finish. Re-run the launcher smoke matrix against the final artifacts.

Tag Docker with the PAA version and update `latest` only after validation. Record the Docker digest and SIF SHA-256 before pushing. Registry login and pushes are explicit, separately authorized steps.

## 5. Submit the Bioconda update

Prepare the sparse checkout and dry-run before release day, but calculate final
hashes only after the GitHub tags exist. Use
`conda-recipe/update_bioconda_recipe.py` as described in
`conda-recipe/README.md`, inspect the `meta.yaml`, `build.sh`, and BFBA solver
patch, run available lint/build checks, then open the manual PR. Verify that the
single `ampliconsuite` recipe embeds the exact released versions of all four
components and resets its build number appropriately. Do not create a separate
Bioconda BFBArchitect package for this release bundle. Treat any AC or BFBA
installation failure as a package-build failure; never substitute an older AC
release or omit one of the four components.

Audit new transitive dependencies before updating versions. AC 2.0 makes
BFBArchitect a core dependency, so the recipe installs the pinned BFBA source
with `--no-deps` and represents its open-source dependencies explicitly. The
Conda compatibility patch allows PuLP 2.8 or newer and selects external CBC
when PuLP has no bundled executable; it must continue to retain BFBA's required
`gurobipy` metadata. Because Bioconda cannot resolve the Gurobi and Mosek vendor
channels, install those bindings in the clean consumer environment before
running `install.sh --finalize_only` and the final solver tests.

Finish by recording release URLs, container digests, the Bioconda PR, and any intentionally deferred issues in the release notes.
