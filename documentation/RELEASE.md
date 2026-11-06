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
3. Convert the candidate Docker image to a SIF with `apptainer build <candidate>.sif docker-daemon://<candidate tag>`, then launch standalone and grouped runs through the Python wrappers. The wrappers use `apptainer` if it is available and `singularity` otherwise.

`singularity/ampliconsuite-pipeline.def` remains available for building a SIF without Docker, but published SIFs come from the Docker image so that both distributions have the same contents. If the definition changed, build it and check that it runs.

The Docker and Apptainer recipes download PAA and AA from `master` and AC from `main`. From the local checkout they take only the requirements file and the internal scripts, and the launchers run on the host. An image built from a feature branch therefore has that branch's pins but `master`'s code. Merge candidate code before testing it in an image (no tag is needed), and use `--no-cache --pull` so that cached layers can't hide a change. Confirm what was tested by hashing a changed file inside the image and comparing it with `git show <sha>:<path> | sha256sum`.

Published images must contain samtools >=1.21. Check the installed version inside
each final Docker/SIF artifact; both build definitions enforce this minimum.
Include a CRAM 3.1 launcher smoke test that populates a persistent reference cache,
then repeats with the populated cache offline. Verify standalone and grouped
launchers reuse the cache and stop on reference-preparation failures.

Cover Clarabel without a Mosek license, a licensed Mosek run, and
AC/BFBArchitect with Gurobi, Mosek, and CBC available. The Gurobi and Mosek
bindings are checked by the installer; a missing Gurobi binding fails the
installation, while a Mosek binding that cannot be installed or loaded produces
a warning and AA can fall back to Clarabel. Unrestricted commercial licenses
are optional. Also test the no-license fallback path and confirm the logs
identify the solver actually used. Gurobi, Mosek, and CBC can legitimately give
different BFB scores, so compare each solver only against its own earlier
results. Include
at least two distinct samples in grouped testing and verify that AC runs once
across the cohort, writes similarity scores, and respects the shared thread
budget.

Record commands, component versions, image IDs or SIF checksum, input sample/reference, solver, exit status, and output location. Keep validation outputs in an external directory on real disk, and never add BAMs, FASTQs, licenses, or generated results to Git. `/tmp` is RAM-backed on many systems, so large test files there compete with the running jobs for memory. Run heavy jobs such as AA matrices and Conda installs one at a time, and check available memory before starting each one.

## 3. Publish GitHub releases

After candidate tests pass and all repositories are clean, create annotated version tags and GitHub releases in dependency order: BFBArchitect, AA, AC, then PAA. Confirm that each tag resolves to the tested SHA and that the release archive reports the expected CLI version. Do not move or recreate published tags.

## 4. Build and publish containers

Build the final Docker image only after the GitHub releases exist, then convert it to the final SIF. The Docker recipe downloads the PAA `master`, AA `master`, and AC `main` branch archives, so build immediately from verified release branch heads and do not merge later work until the build finishes. Re-run the launcher smoke matrix against the final artifacts.

Tag Docker with the PAA version and update `latest` only after validation. Record the Docker digest and SIF SHA-256 before pushing. Registry login and pushes are explicit, separately authorized steps.

Docker images go to `jluebeck/ampliconsuite-pipeline` on Docker Hub. SIFs go to the GitHub Container Registry as ORAS artifacts, under the same version and `latest` tags:

```bash
apptainer build ampliconsuite-pipeline_<version>.sif docker-daemon://jluebeck/ampliconsuite-pipeline:<version>
apptainer registry login --username <github-user> oras://ghcr.io  # token with write:packages
apptainer push ampliconsuite-pipeline_<version>.sif oras://ghcr.io/ampliconsuite/ampliconsuite-pipeline:<version>
apptainer push ampliconsuite-pipeline_<version>.sif oras://ghcr.io/ampliconsuite/ampliconsuite-pipeline:latest
```

New GHCR packages start private, so set the package to public after its first push. Then confirm that the commands in the README work without logging in, for `docker pull` and for both the `oras://` and `docker://` Apptainer routes. Use an empty `DOCKER_CONFIG` and `HOME`, and remove the local Docker image first so that the layers are really downloaded.

The GHCR package holds a SIF, not a Docker image, so `docker pull` on the GHCR address fails with `unsupported media type`. Docker users pull from Docker Hub.

`apptainer push` doesn't set the description on the GHCR package page. GHCR reads it from the manifest's `org.opencontainers.image.description` annotation, which is limited to 512 characters. Add it after pushing, using a GitHub token with `write:packages`:

```bash
REPO=ampliconsuite/ampliconsuite-pipeline
TYPE=application/vnd.oci.image.manifest.v1+json
DESC='AmpliconSuite-pipeline Apptainer/Singularity image (SIF). Pull with: apptainer pull ampliconsuite-pipeline.sif oras://ghcr.io/ampliconsuite/ampliconsuite-pipeline:latest . For Docker, use the Docker Hub image: docker pull jluebeck/ampliconsuite-pipeline . Running docker pull on this GHCR address will not work.'
TOKEN=$(curl -s -u "<github-user>:$(gh auth token)" "https://ghcr.io/token?service=ghcr.io&scope=repository:$REPO:pull,push" | jq -r .token)
for tag in <version> latest; do
  curl -s -H "Authorization: Bearer $TOKEN" -H "Accept: $TYPE" "https://ghcr.io/v2/$REPO/manifests/$tag" \
    | jq -cj --arg d "$DESC" '.annotations["org.opencontainers.image.description"] = $d' > manifest.json
  curl -s -o /dev/null -w "$tag HTTP %{http_code}\n" -X PUT -H "Authorization: Bearer $TOKEN" -H "Content-Type: $TYPE" \
    --data-binary @manifest.json "https://ghcr.io/v2/$REPO/manifests/$tag"
done
```

This leaves the SIF layer unchanged but gives the tags a new manifest digest, so record the digest after this step. The earlier manifest stays in the package, untagged. Pull the tag again and confirm that the SIF's SHA-256 still matches.

## 5. Submit the Bioconda update

Prepare the sparse checkout and dry-run before release day, but calculate final
hashes only after the GitHub tags exist. Use
`conda-recipe/update_bioconda_recipe.py` as described in
`conda-recipe/README.md`, inspect the `meta.yaml` and `build.sh`, run available
lint/build checks, then open the manual PR. Verify that the
single `ampliconsuite` recipe embeds the exact released versions of all four
components and resets its build number appropriately. Do not create a separate
Bioconda BFBArchitect package for this release bundle. Treat any AC or BFBA
installation failure as a package-build failure; never substitute an older AC
release or omit one of the four components.

Audit new transitive dependencies before updating versions. AC 2.0 makes
BFBArchitect a core dependency, so the recipe installs the pinned BFBA source
with `--no-deps` and represents its open-source dependencies explicitly. BFBA
selects external CBC when PuLP has no bundled executable. Because Bioconda
cannot obtain Gurobi or Mosek, install Gurobi from its vendor channel in the
clean consumer environment, then
run `install.sh --finalize_only`, which uses pip to add Mosek, before the final
solver tests.

Test the package Bioconda built, not only a local build. Once the PR's CI passes, download the artifact with `gh run download <run-id> -R bioconda/bioconda-recipes -n linux-packages`, then install it in a fresh environment with `-c file://<packages dir> -c conda-forge -c bioconda`. Add the vendor solvers as above and repeat the real-data solver and grouped tests. First unset `AA_SRC`, `AC_SRC`, `REF_CACHE`, `REF_PATH`, `PYTHONPATH`, and any license variables, and use an empty `HOME`, or the run can silently use development checkouts. Afterwards, search the outputs for those checkout paths. The CI also builds a BioContainers image, which is a useful no-license check, but it has no vendor solvers.

To request review, comment `@BiocondaBot please add label` on the PR. Update the checked-in `conda-recipe/meta.yaml` once the GitHub tags exist.

Finish by recording release URLs, container digests, the Bioconda PR, and any intentionally deferred issues in the release notes.
