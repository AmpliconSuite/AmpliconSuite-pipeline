### Standalone AmpliconSuite-pipeline Conda recipe

This version of the conda recipe is for a standalone build and install process, primarily to debug the building when one of the modules is updated. If you want to build this without using
bioconda to fetch the actual bioconda recipe, you can do the following 

```bash
git clone https://github.com/AmpliconSuite/AmpliconSuite-pipeline
conda build -c bioconda -c conda-forge AmpliconSuite-pipeline/conda-recipe/  # add '--python=3.8 --numpy=1.22.4` if needed
conda create -n "ampsuite_local" python=3.10
conda activate ampsuite_local
conda install -n ampsuite_local -c local ampliconsuite
conda install -n ampsuite_local -c mosek mosek
```
When you're finished clear it with

```bash
conda uninstall -n ampsuite_local -c local -c mosek ampliconsuite mosek
```

### Bioconda recipe updates

The helper script in this directory updates the AmpliconSuite Bioconda recipe
with new component versions, downloads the corresponding GitHub release
tarballs, recalculates their `sha256` values, and resets the build number to
`0` by default.

You do not need a full clone of the Bioconda repository for the usual update
PR. Use a sparse, blobless clone of your fork and check out only this recipe:

```bash
git clone --depth=1 --filter=blob:none --sparse https://github.com/AmpliconSuite/bioconda-recipes.git

cd bioconda-recipes
git sparse-checkout set recipes/ampliconsuite
git switch -c update-ampliconsuite-1.5.3
```

From this repository, dry-run the update first:

```bash
python conda-recipe/update_bioconda_recipe.py \
  --meta-yaml /path/to/bioconda-recipes/recipes/ampliconsuite/meta.yaml \
  --as-version 1.5.3 \
  --aa-version 1.5.r7 \
  --ac-version 1.5.2 \
  --dry-run
```

If the diff looks right, rerun the same command without `--dry-run` to update
the Bioconda checkout:

```bash
python conda-recipe/update_bioconda_recipe.py \
  --meta-yaml /path/to/bioconda-recipes/recipes/ampliconsuite/meta.yaml \
  --as-version 1.5.3 \
  --aa-version 1.5.r7 \
  --ac-version 1.5.2
```

Then commit and push from the Bioconda checkout:

```bash
cd /path/to/bioconda-recipes
git diff recipes/ampliconsuite/meta.yaml
git add recipes/ampliconsuite/meta.yaml
git commit -m "update ampliconsuite to 1.5.3"
git push -u origin update-ampliconsuite-1.5.3
```

Open a PR from your fork to `bioconda/bioconda-recipes`.

The update helper only needs the Python standard library. To run local Bioconda
validation before opening the PR, create a separate environment with
`bioconda-utils`:

```bash
conda create -n bioconda-build -c conda-forge -c bioconda \
  python=3.11 conda-build bioconda-utils
conda activate bioconda-build

cd /path/to/bioconda-recipes
bioconda-utils lint recipes config.yml --packages ampliconsuite
bioconda-utils build recipes config.yml --packages ampliconsuite
```

Full local validation may require more of the Bioconda repository than the
sparse checkout. For a normal recipe update PR, it is reasonable to push the
one-file change and let Bioconda CI run the full checks.
