# AmpliconSuite Bioconda recipe

The `ampliconsuite` Conda package is one coordinated four-tool release. Its
recipe downloads and installs exact versions of AmpliconSuite-pipeline,
AmpliconArchitect, AmpliconClassifier, and BFBArchitect. These four versions
must be updated and tested together; general runtime libraries may move within
their declared compatibility constraints.

BFBArchitect is embedded in this recipe rather than published as a separate
Bioconda package. It uses Conda's external CBC executable when PuLP has no
bundled CBC. Bioconda cannot obtain Gurobi or Mosek during its build and neither
is on conda-forge, so they are installed afterward: users add Gurobi from its
vendor channel, and `install.sh --finalize_only` uses pip to add Mosek. The
unrestricted license files remain optional.

## Build and test locally

```bash
conda build -c conda-forge -c bioconda AmpliconSuite-pipeline/conda-recipe/
conda create -n ampsuite_local python=3.10
conda activate ampsuite_local
conda install -n ampsuite_local -c local -c conda-forge -c bioconda ampliconsuite
conda install -n ampsuite_local -c gurobi gurobi
./install.sh --finalize_only
```

Installing Gurobi from its vendor channel and Mosek through pip via
`install.sh --finalize_only` makes both commercial solver paths available;
providing a license determines whether unrestricted models can use them.
Clarabel and CBC are installed by the AmpliconSuite recipe and require no
license. Validate a real BFB-positive example in this complete consumer
environment before release.

## Update the Bioconda recipe

Release archives must exist before calculating hashes. A sparse checkout is
enough for the normal update:

```bash
git clone --depth=1 --filter=blob:none --sparse \
  https://github.com/AmpliconSuite/bioconda-recipes.git
cd bioconda-recipes
git sparse-checkout set recipes/ampliconsuite
git switch -c update-ampliconsuite-1.6.2
```

The published recipe already contains the four-tool structure. From the
AmpliconSuite-pipeline checkout, dry-run the version and hash update:

```bash
python conda-recipe/update_bioconda_recipe.py \
  --meta-yaml conda-recipe/meta.yaml \
  --as-version 1.6.2 \
  --aa-version 1.6.r2 \
  --ac-version 2.0.1 \
  --bfba-version 1.0.2 \
  --dry-run
```

Rerun without `--dry-run`, then synchronize all reviewed recipe files as
needed. Keep `build.sh` identical to the version maintained in this repository:

```bash
cp conda-recipe/meta.yaml \
  /path/to/bioconda-recipes/recipes/ampliconsuite/meta.yaml
cp conda-recipe/build.sh \
  /path/to/bioconda-recipes/recipes/ampliconsuite/build.sh
```

The update helper does not change anything else in the recipe. The published
recipe for BFBArchitect 1.0.1 has a solver patch and a PuLP limit, so on the
first update to BFBArchitect 1.0.2, delete `bfba-conda-solvers.patch` and its
`patches:` entry, and make the requirements match `conda-recipe/meta.yaml`.

Inspect both files and reset the recipe build number to `0` for a new
AmpliconSuite version. The update helper uses only the Python standard library;
it updates the four exact versions, downloads all four release archives,
recalculates their SHA-256 hashes, and updates the build number.

The build intentionally fails if either the pinned AC or BFBA source cannot be
installed. It must never substitute an older AC layout or silently omit BFBA;
the four requested component versions are an atomic release unit.

For local Bioconda validation, use a separate `bioconda-utils` environment:

```bash
bioconda-utils lint recipes config.yml --packages ampliconsuite
bioconda-utils build recipes config.yml --packages ampliconsuite
```

After CI succeeds, install the produced package into a fresh environment, add
Gurobi from its vendor channel, run `install.sh --finalize_only` to install
Mosek through pip, and repeat the version, CBC, Gurobi, and real-data checks
before requesting merge.
