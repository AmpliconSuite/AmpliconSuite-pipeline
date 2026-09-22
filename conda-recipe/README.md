# AmpliconSuite Bioconda recipe

The `ampliconsuite` Conda package is one coordinated four-tool release. Its
recipe downloads and installs exact versions of AmpliconSuite-pipeline,
AmpliconArchitect, AmpliconClassifier, and BFBArchitect. These four versions
must be updated and tested together; general runtime libraries may move within
their declared compatibility constraints.

BFBArchitect is embedded in this recipe rather than published as a separate
Bioconda package. The `bfba-conda-solvers.patch` file permits PuLP 2.8 or newer
and selects Conda's external CBC executable when PuLP has no bundled CBC. It
retains BFBArchitect's mandatory `gurobipy` metadata. Bioconda cannot obtain
Gurobi or Mosek from vendor channels during its build, so users install those
bindings separately; their unrestricted license files remain optional.

## Build and test locally

```bash
conda build -c conda-forge -c bioconda AmpliconSuite-pipeline/conda-recipe/
conda create -n ampsuite_local python=3.10
conda activate ampsuite_local
conda install -n ampsuite_local -c local -c conda-forge -c bioconda ampliconsuite
conda install -n ampsuite_local -c gurobi gurobi
conda install -n ampsuite_local -c mosek mosek
./install.sh --finalize_only
```

Installing the Gurobi and Mosek packages wires both commercial solver paths;
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
git switch -c update-ampliconsuite-1.6.1
```

The published recipe already contains the four-tool structure. From the
AmpliconSuite-pipeline checkout, dry-run the version and hash update:

```bash
python conda-recipe/update_bioconda_recipe.py \
  --meta-yaml conda-recipe/meta.yaml \
  --as-version 1.6.1 \
  --aa-version 1.6.r1 \
  --ac-version 2.0.0 \
  --bfba-version 1.0.1 \
  --dry-run
```

Rerun without `--dry-run`, then synchronize all reviewed recipe files as
needed. Keep `build.sh` and the BFBA patch identical to the versions maintained
in this repository:

```bash
cp conda-recipe/meta.yaml \
  /path/to/bioconda-recipes/recipes/ampliconsuite/meta.yaml
cp conda-recipe/build.sh \
  /path/to/bioconda-recipes/recipes/ampliconsuite/build.sh
cp conda-recipe/bfba-conda-solvers.patch \
  /path/to/bioconda-recipes/recipes/ampliconsuite/bfba-conda-solvers.patch
```

Inspect all three files and reset the recipe build number to `0` for a new
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
the two vendor solver packages, run `install.sh --finalize_only`, and repeat the
version, CBC, Gurobi, and real-data checks before requesting merge.
