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
