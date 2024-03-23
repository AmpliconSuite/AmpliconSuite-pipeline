### Standalone AmpliconSuite-pipeline Conda recipe

This version of the conda recipe is for a standalone build and install process, primarily to debug the building when one of the modules is updated. If you want to build this without using
bioconda to fetch the actual bioconda recipe, you can do the following 

```bash
conda build AmpliconSuite-pipeline/conda-recipe/  # add '--python=3.8 --numpy=1.22.4` if needed for older systems
conda create -n "ampsuite" python>=3.8.0
conda activate ampsuite
conda install -n ampsuite -c local -c mosek ampliconsuite mosek
```
When you're finished clear it with

```bash
conda uninstall -n ampsuitetest -c local -c mosek ampliconsuite mosek
```
