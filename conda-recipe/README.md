### Standalone AmpliconSuite-pipeline Conda recipe

This version of the conda recipe is for a standalone build and install process. If you want to build this without using
bioconda to fetch the actual bioconda recipe, you can do the following 

```bash
conda build AmpliconSuite-pipeline/conda-recipe/ 
conda create -n "ampsuite" python>=3.8.0
conda activate ampsuite
conda install -n ampsuite -c local -c mosek ampliconsuite mosek
```