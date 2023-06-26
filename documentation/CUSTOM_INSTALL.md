## Performing a standalone custom install from each individual module.

The easiest way to do a standalone install is to follow the instructions from the README and use the `install.sh` script.
If that or none of the other installation options are possible (e.g. you need to use `python2`, you can install each sub-module and dependency individually.

1. Clone the AmpliconSuite-pipeline git rep:

`git clone https://github.com/jluebeck/AmpliconSuite-pipeline.git`

2. Individually install other prerequisites from the section below following the install instrucitons on each.
3. Set the location where you would like to store the `$AA_DATA_REPO`
```bash
        mkdir data_repo && cd data_repo
        # copy or download files into data_repo directory
        wget [url for data repo [hg19/GRCh37/GRCh38/mm10].tar.gz]
        tar -xzf [hg19/GRCh37/GRCh38/mm10].tar.gz

        echo export AA_DATA_REPO=$PWD >> ~/.bashrc
        touch coverage.stats && chmod a+r coverage.stats
        source ~/.bashrc
```
4. Run `install.sh --finalize_only` script from AmpliconSuite-pipeline.

## Prerequisites for standalone installation:
AmpliconSuite-pipeline supports both `python2` and `python3`, however CNVkit requires `python3`. `Python3` support for AmpliconArchitect was added in version 1.3. 

Unless you are using a containerized version, and depending on what input data you are starting from, AmpliconSuite-pipeline may require the following tools to be installed beforehand:
- (required) The [AmpliconSuite/AmpliconArchictect fork](https://github.com/AmpliconSuite/AmpliconArchitect) must be installed. Instructions for that are available [here](https://github.com/AmpliconSuite/AmpliconArchitect/blob/master/docs/standalone_usage.md).
- (required) The latest AmpliconArchitect [data repo](https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/).
  - versions of the data repos containing bwa index files are also provided [here](https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/). Indexed version recommended if starting from unaligned fastq reads. Instructions for setting up the AA data repo are available [here](https://github.com/AmpliconSuite/AmpliconArchitect/blob/master/docs/standalone_usage.md).
- (recommended) [AmpliconClassifier](https://github.com/jluebeck/AmpliconClassifier) to generate classifications of AmpliconArchitect outputs.
- (recommended) [CNVkit](https://github.com/etal/cnvkit) to generate CNV calls for focal amplification seed region identification.
- (optional) [bwa mem](https://github.com/lh3/bwa) (unless supplying your own BAM file)
- (optional) [samtools](http://www.htslib.org/) (unless you already have a coordinate-sorted and indexed BAM file).
- Scripts packaged with AmpliconSuite-pipeline require the `numpy`, `matplotlib` and `intervaltree` python packages. Those packages can be installed with `pip`, `conda` or similar.

AmpliconSuite-pipeline assumes both `samtools` and `bwa` executables are on the system path and can be directly invoked from bash without pathing to the executables. AmpliconSuite-pipeline will generate a BWA index for the reference genome if one is not yet in place. This adds >1hr to running time for the first use only when alignment is performed. Data repos with BWA index pre-generated are available [here](https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/). AmpliconSuite-pipeline will also function on coordinate-sorted CRAM files, [provided that the CRAM reference is in place](http://www.htslib.org/workflow/#:~:text=One%20of%20the%20key%20concepts,genome%20used%20to%20generate%20it.).

AmpliconSuite-pipeline has been tested with Ubuntu (16.04 and above) and CentOS 7. AmpliconSuite-pipeline's optional dependencies related to CNV calling will not work on CentOS 6.

**Note on using CNVkit**: We currently recommend using CNVkit for identification of AA seeds. CNVkit requires
`python3`. It also requires `R` version >= 3.5, which is non-standard on Ubuntu 16.04/14.04.

## Getting `mscorefonts` onto your system.
AmpliconArchitect figures will attempt to use the Arial font, and will fall back to the default `matplotlib` "Deja Vu Sans" font. On macOS, Arial will likely already be present. 
Install the `mscorefonts` package one of two ways:

a) First run `conda install mscorefonts` then launch python and do
```python
import matplotlib.font_manager
matplotlib.font_manager._load_fontmanager(try_read_cache=False)
```

b) (Ubuntu) `sudo apt update && sudo apt install ttf-mscorefonts-installer`. Then do `sudo fc-cache -f -v` to rebuild the font cache.