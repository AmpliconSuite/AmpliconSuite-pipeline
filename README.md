# AmpliconSuite-pipeline
![GitHub](https://img.shields.io/github/license/jluebeck/AmpliconSuite-pipeline)

A multithread-enabled end-to-end wrapper for [AmpliconArchitect](https://github.com/jluebeck/AmpliconArchitect) and the associate tools for data preparation and interpretation. 
Performs preliminary steps (alignment, seed detection, & seed filtering) required prior to running AmpliconArchitect. AmpliconSuite-pipeline can be invoked to begin at any intermediate stage of the data preparation process and can itself invoke both AmpliconArchitect and the downstream tool AmpliconClassifier, which is used to classify ecDNAs and BFB. AmpliconSuite-pipeline was formerly called "PrepareAA".

AmpliconSuite-pipeline supports hg19, GRCh37, GRCh38 (hg38), and mouse genome mm10 (GRCm38). The tool also supports analysis with a human-viral hybrid reference genome we provide, "GRCh38_viral", which can be used to detect oncoviral hybrid focal amplifications and ecDNA in cancers with oncoviral infections such as HPV and HBV.

**Current version: 0.1458.0**

[comment]: # (Versioning based on major_version.days_since_initial_commit.minor_version. Initial commit: March 5th, 2019)

Please check out our [**detailed guide**](https://github.com/jluebeck/PrepareAA/blob/master/GUIDE.md) on running to learn about best practices and see some FAQs.


## Prerequisites:
AmpliconSuite-pipeline supports both `python2` and `python3`, however CNVKit requires `python3`. `Python3` support for AmpliconArchitect was added in version 1.3. 

Unless you are using a containerized version, and depending on what input data you are starting from, AmpliconSuite-pipeline may require the following tools to be installed beforehand:
- (required) The [jluebeck/AmpliconArchictect fork](https://github.com/jluebeck/AmpliconArchitect) must be installed.
- (required) The latest AmpliconArchitect [data repo](https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/).
  - versions of the data repos containing bwa index files are also provided [here](https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/). Indexed version recommended if starting from unaligned fastq reads.
- (recommended) [AmpliconClassifier](https://github.com/jluebeck/AmpliconClassifier) to generate classifications of AmpliconArchitect outputs.
- (recommended) [CNVkit](https://github.com/etal/cnvkit) to generate CNV calls for focal amplification seed region identification.
- (optional) [bwa mem](https://github.com/lh3/bwa) (unless supplying your own BAM file)
- (optional) [samtools](http://www.htslib.org/) (unless you already have a coordinate-sorted and indexed BAM file).
- Scripts packaged with AmpliconSuite-pipeline require the `numpy`, `matplotlib` and `intervaltree` python packages. Those packages can be installed with `pip`, `conda` or similar.

AmpliconSuite-pipeline assumes both `samtools` and `bwa` executables are on the system path and can be directly invoked from bash without pathing to the executables. AmpliconSuite-pipeline will generate a BWA index for the reference genome if one is not yet in place. This adds >1hr to running time for the first use only when alignment is performed. Data repos with BWA index pre-generated are available [here](https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/). AmpliconSuite-pipeline will also function on coordinate-sorted CRAM files, [provided that the CRAM reference is in place](http://www.htslib.org/workflow/#:~:text=One%20of%20the%20key%20concepts,genome%20used%20to%20generate%20it.).

AmpliconSuite-pipeline has been tested with Ubuntu (16.04 and above) and CentOS 7. AmpliconSuite-pipeline's optional dependencies related to CNV calling will not work on CentOS 6.

**Note on using CNVKit**: We currently recommend using CNVKit for identification of AA seeds. CNVKit requires
`python3`. It also requires `R` version >= 3.5, which is non-standard on Ubuntu 16.04/14.04.


[//]: # (**Note on using Canvas**: If using Canvas, please make sure the Canvas reference genome files are located in the expected location for Canvas. To do this, you can follow instructions on the Canvas Github page. We also provide a script `$ install_canvas.sh [path/to/installation/directory/`,)

[//]: # (which when run from the AmpliconSuite-pipeline source directory will fetch the Canvas binary and download the `canvasdata` data repository. If installing on your own, create the canvasdata/ reference genome sudirectories in the folder with the Canvas executable. One installation dependency not mentioned explictly on the Canvas Readme is `dotnet-sdk-2.2`, which can be obtained in Ubuntu by running `sudo apt-get install dotnet-sdk-2.2`. )

## Installation

### Standalone installation
1. Clone the AmpliconSuite-pipeline git rep:

`git clone https://github.com/jluebeck/AmpliconSuite-pipeline.git`

2. [Install AmpliconArchitect](https://github.com/jluebeck/AmpliconArchitect), including the data repo and Mosek license.

3. [Install AmpliconClassifier](https://github.com/jluebeck/AmpliconClassifier). 



### AmpliconSuite-pipeline Docker 
A dockerized version of AmpliconSuite-pipeline is [available on dockerhub](https://hub.docker.com/repository/docker/jluebeck/prepareaa) or can be built using the Dockerfile in the `docker/` folder. It will install bwa, CNVKit and AmpliconArchitect inside the docker image. Running this docker image can be done as follows:

1. Docker:
    * Install docker: `https://docs.docker.com/install/`
    * (Optional): Add user to the docker group and relogin:
        `sudo usermod -a -G docker $USER`
   
2. License for Mosek optimization tool:
    * Obtain license file `mosek.lic` (`https://www.mosek.com/products/academic-licenses/` or `https://www.mosek.com/try/`). The license is free for academic use:
    * `mkdir $HOME/mosek`
    * After registering for a Mosek license, download license file `mosek.lic` and place it in the directory `$HOME/mosek/`.
    * If you are not able to place the license in `$HOME/mosek` you can set a custom location by exporting the bash variable `MOSEKLM_LICENSE_FILE=/custom/path/`.
   
3. Download AA data repositories and set environment variable AA_DATA_REPO:
    * Download [here](https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/) to download data repos with (`_indexed`) or
    without the bwa reference index included.
    * Set enviroment variable AA_DATA_REPO to point to the data_repo directory:
        ```bash
        mkdir data_repo && cd data_repo
        wget [url of reference_build]
        tar zxf [reference_build].tar.gz
        echo export AA_DATA_REPO=$PWD/ >> ~/.bashrc
        touch coverage.stats && chmod a+rw coverage.stats
        source ~/.bashrc
        ```
#### Obtain AmpliconSuite-pipeline image and execution script:
1. Clone GitHub repository to access the runscript
    * `git clone https://github.com/jluebeck/AmpliconSuite-pipeline.git`

2. Run the script `run_paa_docker.py` located in `AmpliconSuite-pipeline/docker`. It uses (most of) the same command line arguments one would pass to `PrepareAA.py`. CNV calling with CNVKit is integrated into the docker image (with help from Owen Chapman).

An example docker command might look like:

`AmpliconSuite-pipeline/docker/run_paa_docker.py -o /path/to/output_dir -s name_of_run -t 8 --bam bamfile.bam --run_AA --run_AC`

**You can opt to run the docker image as your current user by setting `--run_as_user`.** 

### AmpliconSuite-pipeline on Nextflow.
AmpliconSuite-pipeline can also be run through Nextflow, using the [nf-core/circdna pipeline](https://nf-co.re/circdna) constructed by [Daniel Schreyer](https://github.com/DSchreyer).

## Usage
**The main driver script for the pipeline is called `PrepareAA.py`.** Example AmpliconSuite-pipeline commands are given below.

#### Example 1: Starting from .fastq files, using CNVkit for seed generation.
`
/path/to/AmpliconSuite-pipeline/PrepareAA.py -s sample_name  -t number_of_threads --cnvkit_dir /path/to/cnvkit.py --fastqs sample_r1.fq.gz sample_r2.fq.gz --ref hg38 [--run_AA] [--run_AC]
`

`--run_AA` will invoke AmpliconArchitect directly at the end of the data preparation.
`--run_AC` will invoke AmpliconClassifier on the AmpliconArchitect outputs.

#### Example 2: Starting from sorted .bam, using CNVkit for seed generation
`
/path/to/AmpliconSuite-pipeline/PrepareAA.py -s sample_name -t n_threads --cnvkit_dir /path/to/cnvkit.py --bam sample.bam [--run_AA] [--run_AC]
`

##### Example 3: Starting from BAM and your own CNV calls (or recycled AA_CNV_SEEDS.bed)
* If using your own CNV calls:

`
/path/to/AmpliconSuite-pipeline/PrepareAA.py -s sample_name -t number_of_threads --cnv_bed your_cnvs.bed (--fastqs sample_r1.fq sample_r2.fq | --bam sample.bam) [--run_AA] [--run_AC]
`

Where the CNV bed file is formatted as (**without a header present**):

`chr    start        end       copy_number`

Additional fields between `end` and `copy_number` may exist, but `copy_number` must always be the last column.

* Note: You can also use a CNVKit .cns file instead of .bed for this argument.

* Note: CNVkit requires R version 3.5 or greater. This is not standard on older Linux systems. Specify `--rscript_path /path/to/Rscript` with your locally installed current R version if needed. 

#### Example 4: Analyzing an oncoviral sample for human-viral hybrid ecDNA detection
Note that users must start with fastq files so that the reads can also be aligned to viral genomes. CNVKit must be used for this mode.

`
/path/to/AmpliconSuite-pipeline/PrepareAA.py -s sample_name  -t n_threads --cnvkit_dir /path/to/cnvkit.py --fastqs sample_r1.fq.gz sample_r2.fq.gz --ref GRCh38_viral --cnsize_min 10000 [--run_AA] [--run_AC]
`

#### Example 5: Starting from completed AA results
If the user has one or more AA results directories inside a directory, the user can use AmpliconSuite-pipeline to call AmpliconClassifier with default settings.

`
/path/to/AmpliconSuite-pipeline/PrepareAA.py -s project_name --completed_AA_runs /path/to/location_of_all_AA_results/ --completed_run_metadata run_metadata_file.json -t 1 --ref hg38
`

Note that when this mode is used all AA results must have been generated with respect to the same reference genome version.

## Command line arguments to AmpliconSuite-pipeline

- `-o | --output_directory [outdir]`: (Optional) Directory where results will be stored. Defaults to current directory.

- `-s | --sample_name [sname]`: (Required) A name for the sample being run.

- `-t | --nthreads [numthreads]`: (Required) Number of threads to use for BWA and freebayes. Recommend 12 or more threads to be used.

- `--bam | --sorted_bam [sample.cs.bam]` **OR** `--fastqs [sample_r1.fq[.gz] sample_r2.fq[.gz]]` (Required) Input files.
Two fastqs (r1 & r2) or a coordinate sorted bam **OR** `--completed_AA_runs [/path/to/some/AA_outputs]`, a directory with 
 AA output files (one or more samples).

[//]: # (- `--canvas_dir [/path/to/Canvas_files/]` &#40;Required if not `--reuse_canvas` and not `--cnv_bed [cnvfile.bed]` and not `--cnvkit_dir`&#41; Path to directory containing the Canvas executable and `canvasdata/` subdirectory.)

- `--cnvkit_dir [/path/to/cnvkit.py]` (Required if not `--cnv_bed [cnvfile.bed]`) Path to directory containing `cnvkit.py`.

- `--completed_run_metadata`, (Required if startng with completed results). Specify a run metadata file for previously generated AA results. If you do not have it, set to 'None'." 

- `--rscript_path [/path/to/Rscript]` (Required if system Rscript version < 3.5 and using `--cnvkit_dir`). Specify a path to a local installation of Rscript compatible with CNVkit.

- `--python3_path` (Optional) Specify custom path to python3, if needed when using CNVKit (which requires python3).

- `--aa_python_interpreter` (Optional) By default PrepareAA will use the system's default `python` path. If you would like to use a different python version with AA, set this to either the path to the interpreter or `python3` or `python2` (default `python`)

- `--freebayes_dir` (Currently deprecated) Specify custom path to freebayes installation folder (not path to executable). Assumes freebayes on system path if not set. Please note this flag is currently deprecated.

- `--run_AA`: (Optional) Run AA at the end of the preparation pipeline.

- `--run_AC`: (Optional) Run AmpliconClassifier following AA. No effect if `--run_AA` not set.

- `--ref `: Name of ref genome version ("hg19","GRCh37","GRCh38","GRCh38_viral","mm10","GRCm38"). This will be auto-detected if it is not set.

- `--vcf [your_file.vcf]`: (Currently deprecated). Supply your own VCF to skip the freebayes step. 

- `--cngain [float]`: (Optional) Set a custom threshold for the CN gain considered by AA. Default: 4.5.

- `--cnsize_min [int]`: (Optional) Set a custom threshold for CN interval size considered by AA. Default: 50000.

- `--downsample [float]`: (Optional) Set a custom threshold for bam coverage downsampling during AA. Does not affect coverage in analyses outside of AA. Default: 10.

- `--use_old_samtools`: (Optional) Set this flag if your Samtools version is < 1.0. Default: False.

[//]: # (- `--reuse_canvas` &#40;Optional&#41; Reuse the Canvas results from a previous run. Default: False)

- `--cnv_bed [cnvfile.bed]` (Optional) Supply your own CNV calls, bypasses freebayes and Canvas steps. Bed file with CN estimate in last column or CNVKit .cns file.

- `--no_filter`: (Optional) Do not invoke `amplified_intervals.py` to filter amplified seed regions based on CN, size and ignorefile regions.

- `--no_QC`, (Optional) Skip QC on the BAM file.

- `--sample_metadata`, (Optional) Path to a JSON of sample metadata to build on. See template `sample_metadata_skeleton.json` for example.

- `--normal_bam [matched_normal.bam]` (Optional) Specify a matched normal BAM file for CNVKit. Not used by AA itself.

- `--purity [float between 0 and 1]` (Optional) Specify a tumor purity estimate for CNVKit. Not used by AA itself. 
  Note that specifying low purity may lead to many high copy number seed regions after rescaling is applied consider 
  setting a higher `--cn_gain` threshold for low purity samples undergoing correction.

- `--ploidy [int]` (Optional) Specify a ploidy estimate of the genome for CNVKit. Not used by AA itself.

- `--use_CN_prefilter` (Optional) Pre-filter CNV calls on number of copies gained above median chromosome arm CN. Strongly
recommended if input CNV calls have been scaled by purity or ploidy. This argument is off by default but is automatically set if `--ploidy`
or `--purity` is provided for CNVKit. 

- `--cnvkit_segmentation` Segmentation method for CNVKit (if used), defaults to CNVKit "
                        "default segmentation method (cbs).", choices=['cbs', 'haar', 'hmm', 'hmm-tumor',
                        'hmm-germline', 'none'] 

- `--AA_runmode [FULL, BPGRAPH, CYCLES, SVVIEW]` (Optional, default `FULL`). See AA documentation for more info.

- `--AA_extendmode [EXPLORE/CLUSTERED/UNCLUSTERED/VIRAL]` (Optional, default `EXPLORE`). See AA documentation for more info.

- `-AA_insert_sdevs [float]` (Optional, default 3.0) See AA documentation for more info.
 

## FAQ
Please check out the [guide document](https://github.com/jluebeck/PrepareAA/blob/master/GUIDE.md).

## Citing
If using AmpliconSuite-pipeline in your publication, please cite the [AmpliconArchitect article](https://www.nature.com/articles/s41467-018-08200-y). If using AmpliconSuite-pipeline to wrap other tools (like [CNVkit](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004873)), please cite those tools as well.


## Additional analysis tools and scripts

### - **C**andidate **AM**plicon **P**ath **E**numerato**R** `CAMPER.py`
Exahustively search an AA graph file for longest paths (cyclic and non-cyclic). A median amplicon copy number must be specified, or the script will attempt to estimate on its own.
`CAMPER.py` rescales the copy numbers by the median to estimate the multiplicity of each segment within the amplicon, and then 
searches for plausible longest paths explaining the copy number multiplicities. This is useful for identifiying some candidate ecDNA structures.
The output will be an AA-formatted cycles file with additional annotations for length and quality control filter status. The quality filters take into account root mean square residual of copy numbers ("RMSR", lower score is better), as well as "DBI" representing the Davies-Bouldin index of copy-number to multiplicity clustering. More information on the method can be found in the [methods section of this pre-print](https://www.biorxiv.org/content/10.1101/2021.11.28.470285v1).
The first entry (Cycle1) will be a cyclic path, while the second entry (Cycle2) will be a non-cyclic path. A full explanation of arguments is available with `-h`. Note that this should only be applied to AA amplicons with at most 1 ecDNA present in the AA amplicon (multiple-species reconstruction not supported).

`AmpliconSuite-pipeline/scripts/plausible_paths.py -g sample_amplicon1_graph.txt [--scaling_factor (CN estimate value)] [--remove_short_jumps] [--keep_all_LC] [--max_length (value in kbp)]`

   ### - `breakpoints_to_bed.py`
Requires `intervaltree` python package pre-installed. Write discordant edges (breakpoint junctions) from an AA graph into a pseudo-bed file.

   ### - `convert_cns_to_bed.py`
Many users will choose to run CNVKit outside of AmpliconSuite-pipeline and then want to use the CNVKit calls in AA. We recommend using the `.cns` file as a source for the seeds. 
Note the `.call.cns` file is different and contains more aggressively merged CNV calls, which we do not recommend as a source of seeds. As the `.cns` file specifies a log2 ratio,
we provide the following script to reformat the `.cns` file from CNVKit into a `.bed` file useable with AmpliconSuite-pipeline. 

Usage:
`./scripts/convert_cns_to_bed.py sample.cns`

This will output a bed file which can be fed into AmpliconSuite-pipeline. 

   ### - `cycles_to_bed.py`
Requires `intervaltree` python package pre-installed. Write an AA cycles file as a series of bed files, one for each decomposition. Segments are merged and sorted, and order and orientation of segments is lost.

### - `graph_cleaner.py`
Requires `intervaltree` python package pre-installed. Sequencing artifacts can lead to numerous spurious short breakpoint edges. This script attempts to remove edges which conform to artifactual profiles. 
Namely, very short everted (inside-out read pair) orientation edges. These will appear as numerous short brown 'spikes' in the AA amplicon image.
This script removes them from the graph file.

Usage:

`./scripts/graph_cleaner.py -g /path/to/sample_ampliconx_graph.txt [--max_hop_size 4000] `

or

`./scripts/graph_cleaner.py --graph_list /path/to/list_of_graphfiles.txt [--max_hop_size 4000] `


This will output an AA graph file(s) `/path/to/my_sample_ampliconX_cleaned_graph.txt`.

### - `graph_to_bed.py`
Requires `intervaltree` python package pre-installed. Create a bed file of the graph segments and a bedpe file of the disordant graph edges. Can also filter to only get segments with CN above `--min_cn`. 
Setting `--unmerged` will not merge adjacent graph segments and will print the graph segment CN in the last column.

Usage:

`./scripts/graph_to_bed.py -g sample_amplicon_graph.txt [--unmerged] [--min_cn 0] [--add_chr_tag]`


### - `bfb_foldback_detection.py [deprecated]`
**This script is deprecated and no longer supported, but available for legacy purposes. For more robust BFB detection, please try out [AmpliconClassifier](https://github.com/jluebeck/AmpliconClassifier).**

Requires `intervaltree` python package pre-installed. Script can be used to detect possible BFB-like signatures from AA graph files (documentation below).

To use the `bfb_foldback_detection.py` script on AA output, please create a two column file with the name of the graph file in column 1 and the path to the graph file in column 2. The rest of the command-line arguments are as follows.


##### Required arguments for running on AA results
-  `--exclude [path to $AA_DATA_REPO/[ref]/[mappability excludable file]`

- `-o [output filename prefix]`

- `--ref [hg19, GRCh37, GRCh38]`

- `--AA_graph_list [two-column file listing AA graphs]`

