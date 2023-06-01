# AmpliconSuite-pipeline
![GitHub](https://img.shields.io/github/license/jluebeck/AmpliconSuite-pipeline)

A multithread-enabled end-to-end wrapper for [AmpliconArchitect](https://github.com/jluebeck/AmpliconArchitect) and the associate tools for data preparation and interpretation. 
Performs preliminary steps (alignment, seed detection, & seed filtering) required prior to running AmpliconArchitect. AmpliconSuite-pipeline can be invoked to begin at any intermediate stage of the data preparation process and can itself invoke both AmpliconArchitect and the downstream tool AmpliconClassifier, which is used to classify ecDNAs and BFB. AmpliconSuite-pipeline was formerly called "PrepareAA".

AmpliconSuite-pipeline supports hg19, GRCh37, GRCh38 (hg38), and mouse genome mm10 (GRCm38). The tool also supports analysis with a human-viral hybrid reference genome we provide, "GRCh38_viral", which can be used to detect oncoviral hybrid focal amplifications and ecDNA in cancers with oncoviral infections.

**Current version: 0.1546.0**

[comment]: # (Versioning based on major_version.days_since_initial_commit.minor_version. Initial commit: March 5th, 2019)

We recommend browsing our [**detailed guide**](https://github.com/jluebeck/PrepareAA/blob/master/GUIDE.md) to learn about best practices and to see some FAQs.


[//]: # (**Note on using Canvas**: If using Canvas, please make sure the Canvas reference genome files are located in the expected location for Canvas. To do this, you can follow instructions on the Canvas Github page. We also provide a script `$ install_canvas.sh [path/to/installation/directory/`,)

[//]: # (which when run from the AmpliconSuite-pipeline source directory will fetch the Canvas binary and download the `canvasdata` data repository. If installing on your own, create the canvasdata/ reference genome sudirectories in the folder with the Canvas executable. One installation dependency not mentioned explictly on the Canvas Readme is `dotnet-sdk-2.2`, which can be obtained in Ubuntu by running `sudo apt-get install dotnet-sdk-2.2`. )

## Installation

### Option A: Installation-free methods
The most convenient option, however it is not suitable for analysis of large collections of samples or protected health information (PHI), and may not support more advanced command-line options. An excellent option for most users with small numbers of non-PHI samples.

#### 1. GenePattern Web Interface:
AmpliconSuite-pipeline can be run using the web interface at [GenePatter Web Interface](https://genepattern.ucsd.edu/gp). Simply search the module list for "AmpliconSuite." 
This module was constructed in collaboration with members of the GenePattern team (Edwin Huang, Ted Liefeld, Michael Reich). 

#### 2. AmpliconSuite-pipeline on Nextflow:
AmpliconSuite-pipeline can also be run through Nextflow, using the [nf-core/circdna pipeline](https://nf-co.re/circdna) constructed by [Daniel Schreyer](https://github.com/DSchreyer).


### Option B: `conda install ampliconsuite`


### Option C: Singularity & Docker images 
Containerized versions of AmpliconSuite-pipeline are available for Singularity and Docker.

A dockerized version of AmpliconSuite-pipeline is [available on dockerhub](https://hub.docker.com/repository/docker/jluebeck/prepareaa) or can be built using the Dockerfile in the `docker/` folder. It will install bwa, CNVkit and AmpliconArchitect inside the docker image. Running this docker image can be done as follows:


1. Install the container
  - Option A) Singularity:
    * Singularity installation: https://docs.sylabs.io/guides/3.0/user-guide/installation.html
    * Must have Singularity version 3.6 or higher.
    * Pull the singularity image: `singularity pull library://jluebeck/ampliconsuite-pipeline/ampliconsuite-pipeline`
    

  - Option B) Docker:
    * Docker installation: https://docs.docker.com/install/
    * Pull the docker image: `docker pull jluebeck/prepareaa`
    
    * (Optional): Add user to the docker group (log out and in after performing):
        `sudo usermod -a -G docker $USER`
   
2. License for Mosek optimization tool:
    * Obtain license file `mosek.lic` (`https://www.mosek.com/products/academic-licenses/` or `https://www.mosek.com/try/`). The license is free for academic use:
    * `mkdir $HOME/mosek`
    * After registering for a Mosek license, download license file `mosek.lic` and place it in the directory `$HOME/mosek/`.
    * If you are not able to place the license in `$HOME/mosek` you can set a custom location by exporting the bash variable `MOSEKLM_LICENSE_FILE=/custom/path/`.
   
3. Download AA data repositories and set environment variable AA_DATA_REPO:
   1. Go [here](https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/) to locate data repo(s) of your choice and make note of the URL you want.
   2. `wget`and set a bash environment variable AA_DATA_REPO to point to the data_repo directory:
       ```bash
       mkdir data_repo && cd data_repo
       wget [url of reference_build]
       tar zxf [reference_build].tar.gz
       # command below exports a bash variable which is the parent directory of the individual data repos
       echo export AA_DATA_REPO=$PWD/ >> ~/.bashrc 
       touch coverage.stats && chmod a+rw coverage.stats
       source ~/.bashrc
       ```
#### Obtain AmpliconSuite-pipeline image and execution script:
1. Clone GitHub repository to access the runscript
    * `git clone https://github.com/jluebeck/AmpliconSuite-pipeline.git`

2. Invoke the runscript to launch the container. These scripts use most of the same arguments are the main driver script `PrepareAA.py`
   - Option A) Singularity: `AmpliconSuite-pipeline/singularity/run_paa_singularity.py`
   - Option B) Docker: `AmpliconSuite-pipeline/docker/run_paa_docker.py`.
     * You can opt to run the docker image as your current user (instead of root) by setting `--run_as_user`. 


An example command might look like:

`AmpliconSuite-pipeline/singularity/run_paa_singularity.py -o /path/to/output_dir -s name_of_run -t 8 --bam bamfile.bam --run_AA --run_AC`

### Option D: Standalone installation
1. Clone the AmpliconSuite-pipeline git rep:

`git clone https://github.com/jluebeck/AmpliconSuite-pipeline.git`

2. Install other prerequisites from the section below.

## Prerequisites for standalone installation:
AmpliconSuite-pipeline supports both `python2` and `python3`, however CNVkit requires `python3`. `Python3` support for AmpliconArchitect was added in version 1.3. 

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

**Note on using CNVkit**: We currently recommend using CNVkit for identification of AA seeds. CNVkit requires
`python3`. It also requires `R` version >= 3.5, which is non-standard on Ubuntu 16.04/14.04.

## Usage
The main driver script for the standalone pipeline is called `PrepareAA.py`. 

#### Example 1: Starting from .fastq files, using CNVkit for seed generation.

>`/path/to/AmpliconSuite-pipeline/PrepareAA.py -s sample_name  -t number_of_threads --cnvkit_dir /path/to/cnvkit.py --fastqs sample_r1.fq.gz sample_r2.fq.gz --ref hg38 [--run_AA] [--run_AC]`


`--run_AA` will invoke AmpliconArchitect directly at the end of the data preparation.
`--run_AC` will invoke AmpliconClassifier on the AmpliconArchitect outputs.

#### Example 2: Starting from .bam, using CNVkit for seed generation

>`/path/to/AmpliconSuite-pipeline/PrepareAA.py -s sample_name -t n_threads [--cnvkit_dir /path/to/cnvkit.py] --bam sample.bam [--run_AA] [--run_AC]`

`--cnvkit_dir` is only needed if cnvkit.py is not on the system path (typically if it was a custom install).

#### Example 3: Starting from .bam and your own whole-genome CNV calls, or an existing AA_CNV_SEEDS.bed
* If using your own CNV calls:


>`/path/to/AmpliconSuite-pipeline/PrepareAA.py -s sample_name -t number_of_threads --cnv_bed your_cnvs.bed --bam sample.bam [--run_AA] [--run_AC]`

Where the CNV bed file reports the following four fields:

`chr    start        end       copy_number`

Additional fields between `end` and `copy_number` may exist, but `copy_number` must always be the last column.

* You can also use the CNVkit `sample_name.cns` file instead of .bed for this argument.

* CNVkit requires R version 3.5 or greater. This is not standard on older Linux systems. Specify `--rscript_path /path/to/Rscript` with your locally installed current R version if needed. 

#### Example 4: Analyzing a collection of related samples (same origin)

Have multiple samples from the same patient, cell line, etc.? These should be run as a group to ensure that the same regions are studied across samples.
Please see the `GroupedAnalysis.py` [example below](#--grouped-analysis-of-related-samples-groupedanalysispy) for instructions.

#### Example 5: Analyzing an oncoviral sample for human-viral hybrid ecDNA detection
Note that users must start with fastq files and `--ref GRCh38_viral` or a bam file aligned to the `AA_DATA_REPO/GRCh38_viral` reference.


>`/path/to/AmpliconSuite-pipeline/PrepareAA.py -s sample_name  -t n_threads --fastqs sample_r1.fq.gz sample_r2.fq.gz --ref GRCh38_viral --cnsize_min 10000 [--run_AA] [--run_AC]`

#### Example 6: Starting from completed AA results
If the user has one or more AA results directories inside a directory, the user can use AmpliconSuite-pipeline to call AmpliconClassifier with default settings.


>`/path/to/AmpliconSuite-pipeline/PrepareAA.py -s project_name --completed_AA_runs /path/to/location_of_all_AA_results/ --completed_run_metadata run_metadata_file.json -t 1 --ref hg38`

Note that when this mode is used all AA results must have been generated with respect to the same reference genome version.

## Command line arguments to AmpliconSuite-pipeline

#### Required
- `-o | --output_directory {outdir}`: (Optional) Directory where results will be stored. Defaults to current directory.

- `-s | --sample_name {sname}`: (Required) A name for the sample being run.

- `-t | --nthreads {int}`: (Required) Number of threads to use for BWA and CNVkit. Recommend 12 or more threads to be used.

- One of the following input files:

  * `--bam | --sorted_bam {sample.cs.bam}` Coordinate-sorted bam
  * `--fastqs {sample_r1.fq.gz sample_r2.fq.gz}` Two fastqs (r1 & r2)
  * `--completed_AA_runs {/path/to/some/AA_outputs}`, a directory with AA output files (one or more samples).

#### Optional

- `--cnv_bed {cnvfile.bed}` Supply your own CNV calls. Bed file with CN estimate in last column, or the CNVkit `sample.cns` file. If not specified, CNVkit will be called by the wrapper.

- `--cnvkit_dir {/path/to/cnvkit.py}` Path to directory containing `cnvkit.py`. Required if CNVkit was installed from source and `--cnv_bed [cnvfile.bed]` is not given.

- `--run_AA`: Run AA at the end of the preparation pipeline.

- `--run_AC`: Run AmpliconClassifier following AA. No effect if `--run_AA` not set.

- `--normal_bam {matched_normal.bam}` Specify a matched normal BAM file for CNVkit. Not used by AA itself.

- `--completed_run_metadata {run_metadata.json}`, Required only if starting with completed results (`--completed_AA_runs`). Specify a run metadata file for previously generated AA results. If you do not have it, set to 'None'." 

- `--rscript_path {/path/to/Rscript}` (Relevant if using CNVkit and system Rscript version is < 3.5). Specify a path to a local installation of Rscript.

- `--python3_path {/path/to/python3}` Specify custom path to python3 if needed when using CNVkit.

- `--aa_python_interpreter {/path/to/python}` By default PrepareAA will use the system's default `python` path. If you would like to use a different python version with AA, set this to either the path to the interpreter or `python3` or `python2` (default `python`)

[//]: # (- `--freebayes_dir` &#40;Currently deprecated&#41; Specify custom path to freebayes installation folder &#40;not path to executable&#41;. Assumes freebayes on system path if not set. Please note this flag is currently deprecated.)
- `--ref {ref name} `: Name of ref genome version, one of `"hg19","GRCh37","GRCh38","GRCh38_viral","mm10","GRCm38"`. This will be auto-detected if it is not set.

[//]: # (- `--vcf [your_file.vcf]`: &#40;Currently deprecated&#41;. Supply your own VCF to skip the freebayes step. )
- `--cngain {float}`: Set a custom threshold for the CN gain considered by AA. Default: 4.5.

- `--cnsize_min {int}`: Set a custom threshold for CN interval size considered by AA. Default: 50000.

- `--downsample {float}`: Set a custom threshold for bam coverage downsampling during AA. Does not affect coverage in analyses outside of AA. Default: 10.

- `--use_old_samtools`: Set this flag if your SAMtools version is < 1.0.

- `--no_filter`: Do not invoke `amplified_intervals.py` to filter amplified seed regions based on CN, size and ignorefile regions.

- `--no_QC`: Skip QC on the BAM file.

- `--sample_metadata {sample_metadata.json}`: Path to a JSON of sample metadata to build on. Please expand from the template `sample_metadata_skeleton.json`.

- `--purity {float between 0 and 1}`: Specify a tumor purity estimate for CNVkit (not used by AA). 
  Note that specifying low purity may lead to many high copy-number seed regions after rescaling is applied. Consider 
  setting a higher `--cngain` threshold for low purity samples undergoing correction (e.g. `--cngain 8`).

- `--ploidy {float}`: Specify a ploidy estimate of the genome for CNVkit (not used by AA).

- `--cnvkit_segmentation {str}`: Segmentation method for CNVkit (if used), defaults to `cbs`., choices=`'cbs', 'haar', 'hmm', 'hmm-tumor',
                        'hmm-germline', 'none'` 

- `--AA_runmode {FULL, BPGRAPH, CYCLES, SVVIEW}`: Default `FULL`. See AA documentation for more info.

- `--AA_extendmode {EXPLORE/CLUSTERED/UNCLUSTERED/VIRAL}`: Default `EXPLORE`. See AA documentation for more info.

- `--AA_insert_sdevs {float}`: Default 3.0. Suggest raising to 8 or 9 if library has poorly-controlled insert size (low fraction of properly-paired reads). See AA documentation for more info.

- `--samtools_path`: Path to a specific samtools binary for use (e.g., /path/to/my/samtools). Uses samtools on system path by default.
 

## FAQ
Please check out our [guide document](https://github.com/jluebeck/PrepareAA/blob/master/GUIDE.md).

## Packaging outputs for AmpliconRepository
We will soon release an online platform for storing and sharing your AmpliconSuite-pipeline outputs.

To package a collection of AA outputs for AmpliconRepository, you will need to do the following steps.
1. (Recommended) Before running, using the file `sample_metadata_skeleton.json` as a template, please create a copy of the file for each sample, and fill out the JSON file. Provide this to `PrepareAA.py` using `--sample_metadata {sample_metadata.json}`
2. Run the `make_results_table.py` script from AmpliconClassifier on your AC outputs.
 * `cd [directory of classification files]`
 * `$AC_SRC/make_results_table.py -i samples.input --summary_map samples_summary_map.txt --classification_file samples_amplicon_classification_profiles.tsv --ref [hg19/hg38/...] {--sample_metada_file/--sample_metadata_list}`
   * The .input and \_summary_map.txt files are created by `make_input.sh`.
3. Create a tar.gz file from your AA outputs `tar -czf my_collection.tar.gz /path/to/AA_outputs/` (creating a `.zip` also works)
4. If you have not already, create an account at [GenePattern](https://genepattern.ucsd.edu/gp).
5. Upload your compressed collection of AA output files (one or more samples) to the `AmpliconSuiteAggregator` GenePattern module.
 * Go to https://genepattern.ucsd.edu/gp, and sign in.
 * In the top-left search bar, search for `AmpliconSuiteAggregator`
6. Run the aggregator and download the aggregated `.tar.gz` result file.
7. Upload to AmpliconRepository (coming soon).


## Citing
If using AmpliconSuite-pipeline in your publication, please cite the modules used in the analysis, which are summarized in [CITATIONS.md](https://github.com/jluebeck/AmpliconSuite-pipeline/blob/master/CITATIONS.md).


## Additional analysis tools and scripts

### - Grouped analysis of related samples `GroupedAnalysis.py`
For samples derived from a common origin (longitudinal, multiregional sampling from the same source material), it is advised that the seed intervals be unified before running AA in order to provide the best comparability
between runs. We provide a script `GroupedAnalysis.py` which automates this analysis. `GroupedAnalysis.py` takes almost all the same arguments as `PrepareAA.py`, 
however it requires an additional input file, listing the inputs. This file
is to be formatted as follows

`sample_name` `bamfile` `"tumor"/"normal"` `[CNV_calls]` `[sample_metadata_json]`

Where `CNV_calls` and `sample_metadata_json` are optional. However, they are positional, so if `CNV_calls` is skipped, it should be set as either `NA` or `None`.

AA and AC will be run by default, but can be disabled with `--no_AA`.

Example command:

> `/path/to/PrepareAA/GroupedAnalysis.py -i {inputs.txt} -o {output_dir} -t {num_threads}`

### - **C**andidate **AM**plicon **P**ath **E**numerato**R** `CAMPER.py`
Exahustively search an AA graph file for longest paths (cyclic and non-cyclic). A median amplicon copy number must be specified, or the script will attempt to estimate on its own.
`CAMPER.py` rescales the copy numbers by the median to estimate the multiplicity of each segment within the amplicon, and then 
searches for plausible longest paths explaining the copy number multiplicities. This is useful for identifiying some candidate ecDNA structures.
The output will be an AA-formatted cycles file with additional annotations for length and quality control filter status.
The quality filters take into account root mean square residual of copy numbers ("RMSR", lower score is better), as well as "DBI" representing the Davies-Bouldin index of copy-number to multiplicity clustering. More information on the method can be found in the [methods section of this publication](https://www.nature.com/articles/s41588-022-01190-0).
The first entry (Cycle1) will be a cyclic path, while the second entry (Cycle2) will be a non-cyclic path. A full explanation of arguments is available with `-h`. Note that this should only be applied to AA amplicons with at most 1 ecDNA present in the AA amplicon (multiple-species reconstruction not supported).

> `AmpliconSuite-pipeline/scripts/plausible_paths.py -g sample_amplicon1_graph.txt [--scaling_factor (CN estimate value)] [--remove_short_jumps] [--keep_all_LC] [--max_length (value in kbp)]`

### - `breakpoints_to_bed.py`
Requires `intervaltree` python package pre-installed. Write discordant edges (breakpoint junctions) from an AA graph into a pseudo-bed file.

### - `convert_cns_to_bed.py`
Many users will choose to run CNVkit outside of AmpliconSuite-pipeline and then want to use the CNVkit calls in AA. We recommend using the `.cns` file as a source for the seeds. 
Note the `.call.cns` file is different and contains more aggressively merged CNV calls, which we do not recommend as a source of seeds. As the `.cns` file specifies a log2 ratio,
we provide the following script to reformat the `.cns` file from CNVkit into a `.bed` file useable with AmpliconSuite-pipeline. 

Usage:
>`scripts/convert_cns_to_bed.py sample.cns`

This will output a bed file which can be fed into AmpliconSuite-pipeline. 

### - `cycles_to_bed.py`
Requires `intervaltree` python package pre-installed. Write an AA cycles file as a series of bed files, one for each decomposition. Segments are merged and sorted, and order and orientation of segments is lost.

### - `graph_cleaner.py`
Requires `intervaltree` python package pre-installed. Sequencing artifacts can lead to numerous spurious short breakpoint edges. This script attempts to remove edges which conform to artifactual profiles. 
Namely, very short everted (inside-out read pair) orientation edges. These will appear as numerous short brown 'spikes' in the AA amplicon image.
This script removes them from the graph file.

Usage:

>`scripts/graph_cleaner.py -g /path/to/sample_ampliconx_graph.txt [--max_hop_size 4000] `

or

>`scripts/graph_cleaner.py --graph_list /path/to/list_of_graphfiles.txt [--max_hop_size 4000] `


This will output an AA graph file(s) `/path/to/my_sample_ampliconX_cleaned_graph.txt`.

### - `graph_to_bed.py`
Requires `intervaltree` python package pre-installed. Create a bed file of the graph segments and a bedpe file of the disordant graph edges. Can also filter to only get segments with CN above `--min_cn`. 
Setting `--unmerged` will not merge adjacent graph segments and will print the graph segment CN in the last column.

Usage:

>`scripts/graph_to_bed.py -g sample_amplicon_graph.txt [--unmerged] [--min_cn 0] [--add_chr_tag]`


### - `bfb_foldback_detection.py [deprecated]`
**This script is deprecated and no longer supported, but available for legacy purposes. For more robust BFB detection, please try out [AmpliconClassifier](https://github.com/jluebeck/AmpliconClassifier).**

Requires `intervaltree` python package pre-installed. Script can be used to detect possible BFB-like signatures from AA graph files (documentation below).

To use the `bfb_foldback_detection.py` script on AA output, please create a two column file with the name of the graph file in column 1 and the path to the graph file in column 2. The rest of the command-line arguments are as follows.


##### Required arguments for running on AA results
-  `--exclude [path to $AA_DATA_REPO/[ref]/[mappability excludable file]`

- `-o [output filename prefix]`

- `--ref [hg19, GRCh37, GRCh38]`

- `--AA_graph_list [two-column file listing AA graphs]`

