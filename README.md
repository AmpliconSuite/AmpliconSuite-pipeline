# AmpliconSuite-pipeline
![GitHub release (latest by date)](https://img.shields.io/github/v/release/AmpliconSuite/AmpliconSuite-pipeline)
[![Docker Image Version (latest by date)](https://img.shields.io/docker/v/jluebeck/prepareaa?logo=docker)](https://hub.docker.com/r/jluebeck/prepareaa)
[![Docker pulls](https://img.shields.io/docker/pulls/jluebeck/prepareaa?logo=docker)]((https://hub.docker.com/r/jluebeck/prepareaa))
[![Singularity](https://img.shields.io/badge/singularity-available-blue)](https://cloud.sylabs.io/library/jluebeck/ampliconsuite-pipeline/ampliconsuite-pipeline)
[![Conda](https://img.shields.io/conda/dn/bioconda/ampliconsuite?logo=Anaconda)](https://anaconda.org/bioconda/ampliconsuite)


A multithread-enabled end-to-end wrapper for [AmpliconArchitect](https://github.com/jluebeck/AmpliconArchitect) and AmpliconClassifier to enable analysis of focal copy number amplifications such as ecDNA or BFBs from paired-end whole genome sequencing data.

AmpliconSuite-pipeline can be invoked to begin at any intermediate stage of the data preparation process and can itself invoke both AmpliconArchitect and the downstream tool AmpliconClassifier. AmpliconSuite-pipeline was formerly called "PrepareAA".

[comment]: # (Versioning based on major_version.days_since_initial_commit.minor_version. Initial commit: March 5th, 2019)

We recommend browsing our [**detailed guide**](https://github.com/AmpliconSuite/AmpliconSuite-pipeline/blob/master/documentation/GUIDE.md) to learn about best practices and to see some FAQs. 

AmpliconSuite-pipeline supports hg19, GRCh37, GRCh38 (hg38), and mouse genome mm10 (GRCm38). It also supports analysis with a human-viral hybrid reference genome we provide, "GRCh38_viral", which can be used to detect oncoviral hybrid focal amplifications in oncoviral cancers.

## Licenses
The modules wrapped in AmpliconSuite-pipeline use the following licenses. Please note that the AmpliconArchitect license specifies that AmpliconArchitect is for research use and does not give license for commerical for-profit use.
- [AmpliconSuite-pipeline license](https://github.com/AmpliconSuite/AmpliconSuite-pipeline/blob/master/LICENSE) (BSD 2-Clause)
- [AmpliconArchitect license](https://github.com/AmpliconSuite/AmpliconArchitect) (University of California software license)
- [AmpliconClassifier license](https://github.com/AmpliconSuite/AmpliconClassifier/blob/main/LICENSE) (BSD 2-Clause)

Other dependencies used by these modules (e.g. Mosek, samtools, etc.) have their own set of licensing requirements which users should make themselves aware of as needed. The Mosek license requires that users obtain a copy (which is free for academic use) from the Mosek website. More information is available in the installation section.

## Installation

### Option A: Installation-free methods
The most convenient option, however it is not suitable for analysis of large collections of samples or protected health information (PHI), and may not support more advanced command-line options. An excellent option for most users with small numbers of non-PHI samples.

#### GenePattern Web Interface:
AmpliconSuite-pipeline can be run using the web interface at [GenePattern Web Interface](https://genepattern.ucsd.edu/gp). Simply search the module list for "AmpliconSuite." 
This module was constructed in collaboration with members of the GenePattern team (Edwin Huang, Ted Liefeld, Michael Reich). 

#### Nextflow:
AmpliconSuite-pipeline can also be run through Nextflow, using the [nf-core/circdna pipeline](https://nf-co.re/circdna) constructed by [Daniel Schreyer](https://github.com/DSchreyer).

### Option B: Install with Conda or Mamba
```bash 
conda install -c bioconda -c mosek ampliconsuite mosek
wget https://raw.githubusercontent.com/AmpliconSuite/AmpliconSuite-pipeline/bioconda/install.sh
source install.sh --finalize_only  # this will confirm the data repo path and mosek license directory
```

If Conda fails to solve the environment, [Mamba](https://mamba.readthedocs.io/en/latest/installation.html) seems to function more robustly for installing AmpliconSuite. These steps also function on macOS.
```bash
# alternate instructions using Mamba (solves dependencies more effectively on some setups)
mamba create -n ampsuite python=3.10
mamba activate ampsuite
mamba install -c conda-forge -c bioconda -c mosek ampliconsuite mosek
wget https://raw.githubusercontent.com/AmpliconSuite/AmpliconSuite-pipeline/bioconda/install.sh
source install.sh --finalize_only
```

**then proceeed to step 2 of Option C (below) ...**

### Option C: Standalone installation using the installer script
Can be used on most modern Unix systems (e.g. Ubuntu 18.04+, CentOS 7+, macOS). Requires `python>=3.7`.
1. Pull source code and run install script (**skip if installed via Conda**):
    ```bash
    # first install some dependencies (BWA, R, samtools) if you don't already have them
    # for ubuntu:
    sudo apt install bwa r-base samtools
    # or for macOS: brew install bwa r samtools  
   
    git clone https://github.com/AmpliconSuite/AmpliconSuite-pipeline
    cd AmpliconSuite-pipeline
    # To see install  options, consider first doing 
    # source ./install -h
    # The install.sh script will install python dependencies using 'python3 -m pip install' 
    source ./install.sh 
    ```

2. **Start here if you installed with Conda or Mamba**. Populate the AA data repo with required annotations for the reference genomes of interest. 
    - See the list of available AA annotations [here](https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/). Copy the URL of the one(s) you need.
    ```bash
    cd $AA_DATA_REPO
    # go to https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/
    # copy the url of the data repo you need.
    # "_indexed" indicates the bwa index is included - only useful if starting from .fastqs.
    wget [url of reference_build]
    tar -xzf [reference_build].tar.gz
    rm [reference_build.tar.gz]
    ```

3. Obtain the Mosek optimization tool license (free for academic use) and place it in `$HOME/mosek/`. AA will not work without it.

Lastly, as a completely optional step, if you want the Arial font in your AA figures (helpful for publication-quality fonts), but do not have Arial  on your Linux system, please see [these instructions](https://github.com/AmpliconSuite/AmpliconSuite-pipeline/blob/master/documentation/CUSTOM_INSTALL.md#getting-mscorefonts-onto-your-system) for making it available to Matplotlib. 



### Option D: Singularity & Docker images 
Containerized versions of AmpliconSuite-pipeline are available for Singularity and Docker.

1. Obtain the AmpliconSuite-pipeline image from the options below:
   - **Singularity**:
     * Singularity installation: https://docs.sylabs.io/guides/3.0/user-guide/installation.html
     * Must have Singularity version 3.6 or later.
     * Pull the singularity image: `singularity pull library://jluebeck/ampliconsuite-pipeline/ampliconsuite-pipeline`

   - **Docker**:
     * Docker installation: https://docs.docker.com/install/
     * Pull the docker image: `docker pull jluebeck/prepareaa`
    
     * (Optional): Add user to the docker group:
         `sudo usermod -a -G docker $USER` (log out and back in after performing)
   
2. Obtain the execution script and configure the data repo location
    ```bash
    git clone https://github.com/AmpliconSuite/AmpliconSuite-pipeline
    cd AmpliconSuite-pipeline
   # Can use ./install.sh -h to see help before installing
    source ./install.sh --finalize_only
    ```

3. License for Mosek optimization tool:
    * Obtain license file `mosek.lic` (`https://www.mosek.com/products/academic-licenses/`). The license is free for academic use.
    * Place the file in `$HOME/mosek/` (i.e, the `mosek/` folder that now exists in your home directory).
    * If you are not able to place the license in the default location, you can set a custom location by exporting the bash variable `MOSEKLM_LICENSE_FILE=/custom/path/`.

   
4. Download AA data repositories and set environment variable AA_DATA_REPO:
   - Go [here](https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/) to locate data repo(s) of your choice and make note of the URL you want.
      ```bash
      cd $AA_DATA_REPO
      wget [url of reference_build]
      tar zxf [reference_build].tar.gz
      rm [reference_build].tar.gz
      ```
   - If you do not do this process the container runscript attempt to download the files itself before launching the container.

#### Launching the execution script for the container:

These scripts use most of the same arguments are the main driver script `AmpliconSuite-pipeline.py`
   - Singularity: `AmpliconSuite-pipeline/singularity/run_paa_singularity.py`
   - Docker: `AmpliconSuite-pipeline/docker/run_paa_docker.py`.
     * You can opt to run the docker image as your current user (instead of root) by setting `--run_as_user`. 


An example command might look like:

`AmpliconSuite-pipeline/singularity/run_paa_singularity.py -o /path/to/output_dir -s name_of_run -t 8 --bam bamfile.bam --run_AA --run_AC`

### Option E: Standalone installation without installer script
Try this if you are going to use `python2`. Please see [the instructions here](documentation/CUSTOM_INSTALL.md).


## Running AmpliconSuite-pipeline
The main driver script for the standalone pipeline is called `AmpliconSuite-pipeline.py`. 

#### Example 1: Starting from .fastq files, using CNVkit for seed generation.

>`AmpliconSuite-pipeline.py -s sample_name  -t number_of_threads --cnvkit_dir /path/to/cnvkit.py --fastqs sample_r1.fq.gz sample_r2.fq.gz --ref hg38 [--run_AA] [--run_AC]`


`--run_AA` will invoke AmpliconArchitect directly at the end of the data preparation.
`--run_AC` will invoke AmpliconClassifier on the AmpliconArchitect outputs.

#### Example 2: Starting from .bam, using CNVkit for seed generation

>`AmpliconSuite-pipeline.py -s sample_name -t n_threads [--cnvkit_dir /path/to/cnvkit.py] --bam sample.bam [--run_AA] [--run_AC]`

`--cnvkit_dir` is only needed if cnvkit.py is not on the system path (typically if it was a custom install).

#### Example 3: Starting from .bam and your own whole-genome CNV calls, or an existing AA_CNV_SEEDS.bed
* If using your own CNV calls:


>`AmpliconSuite-pipeline.py -s sample_name -t number_of_threads --cnv_bed your_cnvs.bed --bam sample.bam [--run_AA] [--run_AC]`

Where the CNV bed file reports the following four fields:

`chr    start        end       copy_number`

Additional fields between `end` and `copy_number` may exist, but `copy_number` must always be the last column.

* You can also use the CNVkit `sample_name.cns` file instead of .bed for this argument.

* CNVkit requires R version 3.5 or greater. This is not standard on older Linux systems. Specify `--rscript_path /path/to/Rscript` with your locally installed current R version if needed. 

#### Example 4: Analyzing a collection of related samples (same origin)

Have multiple samples from the same patient, cell line, etc.? These should be run as a group to ensure that the same regions are studied across samples.
Please see the `GroupedAnalysisAmpSuite.py` [example below](#--grouped-analysis-of-related-samples-groupedanalysispy) for instructions.

#### Example 5: Analyzing an oncoviral sample for human-viral hybrid ecDNA detection
Note that users must start with fastq files and `--ref GRCh38_viral` or a bam file aligned to the `AA_DATA_REPO/GRCh38_viral` reference.


>`AmpliconSuite-pipeline.py -s sample_name  -t n_threads --fastqs sample_r1.fq.gz sample_r2.fq.gz --ref GRCh38_viral --cnsize_min 10000 [--run_AA] [--run_AC]`

#### Example 6: Starting from completed AA results
If the user has one or more AA results directories inside a directory, the user can use AmpliconSuite-pipeline to call AmpliconClassifier with default settings.


>`AmpliconSuite-pipeline.py -s project_name --completed_AA_runs /path/to/location_of_all_AA_results/ --completed_run_metadata run_metadata_file.json -t 1 --ref hg38`

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

- `--no_filter`: Do not invoke `amplified_intervals.py` to filter amplified seed regions based on CN, size and ignorefile regions. Skipping filtering is not recommended and will harm reliability of results.

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


## Interpreting classification outputs
- Information about the amplicon classification files produced at the end of the workflow are available [here](https://github.com/AmpliconSuite/AmpliconClassifier#3-output).
- Information on interpreting the AA cycles file is available [here](https://github.com/jluebeck/AmpliconArchitect#interpreting-the-aa-cycles-files).

## Testing your installation
To ensure your local installation of, we provide a small test dataset (~3Gb), which users can download from SRA. After 
obtaining either the BAM or FASTQ files for GBM39_FF-8 provided [on SRA](https://trace.ncbi.nlm.nih.gov/Traces/?run=SRR4009231), users 
can run the method and compare their results against the files in the `test_outputs/` directory.

## FAQ
Please check out our [guide document](documentation/GUIDE.md).

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
If using AmpliconSuite-pipeline in your publication, please cite all the relevant modules used in the analysis, which are summarized in [CITATIONS.md](https://github.com/jluebeck/AmpliconSuite-pipeline/blob/master/CITATIONS.md).

## Additional analysis tools and scripts

### - Grouped analysis of related samples `GroupedAnalysisAmpSuite.py`
For samples derived from a common origin (longitudinal, multiregional sampling from the same source material), it is advised that the seed intervals be unified before running AA in order to provide the best comparability
between runs. We provide a script `GroupedAnalysisAmpSuite.py` which automates this analysis. `GroupedAnalysisAmpSuite.py` takes almost all the same arguments as `PrepareAA.py`, 
however it requires an additional input file, listing the inputs. This file is to be formatted as follows

`sample_name` `bamfile` `"tumor"/"normal"` `[CNV_calls]` `[sample_metadata_json]`

Where `CNV_calls` and `sample_metadata_json` are optional. All samples listed in each file should be uniquely named and from the same group of related samples. Do not include different collections of related samples in the same table - make different tables. However, they are positional, so if `CNV_calls` is skipped, it should be set as either `NA` or `None`.

AA and AC will be run by default, but can be disabled with `--no_AA`.

Example command:

> `GroupedAnalysisAmpSuite.py -i {inputs.txt} -o {output_dir} -t {num_threads}`

### - **C**andidate **AM**plicon **P**ath **E**numerato**R** `CAMPER.py`
Exahustively search an AA graph file for longest paths (cyclic and non-cyclic). A median amplicon copy number must be specified, or the script will attempt to estimate on its own.
`CAMPER.py` rescales the copy numbers by the median to estimate the multiplicity of each segment within the amplicon, and then 
searches for plausible longest paths explaining the copy number multiplicities. This is useful for identifiying some candidate ecDNA structures.
The output will be an AA-formatted cycles file with additional annotations for length and quality control filter status.
The quality filters take into account root-mean-square residual of copy numbers ("RMSR", lower score is better), as well as "DBI" representing the Davies-Bouldin index of copy-number to multiplicity clustering. More information on the method can be found in the [methods section of this publication](https://www.nature.com/articles/s41588-022-01190-0).
The first entry (Cycle1) will be a cyclic path, while the second entry (Cycle2) will be a non-cyclic path. A full explanation of arguments is available with `-h`. Note that this should only be applied to AA amplicons with at most 1 ecDNA present in the AA amplicon (multiple-species reconstruction not supported).

> `AmpliconSuite-pipeline/scripts/plausible_paths.py -g sample_amplicon1_graph.txt [--scaling_factor (CN estimate value)] [--remove_short_jumps] [--keep_all_LC] [--max_length (value in kbp)]`

### - `breakpoints_to_bed.py`
Requires `intervaltree` python package pre-installed. Write discordant edges (breakpoint junctions) from an AA graph into a pseudo-bed file.

### - `convert_cns_to_bed.py`
Many users will choose to run CNVkit outside AmpliconSuite-pipeline and then want to use the CNVkit calls in AA. We recommend using the `.cns` file as a source for the seeds. 
Note the `.call.cns` file is different and contains more aggressively merged CNV calls, which we do not recommend as a source of seeds. As the `.cns` file specifies a log2 ratio,
we provide the following script to reformat the `.cns` file from CNVkit into a `.bed` file useable with AmpliconSuite-pipeline. 

Usage:
>`scripts/convert_cns_to_bed.py sample.cns`

This will output a bed file which can be fed into AmpliconSuite-pipeline. 

### - `cycles_to_bed.py`
Requires `intervaltree` python package pre-installed. Write an AA cycles file as a series of bed files, one for each decomposition. Writes two types of bed files, an `unordered_cycle` file, where segments are merged and sorted, and order and orientation of segments is lost, and also writes
an ordered file where the order and orientation of the genome segments comprising the cycle is maintained.

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

