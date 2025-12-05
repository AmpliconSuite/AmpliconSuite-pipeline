![AmpliconSuite-pipeline logo](https://github.com/AmpliconSuite/AmpliconSuite-pipeline/blob/master/images/AmpliconSuite-pipeline.png)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/AmpliconSuite/AmpliconSuite-pipeline)
[![Docker Image Version (latest by date)](https://img.shields.io/docker/v/jluebeck/ampliconsuite-pipeline?logo=docker)](https://hub.docker.com/r/jluebeck/ampliconsuite-pipeline)
[![Docker pulls](https://img.shields.io/docker/pulls/jluebeck/ampliconsuite-pipeline?logo=docker)]((https://hub.docker.com/r/jluebeck/ampliconsuite-pipeline))
[![Singularity](https://img.shields.io/badge/singularity-available-blue)](https://cloud.sylabs.io/library/jluebeck/ampliconsuite-pipeline/ampliconsuite-pipeline)
[![Conda](https://img.shields.io/conda/dn/bioconda/ampliconsuite?logo=Anaconda)](https://anaconda.org/bioconda/ampliconsuite)


A multithread-enabled end-to-end wrapper for [AmpliconArchitect](https://github.com/jluebeck/AmpliconArchitect) and AmpliconClassifier to enable analysis of focal copy number amplifications such as ecDNA or BFBs from paired-end whole genome sequencing data.

AmpliconSuite-pipeline can be invoked to begin at any intermediate stage of the data preparation process and can itself invoke both AmpliconArchitect and the downstream tool AmpliconClassifier. AmpliconSuite-pipeline was formerly called "PrepareAA".

[comment]: # (Initial commit: March 5th, 2019)

We recommend browsing our [**detailed guide**](https://github.com/AmpliconSuite/AmpliconSuite-pipeline/blob/master/documentation/GUIDE.md) to learn about best practices and to see some FAQs. 

AmpliconSuite-pipeline supports hg19, GRCh37, GRCh38 (hg38), and mouse genome mm10 (GRCm38). It also supports analysis with a human-viral hybrid reference genome we provide, "GRCh38_viral", which can be used to detect oncoviral hybrid focal amplifications in oncoviral cancers.

## Licenses
The modules wrapped in AmpliconSuite-pipeline use the following licenses. Please note that the AmpliconArchitect license specifies that AmpliconArchitect is for research use and does not give license for commerical for-profit use.
- [AmpliconSuite-pipeline license](https://github.com/AmpliconSuite/AmpliconSuite-pipeline/blob/master/LICENSE) (BSD 2-Clause)
- [AmpliconArchitect license](https://github.com/AmpliconSuite/AmpliconArchitect) (University of California software license)
- [AmpliconClassifier license](https://github.com/AmpliconSuite/AmpliconClassifier/blob/main/LICENSE) (BSD 2-Clause)

Other dependencies used by these modules (e.g. Mosek, samtools, etc.) have their own set of licensing requirements which users should make themselves aware of as needed. The Mosek license requires that users obtain a copy (which is free for academic use) from the Mosek website. More information is available in the installation section.

## Installation

### Option A: Installation-free platforms

#### [GenePattern](https://cloud.genepattern.org/gp):
AmpliconSuite-pipeline can be run using the web interface at [GenePattern Web Interface](https://cloud.genepattern.org/gp). Search the module list for `AmpliconSuite`. 
Constructed in collaboration with members of the GenePattern team (Edwin Huang, Ted Liefeld, Michael Reich). The most convenient option, but not suitable for analysis of large collections of samples or protected health information (PHI), and may not support more advanced command-line options. An excellent option for most users with small numbers of non-PHI samples.


#### [Nextflow](https://nf-co.re/circdna):
AmpliconSuite-pipeline can also be run through Nextflow, using the [nf-core/circdna pipeline](https://nf-co.re/circdna) constructed by [Daniel Schreyer](https://github.com/DSchreyer).

### Option B: Conda or Mamba
```bash 
conda create -n ampsuite && conda activate ampsuite
conda install -c bioconda -c conda-forge ampliconsuite 
conda install -c mosek mosek

# then run the installer script to finalize the locations of the data repo and mosek license 
wget https://raw.githubusercontent.com/AmpliconSuite/AmpliconSuite-pipeline/master/install.sh
install.sh --finalize_only  # -h to see options
```

Then [obtain the Mosek license](https://www.mosek.com/products/academic-licenses/) (free for academic use) and place it in `$HOME/mosek/`. AA will not work without it.

- If Conda fails to solve the environment, [Mamba](https://mamba.readthedocs.io/en/latest/installation.html) seems to function robustly for installing AmpliconSuite. These steps also function on macOS.
```bash
# alternate instructions using Mamba (solves dependencies more effectively on some setups)
mamba create -n ampsuite python=3.10 && mamba activate ampsuite
mamba install -c conda-forge -c bioconda -c mosek ampliconsuite mosek
# ... proceed with wget and install.sh steps from above
```



### Option C: Standalone installation using the installer script
Can be used on recent Unix systems (e.g. Ubuntu 18.04+, CentOS 7+, macOS). Requires `python>=3.7`.
1. Pull source code and run install script (**skip if installed via Conda**):
    ```bash
    # first install some dependencies (BWA, R, samtools) if you don't already have them
    # for ubuntu:
    sudo apt install bedtools bwa curl git r-base samtools
    # or for macOS: brew install bedtools bwa curl git r samtools  
   
    git clone https://github.com/AmpliconSuite/AmpliconSuite-pipeline
    cd AmpliconSuite-pipeline
    # To see install  options, consider first doing 
    # ./install -h
    # The install.sh script will install python dependencies using 'python3 -m pip install' 
    ./install.sh 
   source ~/.bashrc  # either do this or open a new session to make changes live
    ```
   
Mac users will need to perform one additional installation step:
```bash
brew install coreutils
```

2. [Obtain the Mosek license](https://www.mosek.com/products/academic-licenses/) (free for academic use) and place it in `$HOME/mosek/`. AA will not work without it.


3. (Optional) If you want the Arial font in your AA figures (helpful for publication-quality fonts), but do not have Arial  on your Linux system, please see [these instructions](https://github.com/AmpliconSuite/AmpliconSuite-pipeline/blob/master/documentation/CUSTOM_INSTALL.md#getting-mscorefonts-onto-your-system) for making it available to Matplotlib. 



### Option D: Singularity & Docker images 
Containerized versions of AmpliconSuite-pipeline are available for Singularity and Docker.

1. Obtain the AmpliconSuite-pipeline image from the options below:
   - **Singularity**:
     * Singularity installation: https://docs.sylabs.io/guides/3.0/user-guide/installation.html
     * Must have Singularity version 3.6 or later.
     * Pull the singularity image: `singularity pull library://jluebeck/ampliconsuite-pipeline/ampliconsuite-pipeline`

   - **Docker**:
     * Docker installation: https://docs.docker.com/install/
     * Pull the docker image: `docker pull jluebeck/ampliconsuite-pipeline`
    
     * (Optional): Add user to the docker group:
         `sudo usermod -a -G docker $USER` (log out and back in after performing).
   
2. Obtain the execution script and configure the data repo location
    ```bash
    git clone https://github.com/AmpliconSuite/AmpliconSuite-pipeline
    cd AmpliconSuite-pipeline
   # Can use ./install.sh -h to see help before installing
    ./install.sh --finalize_only
    source ~/.bashrc  # either do this or open a new session to make changes live

    ```

3. License for Mosek dependency:
    * [Obtain Mosek license file](https://www.mosek.com/products/academic-licenses/) `mosek.lic`. The license is free for academic use.
    * Place the file in `$HOME/mosek/` (the `mosek/` folder that now exists in your home directory).
    * If you are not able to place the license in the default location, you can set a custom location by exporting the bash variable `MOSEKLM_LICENSE_FILE=/custom/path/`.

   
4. (Recommended) Pre-download AA data repositories and set environment variable AA_DATA_REPO:
   - See the instructions in the section below on obtaining required reference annotations. 
   - If you do not do this process, the container runscript will attempt to download the files into the container before running. This can add to your compute time, especially if you are running many samples.

#### Launching the execution script for the container:

These scripts use most of the same arguments are the main driver script `AmpliconSuite-pipeline.py`
   - Singularity: `AmpliconSuite-pipeline/singularity/run_paa_singularity.py`
     * The Singularity option also exposes `GroupedAnalysisAmpSuite.py` via `AmpliconSuite-pipeline/singularity/run_ga_singularity.py`
   - Docker: `AmpliconSuite-pipeline/docker/run_paa_docker.py`.
     * You can opt to run the docker image as your current user (instead of root) by setting `--run_as_user`. 

An example command might look like:

`AmpliconSuite-pipeline/singularity/run_paa_singularity.py --sif /path/to/ampliconsuite-pipeline.sif -o /path/to/output_dir -s name_of_run -t 8 --bam bamfile.bam --run_AA --run_AC`

### Option E: Custom installation without installer script
Try this if you are going to use `python2`. Please see [the instructions here](documentation/CUSTOM_INSTALL.md).

## Downloading required reference annotations (AA data repo)
Before running AmpliconSuite-pipeline, populate the data repo with required annotations for the reference genomes of interest.
This can be done using the example command below. Specifying `[ref]_indexed` will download a version that includes the BWA index, which is useful for alignment.

`AmpliconSuite-pipeline.py --download_repo [GRCh38|hg19|mm10|... or GRCh38_indexed|hg19_indexed...]`

- Data repo files can also be downloaded manually and placed in the `$AA_DATA_REPO` directory.
   - Go [here](https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/) to see available data repos and copy the URL of the one you want, then
      ```bash
      cd $AA_DATA_REPO
      wget [url of reference_build]
      tar zxf [reference_build].tar.gz
      rm [reference_build].tar.gz
      ```

## Running AmpliconSuite-pipeline
The main driver script for the standalone pipeline is called `AmpliconSuite-pipeline.py`. 

#### Example: Starting from .fastq files, using CNVkit for seed generation.

>`AmpliconSuite-pipeline.py -s sample_name  -t n_threads --fastqs sample_r1.fq.gz sample_r2.fq.gz --ref GRCh38 [--run_AA] [--run_AC]`


* `--run_AA` will invoke AmpliconArchitect directly at the end of the data preparation.
* `--run_AC` will invoke AmpliconClassifier on the AmpliconArchitect outputs.
* `--cnvkit_dir` is only needed if cnvkit.py is not on the system path (typically if it was a custom install).
* CNVkit requires R version 3.5 or greater. This is not standard on older Linux systems. Specify `--rscript_path /path/to/Rscript` with your locally installed current R version if needed. 



#### Example: Starting from .bam, using CNVkit for seed generation

>`AmpliconSuite-pipeline.py -s sample_name -t n_threads --bam sample.bam [--run_AA] [--run_AC]`


#### Example: Starting from .bam and your own whole-genome CNV calls, or an existing AA_CNV_SEEDS.bed
* If using your own CNV calls:


>`AmpliconSuite-pipeline.py -s sample_name -t 1 --cnv_bed your_cnvs.bed --bam sample.bam [--run_AA] [--run_AC]`

Where the CNV bed file reports the following four fields:

`chr    start        end       copy_number`

Additional fields between `end` and `copy_number` may exist, but `copy_number` must always be the last column.

* You can also use the CNVkit `sample_name.cns` file instead of .bed for this argument.


#### Example: Starting from a completed collection of AA output files (reclassification)

This mode allows reclassification of files or uploading of previously completed runs. `--completed_AA_runs` can be a directory or a .tar.gz/.zip file. 
Note that when this mode is used, all AA results must have been generated with respect to the same reference genome version.

>`AmpliconSuite-pipeline.py -s your_collection_name -o output_location -t 1 --completed_AA_runs directory_of_runs/ --ref GRCh38 [--run_AA]


#### Example: Analyzing a collection of related samples (replicates or multi-region sampling)

If you have multiple samples from the same patient, cell line, etc., these should be run as a group to ensure that the same regions are studied across those samples.
Please see the `GroupedAnalysisAmpSuite.py` [example below](#--grouped-analysis-of-related-samples-groupedanalysisampsuitepy) for instructions.


#### Example: Analyzing an oncoviral sample for human-viral hybrid ecDNA detection
Note that users must start with fastq files and `--ref GRCh38_viral` or a bam file aligned to the `AA_DATA_REPO/GRCh38_viral` reference.

>`AmpliconSuite-pipeline.py -s sample_name  -t n_threads --fastqs sample_r1.fq.gz sample_r2.fq.gz --ref GRCh38_viral --cnsize_min 10000 [--run_AA] [--run_AC]`


## Command line arguments to AmpliconSuite-pipeline

#### Downloading data repo files:
- `--download_repo {ref1, ref2, ...}`: This will populate your `$AA_DATA_REPO` directory with files for your reference genome of choice. 
  - Options for this argument are
  `[hg19, GRCh37, GRCh38, mm10, GRCh38_viral, hg19_indexed, GRCh37_indexed, GRCh38_indexed, mm10_indexed, GRCh38_viral_indexed]`. 
  - `*_indexed` refs include the BWA index, only useful if starting from fastqs.

Otherwise, you will instead need these arguments below:

#### Required

- `-o | --output_directory {outdir}`: Directory where results will be stored. Defaults to current directory.

- `-s | --sample_name {sname}`: A name for the sample being run.

- `-t | --nthreads {int}`: Number of threads to use for BWA and CNVkit. Recommend 12 or more threads to be used.

- One of the following input files:

  * `--bam {sample.cs.bam}` Coordinate-sorted bam
  * `--fastqs {sample_r1.fq.gz sample_r2.fq.gz}` Two fastqs (r1 & r2)
  * `--completed_AA_runs {/path/to/some/AA_outputs}`, a directory with AA output files (one or more samples).

#### To invoke downstream analysis

- `--run_AA`: Run AA at the end of the preparation pipeline.

- `--run_AC`: Run AmpliconClassifier following AA. No effect if `--run_AA` not set.

#### Optional

- `--ref {ref name} `: Name of ref genome version, one of `"hg19","GRCh37","GRCh38","GRCh38_viral","mm10","GRCm38"`. This will be auto-detected if it is not set. Required with fastq input.

- `--cnv_bed {cnvfile.bed}` Supply your own CNV calls. Bed file with CN estimate in last column, or the CNVkit `sample.cns` file. If not specified, CNVkit will be called by the wrapper.

- `--cnvkit_dir {/path/to/cnvkit.py}` Path to directory containing `cnvkit.py`. Required if CNVkit was installed from source and `--cnv_bed [cnvfile.bed]` is not given.

- `--normal_bam {matched_normal.bam}` Specify a matched normal BAM file for CNVkit. Not used by AA itself.

- `--completed_run_metadata {run_metadata.json}`, Required only if starting with completed results (`--completed_AA_runs`). Specify a run metadata file for previously generated AA results. If you do not have it, set to 'None'." 

- `--rscript_path {/path/to/Rscript}` (Relevant if using CNVkit and system Rscript version is < 3.5). Specify a path to a local installation of Rscript.

- `--python3_path {/path/to/python3}` Specify custom path to python3 if needed when using CNVkit.

- `--aa_python_interpreter {/path/to/python}` By default PrepareAA will use the system's default `python` path. If you would like to use a different python version with AA, set this to either the path to the interpreter or `python3` or `python2` (default `python`)

[//]: # (- `--freebayes_dir` &#40;Currently deprecated&#41; Specify custom path to freebayes installation folder &#40;not path to executable&#41;. Assumes freebayes on system path if not set. Please note this flag is currently deprecated.)

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

- `--AA_insert_sdevs {float}`: Default 3.0. Suggest raising to 8 or 9 if library has poorly controlled insert size (low fraction of properly-paired reads). See AA documentation for more info.

- `--pair_support_min {int}`, Default is auto-detected by AA based on downsampling parameter, but 2 for default downsampling. This is the minimum number of reads required for breakpoint support.

- `--foldback_pair_support_min {int}`, Number of read pairs for minimum foldback SV support. Default is the same value as `pair_support_min`, however value will be the maximum of `pair_support_min` and this argument. Raising to 3 will help dramatically in heavily artifacted samples (e.g. FFPE).

- `--samtools_path`: Path to a specific samtools binary for use (e.g., /path/to/my/samtools). Uses samtools on system path by default.

- `--sv_vcf`: Provide a VCF file of externally called SVs to augment SVs identified by AA internally.

- `--sv_vcf_no_filter`: Use all external SV calls from the --sv_vcf arg, even those without 'PASS' in the FILTER column.

#### Optional AmpliconRepository upload arguments
- `--upload`: Sets condition to upload after run completes.

- `--project_uuid {str}`: Project UUID where sample(s) will be uploaded to (required if `--upload` set). Obtained from key icon on AmpRepo project page.

- `--project_key {str}`: Project secret key (required if `--upload` set). Obtained from key icon on AmpRepo project page.

- `--username {str}`: Username on AmpliconRepository

- `--upload_server local|dev|prod`: This is a developer argument to assist with debugging upload functionality. The default is `prod` which is AmpliconRepository.org.


## Interpreting classification outputs
- Information about the amplicon classification files produced at the end of the workflow are available [here](https://github.com/AmpliconSuite/AmpliconClassifier?tab=readme-ov-file#3-outputs).
- Information on interpreting the AA cycles file is available [here](https://github.com/jluebeck/AmpliconArchitect#interpreting-the-aa-cycles-files).

## Testing your installation
To ensure your local installation of, we provide a small test dataset (~3Gb), which users can download from SRA. After 
obtaining either the BAM or FASTQ files for GBM39_FF-8 provided [on SRA](https://trace.ncbi.nlm.nih.gov/Traces/?run=SRR4009231), users 
can run the method and compare their results against the files in the `test_outputs/` directory.

## FAQ
Please check out our [guide document](documentation/GUIDE.md).

## Upload Results to AmpliconRepository.org

AmpliconSuite-pipeline can automatically upload your analysis results to [AmpliconRepository.org](https://ampliconrepository.org), a public repository for sharing amplicon data and results.

### Setup Requirements

1. **Create an account** at [ampliconrepository.org](https://ampliconrepository.org)

2. **Create or access a project**:
   - Create a new project, or 
   - Enter an existing project where you are listed as an owner
   - Project can be either public or private

3. **Get your project credentials**:
   - Click the key icon next to your project name
   - Copy the **Project UUID** and **Project Key** that are displayed

### Upload During Pipeline Execution

To upload results automatically when running the full pipeline, add these arguments:

```bash
AmpliconSuite-pipeline.py \
  --run_AA --run_AC \
  --upload \
  --project_uuid "your-project-uuid-here" \
  --project_key "your-project-key-here" \
  --username "your-username@domain.com" \
  [other pipeline arguments...]
```

**Required arguments for upload:**
- `--upload` - Enable upload functionality
- `--project_uuid` - Your project's UUID from the key icon
- `--project_key` - Your project's key from the key icon  
- `--username` - Your AmpliconRepository account email
- `--run_AA` and `--run_AC` - Both required for upload (ensures complete results)

### Upload Completed Results

You can also upload results from previously completed AmpliconArchitect runs:

```bash
AmpliconSuite-pipeline.py \
  --completed_AA_runs /path/to/AA_results/ \
  --ref GRCh38 \
  --sample_name sample_name \
  --upload \
  --project_uuid "your-project-uuid-here" \
  --project_key "your-project-key-here" \
  --username "your-email@domain.com"
```

This mode will run AmpliconClassifier on the existing results and upload everything to your project.

**Note:** Upload requires both AmpliconArchitect and AmpliconClassifier results to ensure complete analysis data is shared.

The upload function identifies the AA-associated files from the directory and zips them for the site, leaving behind sequencing data files or other archives that may be mixed in.

## Citing
If using AmpliconSuite-pipeline in your publication, please cite all the relevant modules used in the analysis, which are summarized in [CITATIONS.md](https://github.com/jluebeck/AmpliconSuite-pipeline/blob/master/CITATIONS.md).

## Additional analysis tools and scripts

### - Grouped analysis of related samples `GroupedAnalysisAmpSuite.py`
For samples derived from a common origin (longitudinal, multiregional sampling from the same source material), 
it is advised that the seed intervals be unified before running AA in order to provide the best comparability between runs. 
We provide a script `GroupedAnalysisAmpSuite.py` which automates this analysis. 
`GroupedAnalysisAmpSuite.py` takes almost all the same arguments as `AmpliconSuite-pipeline.py`, 
however it requires an additional input file, listing the inputs. This file is to be formatted as follows (columns separated by whitespace)

`sample_name` `bamfile` `"tumor"|"normal"` `[CNV_calls.bed]` `[sample_metadata.json]` `[SV_calls.vcf]`

Where `CNV_calls.bed`, `sample_metadata.json`, `SV_calls.vcf` are all optional. All samples listed in each file should be uniquely named and from the same group of related samples. These arguments are positional, so if `CNV_calls` is skipped, it should be set as either `NA` or `None`.
Do not include different collections of related samples in the same table - make different tables. 

AA and AC will be run by default, but can be disabled with `--no_AA`.

#### Native installation example:
> `GroupedAnalysisAmpSuite.py -i {inputs.txt} -o {output_dir} -t {num_threads} --ref {ref_genome}`

#### Singularity example:
(See installation Option D)
> `AmpliconSuite-pipeline/singularity/run_ga_singularity.py --sif {path_to_sif} -i {inputs.txt} -o {output_dir} -t {num_threads} --ref {ref_genome}`

Note: The Singularity wrapper (`run_ga_singularity.py`) automatically handles file path mounting and uses a hybrid approach where large BAM files are directly mounted while smaller files (CNV calls, metadata, SV calls) are copied to a shared location to minimize mount points.

### - **C**andidate **AM**plicon **P**ath **E**numerato**R** `CAMPER.py`
Exahustively search an AA graph file for longest paths (cyclic and non-cyclic). A median amplicon copy number must be specified, or the script will attempt to estimate on its own.
`CAMPER.py` rescales the copy numbers by the median to estimate the multiplicity of each segment within the amplicon, and then 
searches for plausible longest paths explaining the copy number multiplicities. This is useful for identifiying some candidate ecDNA structures.
The output will be an AA-formatted cycles file with additional annotations for length and quality control filter status.
The quality filters take into account root-mean-square residual of copy numbers ("RMSR", lower score is better), as well as "DBI" representing the Davies-Bouldin index of copy-number to multiplicity clustering. More information on the method can be found in the [methods section of this publication](https://www.nature.com/articles/s41588-022-01190-0).
The first entry (Cycle1) will be a cyclic path, while the second entry (Cycle2) will be a non-cyclic path. A full explanation of arguments is available with `-h`. Note that this should only be applied to AA amplicons with at most 1 ecDNA present in the AA amplicon (multiple-species reconstruction not supported).

> `AmpliconSuite-pipeline/scripts/CAMPER.py -g sample_amplicon1_graph.txt [--scaling_factor (CN estimate value)] [--remove_short_jumps] [--keep_all_LC] [--max_length (value in kbp)]`

### - `breakpoints_to_bed.py`
Requires `intervaltree` python package. Write discordant edges (breakpoint junctions) from an AA graph into a pseudo-bed file.
The `.input` file is automatically produced by AC, but is formatted like so

`samplename /path/to/sample_amplicon1_cycles.txt /path/to/sample_amplicon1_graph.txt`

Usage:
>`scripts/breakpoints_to_bed.py -i (AC.input) [--regions chrA:start-stop chrB:start-stop ...]`

### - `convert_cns_to_bed.py`
Many users will choose to run CNVkit outside AmpliconSuite-pipeline and then want to use the CNVkit calls in AA. We recommend using the `.cns` file as a source for the seeds. 
Note the `.call.cns` file is different and contains more aggressively merged CNV calls, which we do not recommend as a source of seeds. As the `.cns` file specifies a log2 ratio,
we provide the following script to reformat the `.cns` file from CNVkit into a `.bed` file. 

This script should not be needed by most users, as the `*_cnvkit_output/` directory will already contain a 
`.bed` of genome-wide CNV calls produced by CNVkit. 

Usage:
>`scripts/convert_cns_to_bed.py sample.cns`


### - `cycles_to_bed.py`
Requires `intervaltree` python package. Write an AA cycles file as a series of bed files, one for each decomposition. Writes two types of bed files, an `unordered_cycle` file, where segments are merged and sorted, and order and orientation of segments is lost, and also writes
an ordered file where the order and orientation of the genome segments comprising the cycle is maintained.

Usage:

>`scripts/cycles_to_bed.py -c sample_amplicon1_cycles.txt`

### - `graph_cleaner.py`
Requires `intervaltree` python package. Poorly controlled insert size can lead to numerous spurious short breakpoint edges. This enables users to remove SV edges from the AA graph file based on type of SV, support and SV size. 
Namely, for artifact removal it can clean the very short everted (inside-out read pair) duplication-like orientation edges characteristic of sequencing artifact. These will appear as numerous short brown 'spikes' in the AA amplicon image. 
Most artifacts should be filtered though by increasing AA's `--AA_insert_sdevs` paramter to 9 or higher. However users can still use this script for other cleaning purposes (e.g. SVs with weak support).


Usage:

>`scripts/graph_cleaner.py -g /path/to/sample_ampliconx_graph.txt [--max_hop_size 4000] `

or

>`scripts/graph_cleaner.py --graph_list /path/to/list_of_graphfiles.txt [--max_hop_size 4000] `


This will output an AA graph file(s) `/path/to/my_sample_ampliconX_cleaned_graph.txt`.

### - `graph_to_bed.py`
Requires `intervaltree` python package. Create a bed file of the graph segments and a bedpe file of the disordant graph edges. Can also filter to only get segments with CN above `--min_cn`. 
Setting `--unmerged` will not merge adjacent graph segments and will print the graph segment CN in the last column.

Usage:

>`scripts/graph_to_bed.py -g sample_amplicon_graph.txt [--unmerged] [--min_cn 0] [--add_chr_tag]`

