# PrepareAA
![GitHub](https://img.shields.io/github/license/jluebeck/PrepareAA)

A multithread-enabled quickstart tool for [AmpliconArchitect](https://github.com/jluebeck/AmpliconArchitect). 
Performs all preliminary steps (alignment, CNV calling, seed interval detection) required prior to running AmpliconArchitect. 
PrepareAA supports hg19, GRCh37, GRCh38 (hg38) and mouse genome mm10 (GRCm38). PrepareAA can invoked to begin at any intermediate stage of the data preparation process and can invoke both AmpliconArchitect and AmpliconClassifier.
**Current version: 0.1203.9**

Please check out the [**detailed guide**](https://github.com/jluebeck/PrepareAA/blob/master/GUIDE.md) on running AA to learn about best practices and see some FAQs.


## Prerequisites:
PrepareAA supports both `python2` and `python3`, however AmpliconArchitect currently requires `python2` and CNVKit requires `python3`.
`Python3` support for AmpliconArchitect is coming soon. 

Depending on what input data you are starting from, PrepareAA (PAA) may require the following tools to be installed beforehand:
- (required) The [jluebeck/AmpliconArchictect fork](https://github.com/jluebeck/AmpliconArchitect) must be installed.
- (required) The latest AmpliconArchitect [data repo](https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/).
  - mm10 and larger versions of the individual data repos containing bwa index files are also provided [here](https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/). Indexed version recommended if starting from unaligned fastq reads.
- (optional) [AmpliconClassifier](https://github.com/jluebeck/AmpliconClassifier) to generate classifications of AmpliconArchitect outputs.
- (optional) [bwa mem](https://github.com/lh3/bwa) (unless supplying your own BAM file)
- (optional) [samtools](http://www.htslib.org/) (unless you already have a coordinate-sorted and indexed BAM file).
- (optional) [CNVkit](https://github.com/etal/cnvkit) or [Canvas](https://github.com/Illumina/canvas) or  (unless supplying your own CNV calls).
  - (required for Canvas) [freebayes](https://github.com/ekg/freebayes) version 1.3.1 or greater, (unless providing your own VCF calls to Canvas)
- Some optional scripts packaged with PrepareAA require the `numpy`, `matplotlib` and `intervaltree` python packages. Can be installed with `pip`, `conda` or similar. 

PrepareAA assumes both samtools and bwa executables are on the system path and can be directly invoked from bash without pathing to the executables.

PrepareAA has been tested with Ubuntu (16.04 and above) and CentOS 7. PrepareAA's optional dependencies related to CNV calling will not work on CentOS 6.


**Note on using CNVKit**: We currently recommend using CNVKit for identification of AA seeds. Please note that CNVKit requires
`python3`. It also requires `R` version >= 3.5, which is non-standard on Ubuntu 16.04/14.04.

**Note on using Canvas**: If using Canvas, please make sure the Canvas reference genome files are located in the expected location for Canvas. To do this, you can follow instructions on the Canvas Github page. We also provide a script `$ install_canvas.sh [path/to/installation/directory/`,
which when run from the PrepareAA source directory will fetch the Canvas binary and download the `canvasdata` data repository. If installing on your own, create the canvasdata/ reference genome sudirectories in the folder with the Canvas executable. One installation dependency not mentioned explictly on the Canvas Readme is `dotnet-sdk-2.2`, which can be obtained in Ubuntu by running `sudo apt-get install dotnet-sdk-2.2`. 


### Standalone installation

In the directory you want to run AA in, do 

`git clone https://github.com/jluebeck/PrepareAA.git`

Please see the [jluebeck/AmpliconArchitect fork](https://github.com/jluebeck/AmpliconArchitect) for AA installation instructions. AA must be installed to use PAA.

Prepare AA will generate a BWA index for the reference genome if one is not yet in place. This adds >1hr to running time for the first use only when alignment is performed. Data repos with BWA index pre-generated are available [here](https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/).

PrepareAA with CNVKit will also function on coordinate-sorted CRAM files, [provided that the CRAM reference is in place](http://www.htslib.org/workflow/#:~:text=One%20of%20the%20key%20concepts,genome%20used%20to%20generate%20it.).


### PrepareAA Docker 
A dockerized version of PAA is [available on dockerhub](https://hub.docker.com/repository/docker/jluebeck/prepareaa) or can be built using the Dockerfile in the `docker/` folder. It will install bwa, CNVKit and AmpliconArchitect inside the docker image. Running this docker image can be done as follows:

1. Docker:
    * Install docker: `https://docs.docker.com/install/`
    * (Optional): Add user to the docker group and relogin:
        `sudo usermod -a -G docker $USER`
2. License for Mosek optimization tool:
    * Obtain license file `mosek.lic` (`https://www.mosek.com/products/academic-licenses/` or `https://www.mosek.com/try/`)
    * `export MOSEKLM_LICENSE_FILE=<Parent directory of mosek.lic> >> ~/.bashrc && source ~/.bashrc`
3. Download AA data repositories and set environment variable AA_DATA_REPO:
    * Download [here](https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/) to download data repos with (`_indexed`) or
    without the bwa reference index included.
    * Set enviroment variable AA_DATA_REPO to point to the data_repo directory:
        ```bash
        tar zxf data_repo.tar.gz
        echo export AA_DATA_REPO=$PWD/data_repo >> ~/.bashrc
        cd $AA_DATA_REPO && touch coverage.stats && chmod a+rw coverage.stats
        source ~/.bashrc
        ```
#### Obtain PrepareAA image and execution script:
1. Clone GitHub repository to access the runscript
    * `git clone https://github.com/jluebeck/PrepareAA.git`

2. Run the script `run_paa_docker.py` located in `PrepareAA/docker`. It uses (most of) the same command line arguments one would pass to `PrepareAA.py`. CNV calling with CNVKit is integrated into the docker image (with help from Owen Chapman).

An example docker command might look like:

`PrepareAA/docker/run_paa_docker.py -o /path/to/output_dir -s name_of_run -t 8 --bam /path/to/bamfile.bam --run_AA --run_AC`

**You can opt to run the docker image as your current user by setting `--run_as_user`.** 

### Other usage options
PAA can also be run through Nextflow, using the [nf-core/circdna pipeline](https://nf-co.re/circdna) constructed by [Daniel Schreyer](https://github.com/DSchreyer).

## Usage
Two example standard runs of PrepareAA:

#### Starting from .fastq files, using Canvas for seed generation.
```
/path/to/PrepareAA/PrepareAA.py -s sample_name  -t number_of_threads --canvas_dir /path/to/canvas/canvas_exec_dir --fastqs sample_r1.fastq.gz sample_r2.fastq.gz --ref hg19 [--run_AA] [--run_AC]
```

or

#### Starting from sorted .bam, using CNVkit for seed generation
```
/path/to/PrepareAA/PrepareAA.py -s sample_name -t number_of_threads --cnvkit_dir /path/to/cnvkit.py --bam sample.cs.rmdup.bam [--run_AA] [--run_AC]
```

`--run_AA` will invoke AmpliconArchitect directly at the end of the data preparation.
`--run_AC` will invoke AmpliconClassifier on the AmpliconArchitect outputs.


##### Starting from BAM and your own CNV calls (or recycled AA_CNV_SEEDS.bed)
* If you already have your coordinate-sorted bam file, `--fastqs` can be replaced with `--bam`.


* If using your own CNV calls:
```
/path/to/PrepareAA/PrepareAA.py -s sample_name -t number_of_threads --cnv_bed your_cnvs.bed (--fastqs sample_r1.fastq sample_r2.fastq | --bam sample.cs.bam) [--run_AA] [--run_AC]
```
Where the CNV bed file is formatted as (**without a header present**):

`chr    start        end       copy_number`

Additional fields between `end` and `copy_number` may exist, but `copy_number` must always be the last column.

* You can also use a CNVKit .cns file instead of .bed for this argument.

* CNVkit requires R version 3.5 or greater. This is not standard on many Linux systems. Specify `--rscript_path /path/to/Rscript` with your locally installed current R version if needed. 

* If you generated your own VCF but would still like to use Canvas CNV, you can supply `--vcf` to bypass the freebayes step.

* If using your own VCF + Canvas: Canvas only considers sites with "PASS" in the FILTER field of the VCF, so if "." is used, Canvas will fail. If you would like to convert your VCF with "." in the FILTER field to "PASS", you can use the following awk command
```
cat your_file.vcf | "awk '{ if (substr($1,1,1) != \"#\" ) { $7 = ($7 == \".\" ? \"PASS\" : $7 ) }} 1 ' OFS=\"\\t\"" > your_reformatted_file.vcf
```

#### Starting from completed AA results
If the user has one or more AA results directories inside a directory, the user can use PrepareAA to call AmpliconClassifier with default settings.
```
/path/to/PrepareAA/PrepareAA.py -s project_name --completed_AA_runs /path/to/location_of_all_AA_results/ --completed_run_metadata [representative_run_metadata_file].json -t 1 --ref hg38
```

Note that when this mode is used all AA results must have been generated with respect to the same reference genome version.

### Command line arguments to PrepareAA

- `-o | --output_directory [outdir]`: (Optional) Directory where results will be stored. Defaults to current directory.

- `-s | --sample_name [sname]`: (Required) A name for the sample being run.

- `-t | --nthreads [numthreads]`: (Required) Number of threads to use for BWA and freebayes. We do not control thread usage of Canvas. Recommend 12 or more threads to be used.

- `--bam | --sorted_bam [sample.cs.bam]` **OR** `--fastqs [sample_r1.fq[.gz] sample_r2.fq[.gz]]` (Required) Input files. Two fastqs (r1 & r2) or a coordinate sorted bam **OR** .

- `--canvas_dir [/path/to/Canvas_files/]` (Required if not `--reuse_canvas` and not `--cnv_bed [cnvfile.bed]` and not `--cnvkit_dir`) Path to directory containing the Canvas executable and `canvasdata/` subdirectory.

- `--cnvkit_dir [/path/to/cnvkit.py]` (Required if not `--reuse_canvas` and not `--cnv_bed [cnvfile.bed]` and not `--canvas_dir`) Path to directory containing cnvkit.py.

- `--completed_run_metadata`, (Required if startng with completed results). Specify a run metadata file for previously generated AA results. If you do not have it, set to 'None'." 

- `--rscript_path [/path/to/Rscript]` (Required if system Rscript version < 3.5 and using `--cnvkit_dir`). Specify a path to a local installation of Rscript compatible with CNVkit.

- `--python3_path` (Optional) Specify custom path to python3, if needed when using CNVKit (which requires python3).

- `--aa_python_interpreter` (Optional) By default PrepareAA will use the system's default `python` path. If you would like to use a different python version with AA, set this to either the path to the interpreter or `python3` or `python2` (default `python`)

- `--freebayes_dir` (Optional) Specify custom path to freebayes installation folder (not path to executable). Only applied if using Canvas. Assumes freebayes on system path if not set.

- `--run_AA`: (Optional) Run AA at the end of the preparation pipeline.

- `--run_AC`: (Optional) Run AmpliconClassifier following AA. No effect if `--run_AA` not set.

- `--ref `: Name of ref genome version ("hg19","GRCh37","GRCh38","mm10","GRCm38"). This will be auto-detected if it is not set.

- `--vcf [your_file.vcf]`: (Optional) Supply your own VCF to skip the freebayes step.

- `--cngain [float]`: (Optional) Set a custom threshold for the CN gain considered by AA. Default: 4.5.

- `--cnsize_min [int]`: (Optional) Set a custom threshold for CN interval size considered by AA. Default: 50000.

- `--downsample [float]`: (Optional) Set a custom threshold for bam coverage downsampling during AA. Does not affect coverage in analyses outside of AA. Default: 10.

- `--use_old_samtools`: (Optional) Set this flag if your Samtools version is < 1.0. Default: False.

- `--reuse_canvas` (Optional) Reuse the Canvas results from a previous run. Default: False

- `--cnv_bed [cnvfile.bed]` (Optional) Supply your own CNV calls, bypasses freebayes and Canvas steps. Bed file with CN estimate in last column or CNVKit .cns file.

- `--no_filter`: (Optional) Do not invoke `amplified_intervals.py` to filter amplified seed regions based on CN, size and ignorefile regions. 

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
 

### FAQ
Check out the [guide document](https://github.com/jluebeck/PrepareAA/blob/master/GUIDE.md)!

### Citing
If using PrepareAA in your publication, please cite the [AmpliconArchitect article](https://www.nature.com/articles/s41467-018-08200-y). If using PrepareAA to wrap other tools (like [CNVkit](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004873)), please cite those tools as well.


## Additional analysis tools and scripts

### - **C**andidate **AM**plicon **P**ath **E**numerato**R** `CAMPER.py`
Formerly called `CAMPER.py`. Exahustively search an AA graph file for longest paths (cyclic and non-cyclic). A median amplicon copy number must be specified, or the script will attempt to estimate on its own.
`CAMPER.py` rescales the copy numbers by the median to estimate the multiplicity of each segment within the amplicon, and then 
searches for plausible longest paths explaining the copy number multiplicities. This is useful for identifiying some candidate ecDNA structures.
The output will be an AA-formatted cycles file with additional annotations for length and quality control filter status. The quality filters take into account root mean square residual of copy numbers ("RMSR", lower score is better), as well as "DBI" representing the Davies-Bouldin index of copy-number to multiplicity clustering. More information on the method can be found in the [methods section of this pre-print](https://www.biorxiv.org/content/10.1101/2021.11.28.470285v1).
The first entry (Cycle1) will be a cyclic path, while the second entry (Cycle2) will be a non-cyclic path. A full explanation of arguments is available with `-h`.

`PrepareAA/scripts/plausible_paths.py -g sample_amplicon1_graph.txt [--scaling_factor (CN estimate value)] [--remove_short_jumps] [--keep_all_LC] [--max_length (value in kbp)]`

   ### - `breakpoints_to_bed.py`
Requires `intervaltree` python package pre-installed. Write discordant edges (breakpoint junctions) from an AA graph into a pseudo-bed file.

   ### - `convert_cns_to_bed.py`
Many users will choose to run CNVKit outside of PrepareAA and then want to use the CNVKit calls in AA. We recommend using the `.cns` file as a source for the seeds. 
Note the `.call.cns` file is different and contains more aggressively merged CNV calls, which we do not recommend as a source of seeds. As the `.cns` file specifies a log2 ratio,
we provide the following script to reformat the `.cns` file from CNVKit into a `.bed` file useable with PrepareAA. 

Usage:
`./scripts/convert_cns_to_bed.py your_CNVKit_output/your_sample.cns`

This will output a bed file which can be fed into PrepareAA. 

   ### - `cycles_to_bed.py`
Requires `intervaltree` python package pre-installed. Write an AA cycles file as a series of bed files, one for each decomposition. Segments are merged and sorted, and order and orientation of segments is lost.

   ### - `seed_trimmer.py`
AA seeds are not designed to be larger than 10 Mbp - as that passes the upper limit of what is considered a 'focal amplification'.
To pre-filter some of these seeds and break them on regions AA cannot analyze (low mappability, centromeres, segmental duplications), we provide the following script,
which can and should be invoked on any seeds > 10 Mbp. This script should be run prior to running PrepareAA (or `amplified_intervals.py` if not using PrepareAA).

Usage:

`./scripts/seed_trimmer.py --cnv_bed /path/to/my_cnvs.bed --ref hg19/GRCh37/GRCh38/mm10 [--minsize 50000] [--cngain 4.5]`

This will output a bed file `/path/to/my_seeds_trimmed.bed`, which can then be fed into `amplified_intervals.py`. 


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

`./scripts/graph_to_bed.py -g /path/to/sample_amplicon_graph.txt [--unmerged] [--min_cn 0] [--add_chr_tag]`


### - `bfb_foldback_detection.py [deprecated]`
**This script is deprecated and no longer supported, but available for legacy purposes. For more robust BFB detection, please try out [AmpliconClassifier](https://github.com/jluebeck/AmpliconClassifier).**

Requires `intervaltree` python package pre-installed. Script can be used to detect possible BFB-like signatures from AA graph files (documentation below).

To use the `bfb_foldback_detection.py` script on AA output, please create a two column file with the name of the graph file in column 1 and the path to the graph file in column 2. The rest of the command-line arguments are as follows.


##### Required arguments for running on AA results
-  `--exclude [path to $AA_DATA_REPO/[ref]/[mappability excludable file]`

- `-o [output filename prefix]`

- `--ref [hg19, GRCh37, GRCh38]`

- `--AA_graph_list [two-column file listing AA graphs]`

