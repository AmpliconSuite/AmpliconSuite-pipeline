## PrepareAA

A multithread-enabled quickstart tool for [AmpliconArchitect](https://github.com/jluebeck/AmpliconArchitect). 
Performs all preliminary steps (alignment, CNV calling, seed interval detection) required prior to running AmpliconArchitect. 
PrepareAA supports hg19 (or GRCh37) and hg38. PrepareAA can also be invoked to start at intermediate stages of the data preparation process.

Please check out the **detailed guide** on running AA [available here](https://github.com/jluebeck/PrepareAA/blob/master/GUIDE.md) to learn about best practices and see some FAQs.

### Prerequisites:
Depending on what input data you are using, PrepareAA (PAA) may require the following tools to be installed beforehand:
- The [jluebeck/AmpliconArchictect fork](https://github.com/jluebeck/AmpliconArchitect) is recommended for PrepareAA. The development AmpliconArchitect [data repo](https://drive.google.com/drive/folders/18T83A12CfipB0pnGsTs3s-Qqji6sSvCu) must be downloaded and used.
- [bwa mem](https://github.com/lh3/bwa) (unless supplying your own BAM file aligned to the AA reference genome)
- [samtools](http://www.htslib.org/) (unless you already have a coordinate-sorted BAM file. PrepareAA supports versions >= 1.0 and < 1.0)
- [CNVkit](https://github.com/etal/cnvkit) or [Canvas](https://github.com/Illumina/canvas) or  (unless supplying your own CNV calls).
- (If using Canvas) [freebayes](https://github.com/ekg/freebayes) (version 1.3.1 or greater, freebayes is only required if using Canvas - but not if supplying your own VCF calls to Canvas)

**Note on using Canvas**: If using Canvas, please make sure the Canvas reference genome files are located in the expected location for Canvas. To do this, you can follow instructions on the Canvas Github page. We also provide a script `$ install_canvas.sh [path/to/installation/directory/`, which when run from the PrepareAA source directory will fetch the Canvas binary and download the `canvasdata` data repository. If installing on your own, create the canvasdata/ reference genome sudirectories in the folder with the Canvas executable. One installation dependency not mentioned explictly on the Canvas Readme is `dotnet-sdk-2.2`, which can be obtained in Ubuntu by running `sudo apt-get install dotnet-sdk-2.2`. 

**Note on using CNVkit**: Please note that CNVkit requires `R` version >= 3.5, which is non-standard on Ubuntu 16.04/14.04.

PrepareAA assumes both samtools and bwa executables are on the system path and can be directly invoked from bash without pathing to the executables. This is already standard for most users.

### Standalone configuration

In the directory you want to run AA in, do 

`git clone https://github.com/jluebeck/PrepareAA.git`

Please see the [jluebeck/AmpliconArchitect fork]((https://github.com/jluebeck/AmpliconArchitect)) for AA installation instructions. AA must be installed to use PAA.

Prepare AA will generate a BWA index for the reference genome if one is not yet in place. This adds >1hr to running time for the first use only when alignment is performed.

### PrepareAA Docker 
A dockerized version of PAA is available in the docker folder. It will install bwa, CNVKit and AmpliconArchitect inside the docker image. Running this docker image can be done as follows:

1. Docker:
    * Install docker: `https://docs.docker.com/install/`
    * (Optional): Add user to the docker group and relogin:
        `sudo usermod -a -G docker $USER`
2. License for Mosek optimization tool:
    * Obtain license file `mosek.lic` (`https://www.mosek.com/products/academic-licenses/` or `https://www.mosek.com/try/`)
    * `export MOSEKLM_LICENSE_FILE=<Parent directory of mosek.lic> >> ~/.bashrc && source ~/.bashrc`
3. Download AA data repositories and set environment variable AA_DATA_REPO:
    * Download from `https://drive.google.com/drive/folders/18T83A12CfipB0pnGsTs3s-Qqji6sSvCu`
    * Set enviroment variable AA_DATA_REPO to point to the data_repo directory:
        ```bash
        tar zxf data_repo.tar.gz
        echo export AA_DATA_REPO=$PWD/data_repo >> ~/.bashrc
        source ~/.bashrc
        ```
#### Obtain PrepareAA image and execution script:
1. Pull docker image:
    * `docker pull jluebeck/prepareaa`

2. Run the script `run_paa_docker.py` located in `PrepareAA/docker`. It uses (most of) the same command line arguments one would pass to `PrepareAA.py`. CNV calling with CNVKit is integrated into the docker image (with help from Owen Chapman).

An example command might look like:

`PrepareAA/docker/run_paa_docker.py -o /path/to/output_dir -s name_of_run -t 8 --ref hg19 --sorted_bam /path/to/bamfile.bam --run_cnvkit`

**Please make sure your output directory specified for `-o` is writeable.** 

### Usage
Two example standard runs of PrepareAA:

#### Starting from .fastq files, using Canvas for seed generation.
```
/path/to/PrepareAA/PrepareAA.py -s sample_name  -t number_of_threads --canvas_dir /path/to/canvas/canvas_exec_dir --fastqs sample_r1.fastq.gz sample_r2.fastq.gz [--run_AA]
```

or

#### Starting from sorted .bam, using CNVkit for seed generation
```
/path/to/PrepareAA/PrepareAA.py -s sample_name  -t number_of_threads --cnvkit_dir /path/to/cnvkit.py --sorted_bam sample.cs.rmdup.bam [--run_AA]
```

`--run_AA` will invoke AmpliconArchitect directly at the end of the data preparation.

PrepareAA with CNVKit will also function on sorted CRAM files, [provided that the CRAM reference is in place](http://www.htslib.org/workflow/#:~:text=One%20of%20the%20key%20concepts,genome%20used%20to%20generate%20it.).

##### Starting from intermediate steps
* If you already have your coordinate-sorted bam file, `--fastqs` can be replaced with `--sorted_bam`.


* If using your own CNV calls:
```
/path/to/PrepareAA/PrepareAA.py -s sample_name  -t number_of_threads  --cnv_bed your_cnvs.bed  (--fastqs sample_r1.fastq sample_r2.fastq | --sorted_bam sample.cs.bam) [--run_AA]
```
Where the CNV bed file is formatted as (**without a header present**):

`chr    start        end       copy_number`

Additional fields between `end` and `copy_number` may exist, but `copy_number` must always be the last column.

* CNVkit requires R version 3.5 or greater. This is not standard on many Linux systems. Specify `--rscript_path /path/to/Rscript` with your locally installed current R version if needed. 

* If you generated your own VCF but would still like to use Canvas CNV, you can supply `--vcf` to bypass the freebayes step.

* If using your own VCF + Canvas: Canvas only considers sites with "PASS" in the FILTER field of the VCF, so if "." is used, Canvas will fail. If you would like to convert your VCF with "." in the FILTER field to "PASS", you can use the following awk command
```
cat your_file.vcf | "awk '{ if (substr($1,1,1) != \"#\" ) { $7 = ($7 == \".\" ? \"PASS\" : $7 ) }} 1 ' OFS=\"\\t\"" > your_reformatted_file.vcf
```

A description of other command line arguments for PrepareAA is provided below

### Command line arguments to PrepareAA

- `-o | --output_directory [outdir]`: (Optional) Directory where results will be stored. Defaults to current directory.

- `-s | --sample_name [sname]`: (Required) A name for the sample being run.

- `-t | --nthreads [numthreads]`: (Required) Number of threads to use for BWA and freebayes. We do not control thread usage of Canvas. Recommend 12 or more threads to be used.

- `--canvas_dir [/path/to/Canvas_files/]` (Required if not `--reuse_canvas` and not `--cnv_bed [cnvfile.bed]` and not `--cnvkit_dir`) Path to directory containing the Canvas executable and `canvasdata/` subdirectory.

- `--cnvkit_dir [/path/to/cnvkit.py]` (Required if not `--reuse_canvas` and not `--cnv_bed [cnvfile.bed]` and not `--canvas_dir`) Path to directory containing cnvkit.py.

- `--rscript_path [/path/to/Rscript]` (Required if system Rscript version < 3.5 and using `--cnvkit_dir`). Specify a path to a local installation of Rscript compatible with CNVkit.

- `--python3_path` Specify custom path to python3, if needed when using CNVKit (which requires python3)

- `--sorted_bam [sample.cs.bam] | --fastqs [sample_r1.fq[.gz] sample_r2.fq[.gz]]` (Required) Input files. Two fastqs (r1 & r2) or a coordinate sorted bam.

- `--run_AA`: (Optional) Run AA at the end of the preparation pipeline.

- `--ref ["hg19"]`: (Optional) Name of ref genome version ("hg19","GRCh37","GRCh38"). Default: "hg19".

- `--vcf [your_file.vcf]`: (Optional) Supply your own VCF to skip the freebayes step.

- `--cngain [float]`: (Optional) Set a custom threshold for the CN gain considered by AA. Default: 4.5.

- `--cnsize_min [int]`: (Optional) Set a custom threshold for CN interval size considered by AA. Default: 50000.

- `--downsample [float]`: (Optional) Set a custom threshold for bam coverage downsampling during AA. Does not affect coverage in analyses outside of AA. Default: 10.

- `--use_old_samtools`: (Optional) Set this flag if your Samtools version is < 1.0. Default: False.

- `--reuse_canvas` (Optional) Reuse the Canvas results from a previous run. Default: False

- `--cnv_bed [cnvfile.bed]` (Optional) Supply your own CNV calls, bypasses freebayes and Canvas steps.


PrepareAA has been tested with Ubuntu 16.04 and CentOS 7. PrepareAA's dependencies (related to CNV calling) will not work on CentOS 6.


### Citing
If using PrepareAA in your publication, please cite the [AmpliconArchitect article](https://www.nature.com/articles/s41467-018-08200-y). If using PrepareAA to wrap other tools (like [CNVkit](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004873)), please cite those tools as well.


### Additional analysis scripts

   ### - `convert_cns_to_bed.py`
Many users will choose to run CNVKit outside of PrepareAA and then want to use the CNVKit calls in AA. We recommend using the `.cns` file as a source for the seeds. 
Note the `.call.cns` file is different and contains more aggressively merged CNV calls, which we do not recommend as a source of seeds. As the `.cns` file specifies a log2 ratio,
we provide the following script to reformat the `.cns` file from CNVKit into a `.bed` file useable with PrepareAA. 

Usage:
`./scripts/convert_cns_to_bed.py your_CNVKit_output/your_sample.cns`

This will output a bed file which can be fed into PrepareAA. 

   ### - `seed_trimmer.py`
AA seeds are not designed to be larger than 10 Mbp - as that passes the upper limit of what is considered a 'focal amplification'.
To pre-filter some of these seeds and break them on regions AA cannot analyze (low mappability, centromeres, segmental duplications), we provide the following script,
which can and should be invoked on any seeds > 10 Mbp. This script should be run prior to running PrepareAA (or `amplified_intervals.py` if not using PrepareAA).

Usage:

`./scripts/seed_trimmer.py --seeds [/path/to/my_seeds.bed] --ref hg19/GRCh37/GRCh38 [--minsize 50000]`

This will output a bed file `/path/to/my_seeds_trimmed.bed`, which can then be fed into `amplified_intervals.py`. 


### - `graph_cleaner.py`
Sequencing artifacts can lead to numerous spurious short breakpoint edges. This script attempts to remove edges which conform to artifactual profiles. 
Namely, very short everted (inside-out read pair) orientation edges. These will appear as numerous short brown 'spikes' in the AA amplicon image.
This script removes them from the graph file.

Usage:

`./scripts/graph_cleaner.py -g /path/to/sample_ampliconx_graph.txt [--max_hop_size 5000] `

or

`./scripts/graph_cleaner.py --graph_list /path/to/list_of_graphfiles.txt [--max_hop_size 5000] `


This will output an AA graph file(s) `/path/to/my_sample_ampliconX_cleaned_graph.txt`.


   ### - `bfb_foldback_detection.py`
**This script is deprecated, but available for legacy purposes. For more robust BFB detection, please try out [AmpliconClassifier](https://github.com/jluebeck/AmpliconClassifier).**


Script can be used to detect possible BFB-like signatures from AA graph files (documentation below).

To use the `bfb_foldback_detection.py` script on AA output, please create a two column file with the name of the graph file in column 1 and the path to the graph file in column 2. The rest of the command-line arguments are as follows.


##### Required arguments for running on AA results
-  `--exclude [path to $AA_DATA_REPO/[ref]/[mappability excludable file]`

- `-o [output filename prefix]`

- `--ref [hg19, GRCh37, GRCh38]`

- `--AA_graph_list [two-column file listing AA graphs]`

