## PrepareAA

A multithreaded quickstart tool for [AmpliconArchitect](https://github.com/virajbdeshpande/AmpliconArchitect). Performs all preliminary steps (alignment, CNV calling, seed interval detection) required prior to running AmpliconArchitect. PrepareAA supports hg19, GRCh37, and hg38. PrepareAA can also be invoked to start at intermediate stages of the data preparation process.

### Prerequisites:
PrepareAA (PAA) Requires the following tools to be installed beforehand:
- [AmpliconArchitect](https://github.com/virajbdeshpande/AmpliconArchitect)
- [bwa mem](https://github.com/lh3/bwa) (unless supplying your own BAM file aligned to the AA reference genome)
- [samtools](http://www.htslib.org/) (PrepareAA supports both versions >= 1.0 and < 1.0)
- [Canvas](https://github.com/Illumina/canvas) or [CNVkit](https://github.com/etal/cnvkit) (unless supplying your own CNV calls)
- [freebayes](https://github.com/ekg/freebayes) (version 1.3.1 or greater, freebayes is required if using Canvas - unless supplying your own VCF calls)


If using Canvas please make sure the Canvas reference genome files are located in the expected location for Canvas. To do this, you can follow instructions on the Canvas Github page. We also provide a script `$install_canvas.sh [path/to/installation/directory/`, which when run from the PrepareAA source directory will fetch the Canvas binary and download the `canvasdata` data repository. If installing on your own, create the canvasdata/ reference genome sudirectories in the folder with the Canvas executable. One installation dependency not mentioned explictly on the Canvas Readme is `dotnet-sdk-2.2`, which can be obtained in Ubuntu by running `sudo apt-get install dotnet-sdk-2.2`. 

Please note that CNVkit requires `R` version >= 3.5, which is non-standard on Ubuntu 16.04/14.04.

### Installation
Files in PrepareAA/hg19/ directory must be placed into $AA_DATA_REPO/hg19/ prior to using. If using hg38, you will need to download the hg38 AA data repo patch if you have not done so already. If using hg19 or GRCh37, replacing the files `conserved.bed` or `conserved.gain5.bed` in the data repo with the version included here is recommended for compatability with both standard and non-standard reference versions.

Prepare AA will generate a BWA index for the reference genome if one is not yet in place. This adds >1hr to running time for the first use only.

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

##### Starting from intermediate steps
* If you already have your coordinate-sorted bam file, `--fastqs` can be replaced with `--sorted_bam`.


* If using your own CNV calls:
```
/path/to/PrepareAA/PrepareAA.py -s sample_name  -t number_of_threads  --cnv_bed your_cnvs.bed  (--fastqs sample_r1.fastq sample_r2.fastq | --sorted_bam sample.cs.bam) [--run_AA]
```
Where the CNV bed file is formatted as:

`chrN    start        end     some_arbitrary_name      copy_number`


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

- `--ref ["hg19"]`: (Optional) Name of ref genome version ("hg19","GRCh37","hg38"). Only hg19 currently supported. Default: "hg19".

- `--vcf [your_file.vcf]`: (Optional) Supply your own VCF to skip the freebayes step.

- `--cngain [float]`: (Optional) Set a custom threshold for the CN gain considered by AA. Default: 5.

- `--cnsize_min [int]`: (Optional) Set a custom threshold for CN interval size considered by AA. Default: 50000.

- `--downsample [float]`: (Optional) Set a custom threshold for bam coverage downsampling during AA. Does not affect coverage in analyses outside of AA. Default: 5.

- `--use_old_samtools`: (Optional) Set this flag if your Samtools version is < 1.0. Default: False.

- `--reuse_canvas` (Optional) Reuse the Canvas results from a previous run. Default: False

- `--cnv_bed [cnvfile.bed]` (Optional) Supply your own CNV calls, bypasses freebayes and Canvas steps.


PrepareAA has been tested with Ubuntu 16.04 and CentOS 7. PrepareAA's dependencies (related to CNV calling) will not work on CentOS 6.


### Additional analysis scripts
#### BFB detection with `bfb_foldback_detection.py`
This script can be used to detect possible BFB-like signatures from AA graph files (documentation below). A second mode of operation can be invoked using other custom inputs (documentation coming soon).

To use the `bfb_foldback_detection.py` script on AA output, please create a two column file with the name of the graph file in column 1 and the path to the graph file in column 2. The rest of the command-line arguments are as follows.


##### Required arguments for running on AA results
-  `--exclude [path to $AA_DATA_REPO/[ref]/[mappability excludable file]`

- `-o [output filename prefix]`

- `--ref [hg19, GRCh37, GRCh38]`

- `--AA_graph_list [two-column file listing AA graphs]`