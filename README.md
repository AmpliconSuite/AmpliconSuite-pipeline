## PrepareAA

A quickstart tool for [AmpliconArchitect](https://github.com/virajbdeshpande/AmpliconArchitect). Performs all preliminary steps (alignment, CNV calling, seed interval detection) required prior to running AmpliconArchitect.

### Prerequisites:
Requires the following tools to be installed beforehand:
- bwa mem
- samtools
- [freebayes](https://github.com/ekg/freebayes) (unless supplying your own CNV calls)
- [Canvas](https://github.com/Illumina/canvas) (unless supplying your own CNV calls)
- [AmpliconArchitect](https://github.com/virajbdeshpande/AmpliconArchitect)

Canvas can be tricky to install. Unfortunately PrepareAA only supports hg19 currently, as the standalone version of Canvas currently only seems to support hg19.

Please make sure the Canvas hg19 reference genome files are located in the expected location for canvas. To do this, you can follow instructions on the Canvas Github page or more simply, extract canvasdata.tar.gz (available here: https://drive.google.com/open?id=1Wzk7wE6Mk-k8X3XqvZziLySDWyT-tTen) in the folder with the Canvas executable, to create the canvasdata/ directory.


### Installation
Files in hg19/ folder must be placed into $AA_DATA_REPO/hg19 folder prior to using.

### Usage
If using your own CNV calls:
```
/path/to/PrepareAA/PrepareAA.py -s sample_name  -t number_of_threads  --cnv_bed your_cnvs.bed  --fastqs sample_r1.fastq sample_r2.fastq
```
Where the bed file is formatted like:
chrN    start        end        some_arbitrary_name   copy_number

If instead using Canvas for the CNV calls:
```
/path/to/PrepareAA/PrepareAA.py -s sample_name  -t number_of_threads  --canvas_lib_dir /path/to/canvas/canvas_data_dir --fastqs sample_r1.fastq sample_r2.fastq
```

If you already have your coordinate-sorted bam file, `--fastqs` can be replaced with `--sorted_bam`.



