#!/usr/bin/env python2

import sys
import os
import argparse
from subprocess import call


# Parses the command line arguments
parser = argparse.ArgumentParser(
	description="A simple pipeline wrapper for AmpliconArchitect, invoking alignment, variant calling, "
				"and CNV calling prior to AA. The CNV calling is necesary for running AA")
parser.add_argument("-o", "--output_directory", help="output directory names (will create if not already created)")
parser.add_argument("-s", "--sample_name", help="sample name", required=True)
parser.add_argument("-t", "--nthreads", help="Number of threads to use in BWA and CNV calling", required=True)
parser.add_argument("--run_AA", help="Run AA after all files prepared. Default off.", action='store_true')
parser.add_argument("--run_AC", help="Run AmpliconClassifier after all files prepared. Default off.",
					action='store_true')
parser.add_argument("--ref", help="Reference genome version.", choices=["hg19", "GRCh37", "GRCh38", "hg38", "mm10",
																		"GRCm38"])
# parser.add_argument("--vcf", help="VCF (in Canvas format, i.e., \"PASS\" in filter field, AD field as 4th entry of "
# 								  "FORMAT field). When supplied with \"--sorted_bam\", pipeline will start from Canvas CNV stage."
# 					)
parser.add_argument("--cngain", type=float, help="CN gain threshold to consider for AA seeding", default=4.5)
parser.add_argument("--cnsize_min", type=int, help="CN interval size (in bp) to consider for AA seeding",
					default=50000)
parser.add_argument("--downsample", type=float, help="AA downsample argument (see AA documentation)", default=10)
# parser.add_argument("--use_old_samtools", help="Indicate you are using an old build of samtools (prior to version "
# 											   "1.0)", action='store_true', default=False)
# parser.add_argument("--rscript_path", help="Specify custom path to Rscript, if needed when using CNVKit "
# 										   "(which requires R version >3.4)")
# parser.add_argument("--python3_path", help="Specify custom path to python3, if needed when using CNVKit (requires "
# 										   "python3)")
# parser.add_argument("--freebayes_dir", help="Path to directory where freebayes executable exists (not the path to "
# 											"the executable itself). Only needed for Canvas and freebayes is not installed on system path.",
# 					default=None)
# parser.add_argument("--aa_data_repo", help="Specify a custom $AA_DATA_REPO path FOR PRELIMINARY STEPS ONLY(!). Will"
# 										   " not override bash variable during AA")
# parser.add_argument("--aa_src", help="Specify a custom $AA_SRC path. Overrides the bash variable")
# parser.add_argument("--AA_runmode", help="If --run_AA selected, set the --runmode argument to AA. Default mode is "
# 										 "'FULL'", choices=['FULL', 'BPGRAPH', 'CYCLES', 'SVVIEW'], default='FULL')
parser.add_argument("--normal_bam", help="Path to matched normal bam for CNVKit (optional)", default=None)
parser.add_argument("--ploidy", type=int, help="Ploidy estimate for CNVKit (optional)", default=None)
parser.add_argument("--purity", type=float, help="Tumor purity estimate for CNVKit (optional)", default=None)
parser.add_argument("--cnvkit_segmentation", help="Segmentation method for CNVKit (if used), defaults to CNVKit "
												  "default segmentation method (cbs).",
					choices=['cbs', 'haar', 'hmm', 'hmm-tumor','hmm-germline', 'none'], default='cbs')
parser.add_argument("--no_filter", help="Do not run amplified_intervals.py to identify amplified seeds",
					action='store_true')
parser.add_argument("--cnv_bed", help="BED file (or CNVKit .cns file) of CNV changes. Fields in the bed file should"
									  " be: chr start end name cngain", default="")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("--sorted_bam", "--bam", help="Coordinate-sorted BAM file (aligned to an AA-supported reference.)")
group.add_argument("--fastqs", help="Fastq files (r1.fq r2.fq)", nargs=2)

# group2.add_argument("--reuse_canvas", help="Start using previously generated Canvas results. Identify amplified "
# 										   "intervals immediately.", action='store_true')

# group2.add_argument("--canvas_dir", help="Path to folder with Canvas executable and \"/canvasdata\" folder "
# 										 "(reference files organized by reference name).", default="")
# group2.add_argument("--cnvkit_dir", help="Path to cnvkit.py", default="")


args = parser.parse_args()

if args.fastqs and not args.ref:
	sys.stderr.write("Must specify --ref when providing unaligned fastq files.")
	sys.exit(1)

if not args.output_directory:
	args.output_directory = os.getcwd()

try:
	AA_REPO = os.environ['AA_DATA_REPO'] + "/"

except KeyError:
	sys.stderr.write("AA_DATA_REPO bash variable not found. AmpliconArchitect may not be properly installed.\n")
	sys.exit(1)

# try:
# 	AA_SRC = os.environ['AA_SRC']

# except KeyError:
# 	sys.stderr.write("AA_SRC bash variable not found. AmpliconArchitect may not be properly installed.\n")
# 	sys.exit(1)


try:
	# MOSEK LICENSE FILE PATH
	MOSEKLM_LICENSE_FILE = os.environ['MOSEKLM_LICENSE_FILE']
	if not os.path.exists(MOSEKLM_LICENSE_FILE + "/mosek.lic"):
		raise KeyError

except KeyError:
	sys.stderr.write("Mosek license (.lic) file not found. AmpliconArchitect may not be properly installed.\n")
	sys.exit(1)


# attach some directories
cnvdir, cnvname = os.path.split(args.cnv_bed)

# assemble an argstring
argstring = "--ref " + args.ref + " -t " + str(args.nthreads) + " --cngain " + str(args.cngain) + " --cnsize_min " + \
			str(args.cnsize_min) + " --downsample " + str(args.downsample) + " -s " + args.sample_name + \
			" -o /home/output"

if args.sorted_bam:
	bamdir, bamname = os.path.split(args.sorted_bam)
	argstring += " --sorted_bam /home/bam_dir/" + bamname
	if args.normal_bam:
		norm_bamdir, norm_bamname = os.path.split(args.normal_bam)
		argstring += " --normal_bam /home/norm_bam_dir/" + norm_bamname
	else:
		norm_bamdir = bamdir

else:
	_, fq1name = os.path.split(args.fastqs[0])
	bamdir, fq2name = os.path.split(args.fastqs[1])
	norm_bamdir = bamdir
	argstring += " --fastqs /home/bam_dir/" + fq1name + " /home/bam_dir/" + fq2name

if args.ploidy:
	argstring += " --ploidy " + str(args.ploidy)

if args.purity:
	argstring += " --purity " + str(args.purity)

if args.cnv_bed:
	argstring += " --cnv_bed /home/bed_dir/" + cnvname
else:
	argstring += " --cnvkit_dir /home/programs/cnvkit.py"

if args.cnvkit_segmentation:
	argstring += " --cnvkit_segmentation " + args.cnvkit_segmentation

if not args.run_AA:
	argstring += " --run_AA"

if not args.run_AC:
	argstring += " --run_AC"

print("Creating a docker script with the following argstring:")
print(argstring + "\n")
with open("paa_docker.sh", 'w') as outfile:
	outfile.write("#!/bin/bash\n\n")
	outfile.write("export argstring=\"" + argstring + "\"\n")

	# assemble a docker command string
	dockerstring = "docker run --rm -e AA_DATA_REPO=/home/data_repo -e argstring=\"$argstring\"" + \
		" -v $AA_DATA_REPO:/home/data_repo -v " + bamdir + ":/home/bam_dir -v " + norm_bamdir + \
		":/home/norm_bam_dir -v " + cnvdir + ":/home/bed_dir -v " + args.output_directory + ":/home/output -v " + \
		MOSEKLM_LICENSE_FILE + ":/home/programs/mosek/8/licenses jluebeck/prepareaa bash /home/run_paa_script.sh"

	print("\n" + dockerstring + "\n")
	outfile.write(dockerstring)

call("chmod +x ./paa_docker.sh", shell=True)
call("./paa_docker.sh", shell=True)
call("rm paa_docker.sh", shell=True)
