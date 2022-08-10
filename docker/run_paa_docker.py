#!/usr/bin/env python

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
parser.add_argument("--AA_runmode", help="If --run_AA selected, set the --runmode argument to AA. Default mode is "
					"'FULL'", choices=['FULL', 'BPGRAPH', 'CYCLES', 'SVVIEW'], default='FULL')
parser.add_argument("--AA_extendmode", help="If --run_AA selected, set the --extendmode argument to AA. Default "
                    "mode is 'EXPLORE'", choices=["EXPLORE", "CLUSTERED", "UNCLUSTERED", "VIRAL"], default='EXPLORE')
parser.add_argument("--AA_insert_sdevs", help="Number of standard deviations around the insert size. May need to "
                        "increase for sequencing runs with high variance after insert size selection step. (default "
                        "3.0)", type=float, default=3.0)
parser.add_argument("--normal_bam", help="Path to matched normal bam for CNVKit (optional)", default=None)
parser.add_argument("--ploidy", type=int, help="Ploidy estimate for CNVKit (optional)", default=None)
parser.add_argument("--purity", type=float, help="Tumor purity estimate for CNVKit (optional)", default=None)
parser.add_argument("--use_CN_prefilter", help="Pre-filter CNV calls on number of copies gained above median "
                        "chromosome arm CN. Strongly recommended if input CNV calls have been scaled by purity or "
                        "ploidy. This argument is off by default but is set if --ploidy or --purity is provided for"
                        "CNVKit.", action='store_true')
parser.add_argument("--cnvkit_segmentation", help="Segmentation method for CNVKit (if used), defaults to CNVKit default"
					" segmentation method (cbs).", choices=['cbs', 'haar', 'hmm', 'hmm-tumor','hmm-germline', 'none'],
					default='cbs')
parser.add_argument("--no_filter", help="Do not run amplified_intervals.py to identify amplified seeds",
					action='store_true')
parser.add_argument("--align_only", help="Only perform the alignment stage (do not run CNV calling and seeding",
					action='store_true')
parser.add_argument("--cnv_bed", help="BED file (or CNVKit .cns file) of CNV changes. Fields in the bed file should"
									  " be: chr start end name cngain", default="")
parser.add_argument("--run_as_user", help="Run the docker image as the user launching this script. Alternatively, instead of setting this flag"
					" one can also rebuild the docker image using docker build . -t jluebeck/prepareaa:latest --build-arg set_uid=$UID --build-arg set_gid=$(id -g) ",
					action='store_true')
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("--sorted_bam", "--bam", help="Coordinate-sorted BAM file (aligned to an AA-supported reference.)")
group.add_argument("--fastqs", help="Fastq files (r1.fq r2.fq)", nargs=2)
group.add_argument("--completed_AA_runs", help="Path to a directory containing one or more completed AA runs which utilized the same reference genome.")
# group2.add_argument("--reuse_canvas", help="Start using previously generated Canvas results. Identify amplified "
# 										   "intervals immediately.", action='store_true')
# group2.add_argument("--canvas_dir", help="Path to folder with Canvas executable and \"/canvasdata\" folder "
# 										 "(reference files organized by reference name).", default="")
# group2.add_argument("--cnvkit_dir", help="Path to cnvkit.py", default="")
args = parser.parse_args()

if (args.fastqs or args.completed_AA_runs) and not args.ref:
	sys.stderr.write("Must specify --ref when providing unaligned fastq files.")
	sys.exit(1)

if not args.output_directory:
	args.output_directory = os.getcwd()

if args.output_directory == "/":
	sys.stderr.write("Output directory should not be root!\n")
	sys.exit(1)

print("making output directory read/writeable")
cmd = "chmod a+rw {} -R".format(args.output_directory)
print(cmd)
call(cmd, shell=True)

try:
	AA_REPO = os.environ['AA_DATA_REPO'] + "/"

except KeyError:
	sys.stderr.write("AA_DATA_REPO bash variable not found. AmpliconArchitect may not be properly installed.\n")
	sys.exit(1)


if not os.path.exists(os.path.join(AA_REPO, "coverage.stats")):
	print("coverage.stats file not found in " + AA_REPO + "\nCreating a new coverage.stats file.")
	cmd = "touch {}coverage.stats && chmod a+rw {}coverage.stats".format(AA_REPO, AA_REPO)
	print(cmd)
	call(cmd, shell=True)


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
cnvdir = os.path.abspath(cnvdir)

# assemble an argstring
argstring = "-t " + str(args.nthreads) + " --cngain " + str(args.cngain) + " --cnsize_min " + \
			str(args.cnsize_min) + " --downsample " + str(args.downsample) + " -s " + args.sample_name + \
			" -o /home/output" + " --AA_extendmode " + args.AA_extendmode + " --AA_runmode " + args.AA_runmode + \
			" --AA_insert_sdevs " + str(args.AA_insert_sdevs)

if args.ref:
	argstring += " --ref " + args.ref

if args.sorted_bam:
	args.sorted_bam = os.path.abspath(args.sorted_bam)
	bamdir, bamname = os.path.split(args.sorted_bam)
	argstring += " --sorted_bam /home/bam_dir/" + bamname
	if args.normal_bam:
		norm_bamdir, norm_bamname = os.path.split(args.normal_bam)
		argstring += " --normal_bam /home/norm_bam_dir/" + norm_bamname
	else:
		norm_bamdir = bamdir

elif args.fastqs:
	args.fastqs[0], args.fastqs[1] = os.path.abspath(args.fastqs[0]), os.path.abspath(args.fastqs[1])
	_, fq1name = os.path.split(args.fastqs[0])
	bamdir, fq2name = os.path.split(args.fastqs[1])
	norm_bamdir = bamdir
	argstring += " --fastqs /home/bam_dir/" + fq1name + " /home/bam_dir/" + fq2name

else:
	argstring += " --completed_AA_runs /home/bam_dir/ --completed_run_metadata None"
	bamdir = os.path.abspath(args.completed_AA_runs)
	norm_bamdir = bamdir

if args.ploidy:
	argstring += " --ploidy " + str(args.ploidy)

if args.purity:
	argstring += " --purity " + str(args.purity)

if args.use_CN_prefilter:
	argstring += " --use_CN_prefilter"

if args.cnv_bed:
	argstring += " --cnv_bed /home/bed_dir/" + cnvname
elif not args.completed_AA_runs:
	argstring += " --cnvkit_dir /home/programs/cnvkit.py"

if args.cnvkit_segmentation:
	argstring += " --cnvkit_segmentation " + args.cnvkit_segmentation

if args.no_filter:
	argstring += " --no_filter"

if args.align_only:
	argstring += " --align_only"

if args.run_AA:
	argstring += " --run_AA"

if args.run_AC:
	argstring += " --run_AC"

userstring = ""
if args.run_as_user:
	userstring = " -e HOST_UID=$(id -u) -e HOST_GID=$(id -g) -u $(id -u):$(id -g)"

print("Creating a docker script with the following argstring:")
print(argstring + "\n")
with open("paa_docker.sh", 'w') as outfile:
	outfile.write("#!/bin/bash\n\n")
	outfile.write("export argstring=\"" + argstring + "\"\n")

	# assemble a docker command string
	dockerstring = "docker run --rm" + userstring + " -e AA_DATA_REPO=/home/data_repo -e argstring=\"$argstring\"" + \
		" -v $AA_DATA_REPO:/home/data_repo -v " + bamdir + ":/home/bam_dir -v " + norm_bamdir + \
		":/home/norm_bam_dir -v " + cnvdir + ":/home/bed_dir -v " + args.output_directory + ":/home/output -v " + \
		MOSEKLM_LICENSE_FILE + ":/home/programs/mosek/8/licenses jluebeck/prepareaa bash /home/run_paa_script.sh"

	print("\n" + dockerstring + "\n")
	outfile.write(dockerstring)

call("chmod +x ./paa_docker.sh", shell=True)
call("./paa_docker.sh", shell=True)
call("rm paa_docker.sh", shell=True)
