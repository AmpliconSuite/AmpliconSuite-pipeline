#!/usr/bin/env python2

import sys
import os
import argparse
from subprocess import call

# Parses the command line arguments
parser = argparse.ArgumentParser(description="A dockerized pipeline wrapper for AmpliconArchitect. Currently does seed filtering + AA.")
parser.add_argument("-o", "--output_directory", help="output directory names (will create if not already created)")
parser.add_argument("-s", "--sample_name", help="sample name", required=True)
parser.add_argument("-t","--nthreads",help="Number of threads to use in BWA and CNV calling",required=True)
parser.add_argument("--no_AA", help="Do not run AA after all files prepared. Default - run_AA.", action='store_true')
parser.add_argument("--ref", help="Reference genome version.",choices=["hg19","GRCh37","GRCh38"],default="hg19")
# parser.add_argument("--vcf", help="VCF (in Canvas format, i.e., \"PASS\" in filter field, AD field as 4th entry of FORMAT field). When supplied with \"--sorted_bam\", pipeline will start from Canvas CNV stage.")
parser.add_argument("--cngain",type=float,help="CN gain threshold to consider for AA seeding",default=4.999999)
parser.add_argument("--cnsize_min",type=int,help="CN interval size (in bp) to consider for AA seeding",default=50000)
parser.add_argument("--downsample",type=float,help="AA downsample argument (see AA documentation)",default=10)
# parser.add_argument("--use_old_samtools",help="Indicate you are using an old build of samtools (prior to version 1.0)",action='store_true',default=False)
# parser.add_argument("--rscript_path",help="Specify custom path to Rscript, if needed when using CNVKit (which requires R version >3.4)")
# parser.add_argument("--python3_path",help="Specify custom path to python3, if needed when using CNVKit (requires python3)")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("--sorted_bam", help= "Sorted BAM file (aligned to an AA-supported reference.)")
group.add_argument("--fastqs", help="Fastq files (r1.fq r2.fq)", nargs=2)
group2 = parser.add_mutually_exclusive_group(required=True)
# group2.add_argument("--reuse_canvas", help="Start using previously generated Canvas results. Identify amplified intervals immediately.",action='store_true')
group2.add_argument("--cnv_bed",help="BED file of CNV changes. Fields in the bed file should be: chr start end name cngain",default="")
# group2.add_argument("--canvas_dir",help="Path to folder with Canvas executable and \"/canvasdata\" folder (reference files organized by reference name).",default="")
group2.add_argument("--run_cnvkit",help="Use CNVKit 0.9.7 installed in the docker image to generate CNV calls and seed bed file",
					action='store_true')

args = parser.parse_args()

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
	#MOSEK LICENSE FILE PATH
	MOSEKLM_LICENSE_FILE = os.environ['MOSEKLM_LICENSE_FILE']
	if not os.path.exists(MOSEKLM_LICENSE_FILE + "/mosek.lic"):
		raise KeyError

except KeyError:
	sys.stderr.write("Mosek license (.lic) file not found. AmpliconArchitect may not be properly installed.\n")
	sys.exit(1)


#attach some directories
cnvdir,cnvname = os.path.split(args.cnv_bed)

#assemble an argstring
argstring = "--ref " + args.ref + " -t 1 --cngain " + str(args.cngain) + " --cnsize_min " + str(args.cnsize_min) + \
" --downsample " + str(args.downsample) + " -s " + args.sample_name

if args.sorted_bam:
	bamdir, bamname = os.path.split(args.sorted_bam)
	argstring+=" --sorted_bam /home/bam_dir/" + args.sorted_bam

else:
	bamdir, fq1name = os.path.split(args.fastqs[0])
	argstring+=" --fastqs /home/bam_dir/" + args.fastqs[0] + " /home/bam_dir/" + args.fastqs[1]

if args.cnv_bed:
	argstring+=" --cnv_bed /home/bed_dir/" + cnvname + " -o /home/output"

else:
	print("Specifying use of CNVKit from docker image")
	argstring+=" --cnvkit_dir $CNVKIT_PATH" 


if not args.no_AA:
	argstring+=" --run_AA" 

print("Creating a docker script with the following argstring:")
print(argstring + "\n")
with open("paa_docker.sh",'w') as outfile:
	outfile.write("#!/bin/bash\n\n")
	outfile.write("export argstring=\"" + argstring + "\"\n")

	#assemble a docker command string
	dockerstring = "docker run --rm -e AA_DATA_REPO=/home/data_repo -e argstring=\"$argstring\""+ \
	" -v $AA_DATA_REPO:/home/data_repo -v " + bamdir + ":/home/bam_dir -v " + cnvdir + ":/home/bed_dir -v " + \
	args.output_directory + ":/home/output -v " + MOSEKLM_LICENSE_FILE + \
	":/home/programs/mosek/8/licenses jluebeck/prepareaa bash /home/run_paa_script.sh"

	print("\n" + dockerstring + "\n")
	outfile.write(dockerstring)

call("chmod +x ./paa_docker.sh",shell=True)
call("./paa_docker.sh",shell=True)
call("rm paa_docker.sh",shell=True)




