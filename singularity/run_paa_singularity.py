#!/usr/bin/env python

import argparse
import json
import os
import subprocess
from subprocess import call
import sys

# check singularity version
def test_singularity_version():
    singularity_version = subprocess.check_output(['singularity', '--version']).decode().strip().lower().rsplit("version")[1]
    major, minor, patch = map(int, singularity_version.split('.')[0:3])
    assert (major, minor, patch) >= (3, 6, 0),'Singularity version {singularity_version} is not supported. Please upgrade to version 3.6.0 or higher.'


def metadata_helper(metadata_args):
    """
    If metadata provided, this helper is used to parse it.

    input --> Metadata Args
    output --> fp to json file of sample metadata to build on
    """
    keys = "sample_metadata,sample_type,tissue_of_origin, \
            sample_description,run_metadata_file,number_of_AA_amplicons, \
            sample_source,number_of_AA_features".split(',')

    with open(metadata_args[0], 'r') as json_file:
        json_obj = json.load(json_file)

    for key_ind in range(len(keys)):
        key = keys[key_ind]
        json_obj[key] = metadata_args[key_ind]

    with open('/home/metadata.json', 'w') as json_file:
        json.dump(json_obj, json_file, indent=4)
    json_file.close()


# Parses the command line arguments
parser = argparse.ArgumentParser(
    description="A simple pipeline wrapper for AmpliconArchitect, invoking alignment, variant calling, "
    "and CNV calling prior to AA. The CNV calling is necesary for running AA")
parser.add_argument("--sif", help="location of the ampliconsuite-pipeline.sif file.Only required if .sif file not in "
                                  "cache", type=str)
parser.add_argument("-o", "--output_directory",
                    help="output directory names (will create if not already created)")
parser.add_argument("-s", "--sample_name", help="sample name", required=True)
parser.add_argument("-t", "--nthreads",
                    help="Number of threads to use in BWA and CNV calling", required=True)
parser.add_argument("--run_AA", help="Run AA after all files prepared. Default off.", action='store_true')
parser.add_argument("--run_AC", help="Run AmpliconClassifier after all files prepared. Default off.",
                    action='store_true')
parser.add_argument("--ref", help="Reference genome version.",
                    choices=["hg19", "GRCh37", "GRCh38", "GRCh38_viral", "hg38", "mm10", "GRCm38"])
# parser.add_argument("--vcf", help="VCF (in Canvas format, i.e., \"PASS\" in filter field, AD field as 4th entry of "
# 								  "FORMAT field). When supplied with \"--sorted_bam\", pipeline will start from Canvas CNV stage."
# 					)
parser.add_argument("--cngain", type=float,
                    help="CN gain threshold to consider for AA seeding", default=4.5)
parser.add_argument("--cnsize_min", type=int, help="CN interval size (in bp) to consider for AA seeding",
                    default=50000)
parser.add_argument("--downsample", type=float,
                    help="AA downsample argument (see AA documentation)", default=10)

parser.add_argument("--AA_runmode", help="If --run_AA selected, set the --runmode argument to AA. Default mode is "
                    "'FULL'", choices=['FULL', 'BPGRAPH', 'CYCLES', 'SVVIEW'], default='FULL')
parser.add_argument("--AA_extendmode", help="If --run_AA selected, set the --extendmode argument to AA. Default "
                    "mode is 'EXPLORE'", choices=["EXPLORE", "CLUSTERED", "UNCLUSTERED", "VIRAL"], default='EXPLORE')
parser.add_argument("--AA_insert_sdevs", help="Number of standard deviations around the insert size. May need to "
                    "increase for sequencing runs with high variance after insert size selection step. (default "
                    "3.0)", type=float, default=None)
parser.add_argument(
    "--normal_bam", help="Path to matched normal bam for CNVKit (optional)", default=None)
parser.add_argument("--ploidy", type=int,
                    help="Ploidy estimate for CNVKit (optional)", default=None)
parser.add_argument("--purity", type=float,
                    help="Tumor purity estimate for CNVKit (optional)", default=None)
parser.add_argument("--cnvkit_segmentation", help="Segmentation method for CNVKit (if used), defaults to CNVKit default"
                    " segmentation method (cbs).", choices=['cbs', 'haar', 'hmm', 'hmm-tumor', 'hmm-germline', 'none'],
                    default='cbs')
parser.add_argument("--no_filter", help="Do not run amplified_intervals.py to identify amplified seeds",
                                        action='store_true')
parser.add_argument("--align_only", help="Only perform the alignment stage (do not run CNV calling and seeding",
                    action='store_true')
parser.add_argument("--cnv_bed", help="BED file (or CNVKit .cns file) of CNV changes. Fields in the bed file should"
                    " be: chr start end name cngain", default="")
# parser.add_argument("--run_as_user", help="Run the docker image as the user launching this script. Alternatively, instead of setting this flag"
#                     " one can also rebuild the docker image using docker build . -t jluebeck/prepareaa:latest --build-arg set_uid=$UID --build-arg set_gid=$(id -g) ",
#                     action='store_true')
parser.add_argument(
    "--no_QC", help="Skip QC on the BAM file.", action='store_true')
# parser.add_argument(
#     '--ref_path', help="Path to reference Genome. Won't download reference genome if provided.", default="None")
# parser.add_argument(
#     '--AA_seed', help='Seeds that sets randomness for AA', default=0)
parser.add_argument(
    '--metadata', help="Path to a JSON of sample metadata to build on", default="", nargs="+")

# parser.add_argument("--sample_metadata", help="Path to a JSON of sample metadata to build on")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("--sorted_bam", "--bam",
                   help="Coordinate-sorted BAM file (aligned to an AA-supported reference.)")
group.add_argument("--fastqs", help="Fastq files (r1.fq r2.fq)", nargs=2)
group.add_argument("--completed_AA_runs", help="Path to a directory containing one or more completed AA runs which "
                                               "utilized the same reference genome.")

args = parser.parse_args()
test_singularity_version()

if args.sif and not args.sif.endswith("ampliconsuite-pipeline.sif"):
    if not args.sif.endswith("/"): args.sif+="/"
    args.sif+="ampliconsuite-pipeline.sif"

else:
    args.sif = "ampliconsuite-pipeline.sif"

if args.ref == "hg38": args.ref = "GRCh38"
if args.ref == "GRCm38": args.ref = "mm10"
if (args.fastqs or args.completed_AA_runs) and not args.ref:
    sys.stderr.write(
        "Must specify --ref when providing unaligned fastq files.")
    sys.exit(1)

if not args.output_directory:
    args.output_directory = os.getcwd()

args.output_directory = os.path.realpath(args.output_directory)
print("Real path of output directory is set to " + args.output_directory)
if args.output_directory == "/":
    sys.stderr.write("Output directory should not be root!\n")
    sys.exit(1)

if not os.path.exists(args.output_directory):
    print("Output directory did not exist. Made output directory.")
    os.makedirs(args.output_directory)

print("making output directory read/writeable")
cmd = "chmod a+rw {} -R".format(args.output_directory)
print(cmd)
call(cmd, shell=True)

if 'AA_DATA_REPO' in os.environ:
    AA_REPO = os.environ['AA_DATA_REPO'] + "/"
    if not os.path.exists(os.path.join(AA_REPO, "coverage.stats")):
        print("coverage.stats file not found in " + AA_REPO +
              "\nCreating a new coverage.stats file.")
        cmd = "touch {}coverage.stats && chmod a+rw {}coverage.stats".format(
            AA_REPO, AA_REPO)
        print(cmd)
        call(cmd, shell=True)

else:
    AA_REPO = None
    sys.stderr.write("$AA_DATA_REPO bash variable not set. Singularity image will download large (>3 Gb) data repo each"
                     " time it is run. See installation instructions to optimize this process.\n")


if not os.path.exists(os.environ['HOME'] + "/mosek/mosek.lic"):
    sys.stderr.write(
        "Mosek license (mosek.lic) file not found in $HOME/mosek/. Please see README for instructions.\n")
    sys.exit(1)

# attach some directories
cnvdir, cnvname = os.path.split(args.cnv_bed)
cnvdir = os.path.realpath(cnvdir)

# assemble an argstring
argstring = "-t " + str(args.nthreads) + " --cngain " + str(args.cngain) + " --cnsize_min " + \
    str(args.cnsize_min) + " --downsample " + str(args.downsample) + " -s " + args.sample_name + \
    " --AA_extendmode " + args.AA_extendmode + " --AA_runmode " + args.AA_runmode

if args.ref:
    argstring += " --ref " + args.ref

if args.sorted_bam:
    args.sorted_bam = os.path.realpath(args.sorted_bam)
    bamdir, bamname = os.path.split(args.sorted_bam)
    norm_bamdir = bamdir
    argstring += " --sorted_bam /home/bam_dir/" + bamname

elif args.fastqs:
    args.fastqs[0], args.fastqs[1] = os.path.realpath(
        args.fastqs[0]), os.path.realpath(args.fastqs[1])
    _, fq1name = os.path.split(args.fastqs[0])
    bamdir, fq2name = os.path.split(args.fastqs[1])
    norm_bamdir = bamdir
    argstring += " --fastqs /home/bam_dir/" + fq1name + " /home/bam_dir/" + fq2name

else:
    argstring += " --completed_AA_runs /home/bam_dir/ --completed_run_metadata None"
    bamdir = os.path.realpath(args.completed_AA_runs)
    norm_bamdir = bamdir

if args.normal_bam:
    args.normal_bam = os.path.realpath(args.normal_bam)
    norm_bamdir, norm_bamname = os.path.split(args.normal_bam)
    argstring += " --normal_bam /home/norm_bam_dir/" + norm_bamname

if args.sv_vcf:
    vcf_dir, vcf_name = os.path.split(args.sv_vcf)
    argstring += " --sv_vcf /home/vcf_dir/" + vcf_name
    if args.sv_vcf_no_filter:
        argstring += " --sv_vcf_no_filter"

else:
    vcf_dir = bamdir

if args.ploidy:
    argstring += " --ploidy " + str(args.ploidy)

if args.purity:
    argstring += " --purity " + str(args.purity)

if args.cnv_bed:
    argstring += " --cnv_bed /home/bed_dir/" + cnvname

elif args.align_only:
    argstring += " --align_only"

elif not args.completed_AA_runs:
    argstring += " --cnvkit_dir /home/programs/cnvkit.py"

if args.cnvkit_segmentation:
    argstring += " --cnvkit_segmentation " + args.cnvkit_segmentation

if args.no_filter:
    argstring += " --no_filter"

if args.no_QC:
    argstring += " --no_QC"

if args.AA_insert_sdevs:
    argstring += " --AA_insert_sdevs " + str(args.AA_insert_sdevs)

# To use, would need to mount the directory of this file. Users should just modify as needed afterwards.
# if args.sample_metadata:
# 	args.sample_metadata = os.path.abspath(args.sample_metadata)
# 	argstring += " --sample_metadata " + args.sample_metadata

if args.run_AA:
    argstring += " --run_AA"

if args.run_AC:
    argstring += " --run_AC"

if args.metadata != "":
    metadata_helper(args.metadata)
    argstring += " --sample_metadata /home/metadata.json"
#
# os.environ['AA_SEED'] = str(args.AA_seed)

# userstring = ""
# if args.run_as_user:
#     userstring = " -e HOST_UID=$(id -u) -e HOST_GID=$(id -g) -u $(id -u):$(id -g)"

env_outname = "paa_envs_" + args.sample_name + ".txt"
runscript_outname = "paa_singularity_" + args.sample_name + ".sh"

with open(runscript_outname, 'w') as outfile:
    outfile.write("#!/bin/bash\n\n")
    outfile.write("export argstring=\"" + argstring + "\"\n")
    outfile.write("export SAMPLE_NAME=" + args.sample_name + "\n")

    # Download the reference genome if necessary
    no_data_repo = not AA_REPO or (args.ref and not os.path.exists(AA_REPO + args.ref))
    if no_data_repo and args.ref:
        outfile.write('echo DOWNLOADING {} NOW ....\n'.format(args.ref))
        data_repo_d = args.output_directory + '/data_repo'
        outfile.write('mkdir -p ' + data_repo_d + '\n')
        outfile.write('export AA_DATA_REPO=' + data_repo_d + '\n')
        if args.fastqs:
            outfile.write(
                'wget -q -P $AA_DATA_REPO https://datasets.genepattern.org/data/module_support_files/AmpliconArchitect/{}_indexed.tar.gz\n'.format(args.ref))
            outfile.write(
                'wget -q -P $AA_DATA_REPO https://datasets.genepattern.org/data/module_support_files/AmpliconArchitect/{}_indexed_md5sum.txt\n'.format(args.ref))
            outfile.write(
                'tar zxf $AA_DATA_REPO/{}_indexed.tar.gz --directory $AA_DATA_REPO\n'.format(args.ref))

        else:
            outfile.write(
                'wget -q -P $AA_DATA_REPO https://datasets.genepattern.org/data/module_support_files/AmpliconArchitect/{}.tar.gz\n'.format(
                    args.ref))
            outfile.write(
                'wget -q -P $AA_DATA_REPO https://datasets.genepattern.org/data/module_support_files/AmpliconArchitect/{}_md5sum.txt\n'.format(
                    args.ref))
            outfile.write(
                'tar zxf $AA_DATA_REPO/{}.tar.gz --directory $AA_DATA_REPO\n'.format(args.ref))

        outfile.write(
            'touch $AA_DATA_REPO/coverage.stats && chmod a+rw $AA_DATA_REPO/coverage.stats\n')
        outfile.write('echo DOWNLOADING {} COMPLETE\n'.format(args.ref))

    elif no_data_repo and not args.ref:
        sys.stderr.write("Must specify --ref argument!")
        sys.exit(1)

    # write exported envs to env-file
    with open(env_outname, 'w') as env_file:
        env_file.write('argstring="' + argstring + '"\n')
        env_file.write("SAMPLE_NAME=" + args.sample_name)

    # assemble a singularity command string
    sing_string = "singularity exec --no-home --cleanenv --env-file " + env_outname + " --bind $AA_DATA_REPO:" \
                  "/home/data_repo --bind " + bamdir + ":/home/bam_dir --bind " + norm_bamdir + ":/home/norm_bam_dir " \
                  "--bind " + cnvdir + ":/home/bed_dir --bind " + args.output_directory + ":/home/output --bind " \
                  "$HOME/mosek:/home/mosek/ --bind " + vcf_dir + ":/home/vcf_dir " + args.sif + " bash /home/run_paa_script.sh "

    print("\n" + sing_string + "\n")
    outfile.write(sing_string)

outfile.close()

call("chmod +x ./" + runscript_outname, shell=True)
call("./" + runscript_outname, shell=True)
call("rm " + runscript_outname, shell=True)
call("rm -f " + env_outname, shell=True)
if no_data_repo:
    cmd = "rm -rf " + data_repo_d
    print("Cleaning up data repo")
    print(cmd)
    call(cmd, shell=True)
