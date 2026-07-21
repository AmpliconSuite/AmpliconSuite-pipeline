#!/usr/bin/env python3

# This python script is the driver for launching the AmpliconSuite-pipeline singularity image

import argparse
import json
import os
import re
import subprocess
from subprocess import call
import sys


DATA_REPO_BASE_URL = "https://refs.ampliconrepository.org/data/module_support_files/AmpliconArchitect/"


# check singularity version
def test_singularity_version():
    version_output = subprocess.check_output(['singularity', '--version']).decode().strip().lower()
    version_match = re.search(
        r'\b(apptainer|singularity(?:-ce)?)\s+(?:version\s+)?(\d+)\.(\d+)',
        version_output,
    )
    if not version_match:
        sys.stderr.write(
            'Could not identify a Singularity/Apptainer version from: {}\n'.format(version_output))
        sys.exit(1)

    runtime, major, minor = version_match.groups()
    version = (int(major), int(minor))
    if runtime == 'apptainer' and version < (1, 0):
        sys.stderr.write(
            'Apptainer version {}.{} is not supported. Please upgrade to version 1.0 or higher.\n'.format(
                major, minor))
        sys.exit(1)
    if runtime != 'apptainer' and version < (3, 6):
        sys.stderr.write(
            'Singularity version {}.{} is not supported. Please upgrade to version 3.6 or higher.\n'.format(
                major, minor))
        sys.exit(1)


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
    description="A simple pipeline wrapper for AmpliconArchitect and AmpliconClassifier, invoking alignment "
                "and CNV calling prior to AA.")
parser.add_argument("--sif", help="Path of the ampliconsuite-pipeline.sif file.", type=str, required=True)
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
parser.add_argument("--cngain", type=float,
                    help="CN gain threshold to consider for AA seeding", default=4.5)
parser.add_argument("--cnsize_min", type=int, help="CN interval size (in bp) to consider for AA seeding",
                    default=50000)
parser.add_argument("--no_cstats",
                    help="Do not read from or write to the $AA_DATA_REPO/coverage.stats file. (default not set).",
                    action='store_true')
parser.add_argument("--downsample", type=float,
                    help="AA downsample argument (see AA documentation)", default=10)

parser.add_argument("--AA_runmode", help="If --run_AA selected, set the --runmode argument to AA. Default mode is "
                    "'FULL'", choices=['FULL', 'BPGRAPH', 'CYCLES', 'SVVIEW'], default='FULL')
parser.add_argument("--AA_extendmode", help="If --run_AA selected, set the --extendmode argument to AA. Default "
                    "mode is 'EXPLORE'", choices=["EXPLORE", "CLUSTERED", "UNCLUSTERED", "VIRAL"], default='EXPLORE')
parser.add_argument("--AA_insert_sdevs", help="Number of standard deviations around the insert size. May need to "
                    "increase for sequencing runs with high variance after insert size selection step. (default "
                    "3.0)", type=float, default=None)
parser.add_argument('--pair_support_min', dest='pair_support_min', help="Number of read pairs for "
                        "minimum breakpoint support (default 2 but typically becomes higher due to coverage-scaled "
                        "cutoffs)", metavar='INT', action='store', type=int)
parser.add_argument('--foldback_pair_support_min', help="Number of read pairs for minimum foldback SV support "
                        "(default 2 but typically becomes higher due to coverage-scaled cutoffs). Used value will be the maximum"
                        " of pair_support and this argument. Raising to 3 will help dramatically in heavily artifacted samples.",
                        metavar='INT', action='store', type=int)
parser.add_argument("--AA_solver", help="If --run_AA selected, set the copy-number optimizer AA uses via its --solver "
                    "argument. If 'mosek' (default) is unavailable, the pipeline selects 'clarabel'; AA also retries "
                    "with Clarabel if Mosek fails at runtime.", choices=['mosek', 'clarabel'], default='mosek')
parser.add_argument(
    "--normal_bam", help="Path to matched normal bam for CNVKit (optional)", default=None)
parser.add_argument("--ploidy", type=int,
                    help="Ploidy estimate for CNVKit (optional)", default=None)
parser.add_argument("--purity", type=float,
                    help="Tumor purity estimate for CNVKit (optional)", default=None)
parser.add_argument("--cnvkit_segmentation", help="Segmentation method for CNVKit (if used), defaults to CNVKit default"
                    " segmentation method (cbs).", choices=['cbs', 'haar', 'hmm', 'hmm-tumor', 'hmm-germline', 'none'],
                    default='cbs')
parser.add_argument("--no_filter", help="Do not run amplified_intervals.py to remove low confidence candidate seed"
                                            " regions overlapping repetitive parts of the genome", action='store_true')
parser.add_argument("--align_only", help="Only perform the alignment stage (do not run CNV calling and seeding",
                    action='store_true')
parser.add_argument("--cnv_bed", "--bed", help="BED file (or CNVKit .cns file) of CNV changes. Fields in the bed file should"
                    " be: chr start end name cngain", default="")
parser.add_argument(
    "--no_QC", help="Skip QC on the BAM file.", action='store_true')
parser.add_argument("--sv_vcf", help="Provide a VCF file of externally-called SVs to augment SVs identified by AA internally.",
                    metavar='FILE', action='store', type=str)
parser.add_argument("--sv_vcf_no_filter", help="Use all external SV calls from the --sv_vcf arg, even "
                    "those without 'PASS' in the FILTER column.", action='store_true', default=False)
parser.add_argument('--sv_vcf_include_sr',
                    help="Include single-ended reads when counting support for an SV in the provided VCF",
                    action='store_true', default=False)
parser.add_argument('--metadata', help="Path to a JSON of sample metadata to build on", default="", nargs="+")

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("--bam", "--sorted_bam",
                   help="Coordinate-sorted BAM file (aligned to an AA-supported reference.)")
group.add_argument("--fastqs", help="Fastq files (r1.fq r2.fq)", nargs=2)
group.add_argument("--completed_AA_runs", help="Path to a directory containing one or more completed AA runs which "
                                               "utilized the same reference genome.")

args = parser.parse_args()
test_singularity_version()

if args.sif and not args.sif.endswith(".sif"):
    sys.stderr.write("Path of .sif file must go to .sif file!\n")
    sys.exit(1)

if args.ref == "hg38": args.ref = "GRCh38"
if args.ref == "GRCm38": args.ref = "mm10"
if (args.fastqs or args.completed_AA_runs) and not args.ref:
    sys.stderr.write(
        "Must specify --ref when providing unaligned fastq files.\n")
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
    if not args.no_cstats and not os.path.exists(os.path.join(AA_REPO, "coverage.stats")):
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


# Mosek is optional. BFBArchitect may also use it under --run_AC, so an available license is still mounted for AC
# when AA explicitly uses Clarabel; only an AA request for Mosek should produce an AA fallback warning.
MOSEK_DIR = None
mosek_requested_for_aa = args.run_AA and args.AA_solver == "mosek"
if mosek_requested_for_aa or args.run_AC:
    mosek_candidate = os.environ.get('MOSEKLM_LICENSE_FILE')
    if mosek_candidate and mosek_candidate.endswith("mosek.lic"):
        sys.stderr.write(
            "MOSEKLM_LICENSE_FILE should be the path of the directory of the license, not the full path. Please update your .bashrc, and run 'source ~/.bashrc'\n")
        mosek_candidate = None
    if mosek_candidate and os.path.exists(os.path.join(mosek_candidate, "mosek.lic")):
        MOSEK_DIR = os.path.realpath(mosek_candidate)
    else:
        mosek_candidate = os.path.join(os.environ['HOME'], "mosek")
        if os.path.exists(os.path.join(mosek_candidate, "mosek.lic")):
            MOSEK_DIR = mosek_candidate

if mosek_requested_for_aa and MOSEK_DIR is None:
    sys.stderr.write(
        "Mosek is unavailable; AmpliconSuite-pipeline will use the license-free 'clarabel' solver for AA. "
        "If Mosek fails later during optimization, AA also retries with Clarabel.\n")

# Gurobi is an optional BFBArchitect accelerator under --run_AC. Discover it only when AC will run, and remain
# silent when it is absent because the built-in solver is the normal choice for most analyses.
GRB_LICENSE_FILE = None
if args.run_AC:
    for _grb_cand in [os.environ.get('GRB_LICENSE_FILE'), os.path.join(os.environ['HOME'], "gurobi.lic")]:
        if _grb_cand and os.path.isfile(_grb_cand):
            GRB_LICENSE_FILE = os.path.realpath(_grb_cand)
            break

# attach some directories
cnvdir, cnvname = os.path.split(args.cnv_bed)
cnvdir = os.path.realpath(cnvdir)

# assemble an argstring
argstring = "-t " + str(args.nthreads) + " --cngain " + str(args.cngain) + " --cnsize_min " + \
    str(args.cnsize_min) + " --downsample " + str(args.downsample) + " -s " + args.sample_name + \
    " --AA_extendmode " + args.AA_extendmode + " --AA_runmode " + args.AA_runmode

if args.ref:
    argstring += " --ref " + args.ref

if args.bam:
    args.bam = os.path.realpath(args.bam)
    bamdir, bamname = os.path.split(args.bam)
    norm_bamdir = bamdir
    argstring += " --bam /home/bam_dir/" + bamname

elif args.fastqs:
    args.fastqs[0], args.fastqs[1] = os.path.realpath(
        args.fastqs[0]), os.path.realpath(args.fastqs[1])
    _, fq1name = os.path.split(args.fastqs[0])
    bamdir, fq2name = os.path.split(args.fastqs[1])
    norm_bamdir = bamdir
    argstring += " --fastqs /home/bam_dir/" + fq1name + " /home/bam_dir/" + fq2name

else:
    # Resolve the full path first to handle symlinks
    resolved_input_path = os.path.realpath(args.completed_AA_runs)

    # Check if completed_AA_runs is a file or directory on the host
    if os.path.isfile(resolved_input_path):
        # It's a file - mount parent directory and specify the file path in container
        host_dir = os.path.dirname(resolved_input_path)
        filename = os.path.basename(resolved_input_path)
        container_file_path = "/home/bam_dir/{}".format(filename)

        argstring += " --completed_AA_runs {} --completed_run_metadata None".format(container_file_path)
        bamdir = host_dir  # Already resolved, no symlinks

    else:
        # It's a directory - mount the whole directory
        argstring += " --completed_AA_runs /home/bam_dir/ --completed_run_metadata None"
        bamdir = resolved_input_path  # Already resolved, no symlinks

    norm_bamdir = bamdir

if args.normal_bam:
    args.normal_bam = os.path.realpath(args.normal_bam)
    norm_bamdir, norm_bamname = os.path.split(args.normal_bam)
    argstring += " --normal_bam /home/norm_bam_dir/" + norm_bamname

if args.sv_vcf:
    args.sv_vcf = os.path.realpath(args.sv_vcf)
    vcf_dir, vcf_name = os.path.split(args.sv_vcf)
    argstring += " --sv_vcf /home/vcf_dir/" + vcf_name
    if args.sv_vcf_no_filter:
        argstring += " --sv_vcf_no_filter"
    if args.sv_vcf_include_sr:
        argstring += " --sv_vcf_include_sr"

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

if args.no_cstats:
    argstring += " --no_cstats"

if args.no_QC:
    argstring += " --no_QC"

if args.AA_insert_sdevs:
    argstring += " --AA_insert_sdevs " + str(args.AA_insert_sdevs)

if args.pair_support_min:
    argstring += " --pair_support_min " + str(args.pair_support_min)

if args.foldback_pair_support_min:
    argstring += " --foldback_pair_support_min " + str(args.foldback_pair_support_min)

if args.AA_solver:
    argstring += " --AA_solver " + args.AA_solver

if args.run_AA:
    argstring += " --run_AA"

if args.run_AC:
    argstring += " --run_AC"

if args.metadata != "":
    metadata_helper(args.metadata)
    argstring += " --sample_metadata /home/metadata.json"

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
                'wget -q -P $AA_DATA_REPO {}{}_indexed.tar.gz\n'.format(DATA_REPO_BASE_URL, args.ref))
            outfile.write(
                'wget -q -P $AA_DATA_REPO {}{}_indexed_md5sum.txt\n'.format(DATA_REPO_BASE_URL, args.ref))
            outfile.write(
                'tar zxf $AA_DATA_REPO/{}_indexed.tar.gz --directory $AA_DATA_REPO\n'.format(args.ref))

        else:
            outfile.write(
                'wget -q -P $AA_DATA_REPO {}{}.tar.gz\n'.format(DATA_REPO_BASE_URL, args.ref))
            outfile.write(
                'wget -q -P $AA_DATA_REPO {}{}_md5sum.txt\n'.format(DATA_REPO_BASE_URL, args.ref))
            outfile.write(
                'tar zxf $AA_DATA_REPO/{}.tar.gz --directory $AA_DATA_REPO\n'.format(args.ref))

        if not args.no_cstats:
            outfile.write(
                'touch $AA_DATA_REPO/coverage.stats && chmod a+rw $AA_DATA_REPO/coverage.stats\n')
        outfile.write('echo DOWNLOADING {} COMPLETE\n'.format(args.ref))

    elif no_data_repo and not args.ref:
        sys.stderr.write("Must specify --ref argument!\n")
        sys.exit(1)

    # write exported envs to env-file
    with open(env_outname, 'w') as env_file:
        env_file.write('argstring="' + argstring + '"\n')
        env_file.write("SAMPLE_NAME=" + args.sample_name)
        # Only export GRB_LICENSE_FILE when a license is actually mounted, so gurobipy never sees a missing path.
        if GRB_LICENSE_FILE:
            env_file.write("\nGRB_LICENSE_FILE=/home/gurobi/gurobi.lic")

    # Optional license binds. Mosek and Gurobi are both optional (AA falls back to clarabel, BFBArchitect to CBC),
    # so only bind a license when it is actually present on the host.
    license_binds = ""
    if MOSEK_DIR:
        license_binds += "--bind " + MOSEK_DIR + ":/home/mosek/:ro "
    if GRB_LICENSE_FILE:
        license_binds += "--bind " + GRB_LICENSE_FILE + ":/home/gurobi/gurobi.lic:ro "

    # assemble a singularity command string
    sing_string = "singularity exec --no-home --cleanenv --env-file " + env_outname + " --bind $AA_DATA_REPO:" \
                  "/home/data_repo --bind " + bamdir + ":/home/bam_dir --bind " + norm_bamdir + ":/home/norm_bam_dir " \
                  "--bind " + cnvdir + ":/home/bed_dir --bind " + args.output_directory + ":/home/output " \
                  + license_binds + "--bind " + vcf_dir + ":/home/vcf_dir " + args.sif + " bash /home/internal_singularity_script.sh "

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
