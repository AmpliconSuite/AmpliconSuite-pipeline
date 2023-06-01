#!/usr/bin/env python3

# author: Jens Luebeck (jluebeck [at] ucsd.edu)

import argparse
from datetime import datetime
import json
import os
import random
from subprocess import *
import sys
import time
import threading

PAA_PATH = os.path.dirname(os.path.realpath(__file__)) + "/PrepareAA.py"


def generate_individual_seeds(cmd_dict, aa_py, odir, cnv_bed_dict):
    individual_seed_dct = {}
    print('Generating individual seeds')
    for sname, argstring in cmd_dict.items():
        with open(sname + "_CNV_out.txt", 'w') as outfile:
            cmd = '{} {}{}'.format(aa_py, PAA_PATH, argstring)
            print(sname)
            print(cmd + "\n")
            call(cmd, stdout=outfile, stderr=outfile, shell=True)

        # if it was a seeds file, PAA won't modify, so move it into the right location
        if sname in cnv_bed_dict and cnv_bed_dict[sname].endswith("AA_CNV_SEEDS.bed"):
            cmd = "cp {} {}/".format(cnv_bed_dict[sname], odir)
            call(cmd, shell=True)

        # store the name of the path of the seeds file
        individual_seed_dct[sname] = '{}/{}_AA_CNV_SEEDS.bed'.format(odir, sname)

    return individual_seed_dct


def group_seeds(individual_seed_dct, odir):
    samplist = list(individual_seed_dct.keys())
    outname = odir + "_".join(samplist[:2])
    if len(samplist) > 2:
        outname += "_etc_n" + str(len(samplist))

    outname+="_merged_AA_CNV_SEEDS.bed"

    bedlist = " ".join(individual_seed_dct.values())
    print("Merging seeds")
    cmd = "sort -k1,1 -k2,2n {} | bedtools merge -i - > {}".format(bedlist, outname)
    print(cmd)
    call(cmd, shell=True)
    return outname


def launch_AA_AC(jobq, aa_py, PAA_PATH):
    try:
        sname, arg_string = jobq.pop()

    except IndexError:
        return

    with open(sname + "_AA_AC_out.txt", 'w') as outfile:
        time.sleep(random.uniform(0, 0.75))
        cmd = "{} {}{}".format(aa_py, PAA_PATH, arg_string)
        print("\nLaunching AA+AC job for " + sname + "\n" + cmd)
        call(cmd, stdout=outfile, stderr=outfile, shell=True)


def create_AA_AC_cmds(tumor_lines, base_argstring, grouped_seeds):
    cmd_dict = dict()
    for tf in tumor_lines:
        curr_argstring = "{} --run_AA --run_AC -s {} --bam {} --bed {}".format(base_argstring, tf[0], tf[1],
                                                                               grouped_seeds)

        optionals = zip(["--sample_metadata", ], tf[4:])
        for k, v in optionals:
            if v:
                curr_argstring += " {} {}".format(k, v)

        cmd_dict[tf[0]] = curr_argstring

    return cmd_dict


# convert the parsed group input data to PrepareAA commands
def create_CNV_cmds(tumor_lines, normal_lines, base_argstring, cnvkit_dir):
    if not normal_lines:
        normalbam = None

    else:
        normalbam = normal_lines[0]
        if len(normal_lines) > 1:
            print("More than one normal sample specified. Only the first will be used: " + normalbam[0])

    cmd_dict = dict()
    cnv_bed_dict = dict()
    for tf in tumor_lines:
        curr_argstring = "{} -s {} --bam {}".format(base_argstring, tf[0], tf[1])
        if normalbam:
            curr_argstring += " --normal_bam {}".format(normalbam[1])

        optionals = zip(["--cnv_bed", "--sample_metadata"], tf[3:])
        for k, v in optionals:
            if v:
                curr_argstring += " {} {}".format(k, v)
                if k == "--cnv_bed":
                    cnv_bed_dict[tf[0]] = v

        if "--cnv_bed" not in curr_argstring and cnvkit_dir:
            curr_argstring+=" --cnvkit_dir " + cnvkit_dir

        # if QC is desired it will be done during stage 3
        if "--no_QC" not in curr_argstring:
            curr_argstring+=" --no_QC"

        cmd_dict[tf[0]] = curr_argstring

    return cmd_dict, cnv_bed_dict


def make_base_argstring(arg_dict, stop_at_seeds=False):
    base_argstring = ""
    for k, v in arg_dict.items():
        if v is True:
            if k != "no_AA":
                arg = " --" + k
                base_argstring+=arg

        elif v is not False and not k == "input" and not k == "cnvkit_dir":
            arg = " --{} {}".format(k, str(v))
            base_argstring+=arg

    return base_argstring


# read a file providing the group data
def read_group_data(input_file):
    """
    group data is formatted as follows:
    sample_name  bam_file  sample_type
    where 'sample_type' is either 'tumor' or 'normal'
    additional optional fields are as follows:
    cnv_calls  sample_metadata_json
    """
    tumor_lines = []
    normal_lines = []
    with open(input_file) as infile:
        for line in infile:
            if line.startswith("#"):
                continue

            fields = line.rstrip().rsplit()
            if not fields:
                continue

            for ind, v in enumerate(fields):
                if v.upper() == "NA" or v.upper() == "NONE":
                    fields[ind] = None

            if fields[2].lower() == "tumor":
                tumor_lines.append(fields)

            elif fields[2].lower() == "normal":
                normal_lines.append(fields)

            else:
                sys.stderr.write("Input formatting error! Column 3 must either be 'tumor' or 'normal'.\nSee README for "
                                 "group input formatting instructions.\n\n")
                sys.exit(1)


    return tumor_lines, normal_lines


def get_argdict(args):
    arg_dict = dict()
    for arg in vars(args):
        value = getattr(args, arg)
        if value is not None and value != "":
            arg_dict[arg] = value

    return arg_dict


# MAIN #
if __name__ == '__main__':
    # Parses the command line arguments
    parser = argparse.ArgumentParser(
        description="A pipeline wrapper for AmpliconArchitect, invoking alignment CNV calling and CNV filtering prior. "
                    "Can launch AA, as well as downstream amplicon classification.")
    parser.add_argument("-i", "--input", help="Input file providing the multi-sample information. See README for "
                                              "information on how to format the input file.", required=True)
    parser.add_argument("-o", "--output_directory", help="output directory names (will create if not already created)",
                        required=True)
    # parser.add_argument("-s", "--sample_name", help="sample name", required=True)
    parser.add_argument("-t", "--nthreads", help="Number of threads to use in BWA, CNV calling and concurrent "
                                                 "instances of PAA", type=int, required=True)
    parser.add_argument("--no_AA", help="Only produce the union of seeds for the group. Do not run AA/AC",
                        action='store_true')
    # parser.add_argument("--run_AA", help="Run AA after all files prepared. Default off.", action='store_true')
    # parser.add_argument("--run_AC", help="Run AmpliconClassifier after all files prepared. Default off.",
    #                     action='store_true')
    parser.add_argument("--ref", help="Reference genome version.", choices=["hg19", "GRCh37", "GRCh38", "hg38", "mm10",
                                                                            "GRCm38", "GRCh38_viral"])
    parser.add_argument("--cngain", type=float, help="CN gain threshold to consider for AA seeding", default=4.5)
    parser.add_argument("--cnsize_min", type=int, help="CN interval size (in bp) to consider for AA seeding",
                        default=50000)
    parser.add_argument("--downsample", type=float, help="AA downsample argument (see AA documentation)", default=10)
    parser.add_argument("--use_old_samtools", help="Indicate you are using an old build of samtools (prior to version "
                                                   "1.0)", action='store_true', default=False)
    parser.add_argument("--rscript_path", help="Specify custom path to Rscript, if needed when using CNVKit "
                                               "(which requires R version >3.4)")
    parser.add_argument("--python3_path", help="If needed, specify a custom path to python3.")
    parser.add_argument("--aa_python_interpreter",
                        help="By default PrepareAA will use the system's default python path. If you would like to use "
                             "a different python version with AA, set this to either the path to the interpreter or "
                             "'python3' or 'python2'", type=str, default='python')
    # parser.add_argument("--freebayes_dir",
    #                     help="Path to directory where freebayes executable exists (not the path to the executable "
    #                          "itself). Only needed if using Canvas and freebayes is not installed on system path.")
    # parser.add_argument("--vcf", help="VCF (in Canvas format, i.e., \"PASS\" in filter field, AD field as 4th entry of "
    #                     "FORMAT field). When supplied with \"--sorted_bam\", pipeline will start from Canvas CNV stage."
    #                     )
    parser.add_argument("--AA_src", help="Specify a custom $AA_SRC path. Overrides the bash variable")
    parser.add_argument("--AA_runmode", help="If --run_AA selected, set the --runmode argument to AA. Default mode is "
                                             "'FULL'", choices=['FULL', 'BPGRAPH', 'CYCLES', 'SVVIEW'], default='FULL')
    parser.add_argument("--AA_extendmode", help="If --run_AA selected, set the --extendmode argument to AA. Default "
                                                "mode is 'EXPLORE'",
                        choices=["EXPLORE", "CLUSTERED", "UNCLUSTERED", "VIRAL"],
                        default='EXPLORE')
    parser.add_argument("--AA_insert_sdevs", help="Number of standard deviations around the insert size. May need to "
                                                  "increase for sequencing runs with high variance after insert size "
                                                  "selection step. (default 3.0)", type=float, default=3.0)
    # parser.add_argument("--normal_bam", help="Path to matched normal bam for CNVKit (optional)")
    # parser.add_argument("--ploidy", type=float, help="Ploidy estimate for CNVKit (optional). This is not used outside "
    #                                                  "of CNVKit.", default=None)
    # parser.add_argument("--purity", type=float, help="Tumor purity estimate for CNVKit (optional). This is not used "
    #                                                  "outside of CNVKit.", default=None)
    parser.add_argument("--cnvkit_segmentation", help="Segmentation method for CNVKit (if used), defaults to CNVKit "
                                                      "default segmentation method (cbs).",
                        choices=['cbs', 'haar', 'hmm', 'hmm-tumor',
                                 'hmm-germline', 'none'], default='cbs')
    parser.add_argument("--no_filter", help="Do not run amplified_intervals.py to identify amplified seeds",
                        action='store_true')
    parser.add_argument("--no_QC", help="Skip QC on the BAM file.", action='store_true')
    parser.add_argument("--skip_AA_on_normal_bam", help="Skip running AA on the normal bam", action='store_true')
    # parser.add_argument("--sample_metadata", help="Path to a JSON of sample metadata to build on")

    # group = parser.add_mutually_exclusive_group(required=True)
    # group.add_argument("--sorted_bam", "--bam", help="Coordinate sorted BAM file (aligned to an AA-supported "
    #                                                  "reference.)")
    # group.add_argument("--fastqs", help="Fastq files (r1.fq r2.fq)", nargs=2)
    # group.add_argument("--completed_AA_runs",
    #                    help="Path to a directory containing one or more completed AA runs which utilized the same reference genome.")

    # group2 = parser.add_mutually_exclusive_group()
    # group2.add_argument("--cnv_bed", "--bed",
    #                     help="BED file (or CNVKit .cns file) of CNV changes. Fields in the bed file should"
    #                          " be: chr start end name cngain")
    parser.add_argument("--cnvkit_dir", help="Path to cnvkit.py. Assumes CNVKit is on the system path if not set. "
                                             "Not needed if --bed is given.")
    # group2.add_argument("--completed_run_metadata",
    #                     help="Run metadata JSON to retroactively assign to collection of samples", default="")
    # group2.add_argument("--align_only", help="Only perform the alignment stage (do not run CNV calling and seeding",
    #                     action='store_true')

    args = parser.parse_args()

    if args.output_directory and not args.output_directory.endswith('/'):
        args.output_directory+='/'

    if not args.aa_python_interpreter:
        args.aa_python_interpreter = 'python'

    arg_dict = get_argdict(args)
    tumor_lines, normal_lines = read_group_data(args.input)
    print("Found {} tumor samples and {} normals\n".format(str(len(tumor_lines)), str(len(normal_lines))))

    # Stage 1: iterate over and launch each that needs CN calling. collect CN seeds files
    base_argstring = make_base_argstring(arg_dict, stop_at_seeds=True)
    print("Setting base argstring for Stage 1 as:")
    print(base_argstring + "\n")
    cmd_dict, cnv_bed_dict = create_CNV_cmds(tumor_lines, normal_lines, base_argstring, args.cnvkit_dir)
    individual_seed_dct = generate_individual_seeds(cmd_dict, args.aa_python_interpreter, args.output_directory,
                                                    cnv_bed_dict)

    # Stage 2: merge seeds (bedtools - gotta sort and merge), and get new args
    grouped_seeds = group_seeds(individual_seed_dct, args.output_directory)

    # Stage 3: launch each AA job in parallel
    if not args.no_AA:
        if args.skip_AA_on_normal_bam:
            normal_lines = []

        all_lines = normal_lines + tumor_lines
        cmd_dict = create_AA_AC_cmds(all_lines, base_argstring, grouped_seeds)
        threadL = []
        paa_threads = min(args.nthreads, len(all_lines))
        print("\nQueueing " + str(len(all_lines)) + " PAA jobs")
        jobq = []
        for i in range(len(all_lines)):
            sname = all_lines[i][0]
            cmd_string = cmd_dict[sname]
            jobq.append((sname, cmd_string))

        for i in range(paa_threads):
            threadL.append(threading.Thread(target=launch_AA_AC, args=(jobq, args.aa_python_interpreter, PAA_PATH)))
            # threadL.append(workerThread(i, launch_AA_AC, cmd_string, args.aa_python_interpreter, PAA_PATH, sname))
            threadL[i].start()

        for t in threadL:
            t.join()

        print("All jobs completed")
