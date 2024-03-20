#!/usr/bin/env python3

# author: Jens Luebeck (jluebeck [at] ucsd.edu)

import argparse
import os
import random
from subprocess import *
import sys
import time
import threading

from paalib._version import __ampliconsuitepipeline_version__

PAA_PATH = os.path.dirname(os.path.realpath(__file__)) + "/AmpliconSuite-pipeline.py"
PY3_PATH = "python3"  # updated by command-line arg if specified

try:
    AC_SRC = os.environ['AC_SRC']
except KeyError:
    try:
        import ampclasslib
        ac_path = check_output("which amplicon_classifier.py", shell=True).decode("utf-8")
        AC_SRC = ac_path.rsplit("/amplicon_classifier.py")[0]
    except Exception as e:
        sys.stderr.write(str(e) + "\n")
        sys.stderr.write(
            "\nAC_SRC bash variable or library files not found. AmpliconClassifier may not be properly installed.\n")
        sys.exit(1)

def generate_individual_seeds(cmd_dict, aa_py, parent_odir, cnv_bed_dict):
    individual_seed_dct = {}
    print('Generating individual seeds')
    for sname, argstring in cmd_dict.items():
        odir = "{}{}/".format(parent_odir, sname)
        if not os.path.exists(odir):
            os.makedirs(odir)

        with open("{}{}_CNV_out.txt".format(odir, sname), 'w') as outfile:
            cmd = '{} {}{}'.format(aa_py, PAA_PATH, argstring)
            print(sname)
            print(cmd + "\n")
            call(cmd, stdout=outfile, stderr=outfile, shell=True)

        # if it was a seeds file, PAA won't modify, so copy it into the right location
        if sname in cnv_bed_dict and cnv_bed_dict[sname].endswith("AA_CNV_SEEDS.bed"):
            if not os.path.dirname(os.path.realpath(cnv_bed_dict[sname])) == os.path.realpath(odir):
                cmd = "cp {} {}/".format(cnv_bed_dict[sname], odir)
                call(cmd, shell=True)

            individual_seed_dct[sname] = cnv_bed_dict[sname]

        else:
            # store the name of the path of the seeds file
            individual_seed_dct[sname] = '{}/{}_AA_CNV_SEEDS.bed'.format(odir, sname)

    return individual_seed_dct


def group_seeds(individual_seed_dct, odir):
    samplist = list(individual_seed_dct.keys())
    all_ind_seeds = set(individual_seed_dct.values())
    if len(all_ind_seeds) > 1:
        outname = odir + "_".join(samplist[:2])
        if len(samplist) > 2:
            outname += "_etc_n" + str(len(samplist))

        outname += "_merged_AA_CNV_SEEDS.bed"
        bedlist = " ".join(all_ind_seeds)
        print("Merging seeds")
        cmd = "sort -k1,1 -k2,2n {} | bedtools merge -i - > {}".format(bedlist, outname)
        print(cmd)
        call(cmd, shell=True)

    else:
        outname = all_ind_seeds.pop()

    gs_dict = {x: outname for x in samplist}
    return gs_dict


def launch_AA_AC(jobq, aa_py, PAA_PATH, parent_odir):
    try:
        sname, arg_string = jobq.pop()

    except IndexError:
        return

    odir = parent_odir + sname
    with open("{}/{}_AA_AC_out.txt".format(odir, sname), 'w') as outfile:
        time.sleep(random.uniform(0, 0.75))
        cmd = "{} {}{}".format(aa_py, PAA_PATH, arg_string)
        print("\nLaunching AA+AC job for " + sname + "\n" + cmd)
        call(cmd, stdout=outfile, stderr=outfile, shell=True)


def create_AA_AC_cmds(tumor_lines, base_argstring, grouped_seeds, parent_odir):
    cmd_dict = dict()
    for tf in tumor_lines:
        odir = parent_odir + tf[0]
        curr_seeds = grouped_seeds[tf[0]]
        curr_argstring = "{} --run_AA --run_AC -s {} --bam {} --bed {} -o {}".format(base_argstring, tf[0], tf[1],
                                                                               curr_seeds, odir)

        optionals = zip(["--sample_metadata", ], tf[4:])
        for k, v in optionals:
            if v:
                curr_argstring += " {} {}".format(k, v)

        cmd_dict[tf[0]] = curr_argstring

    return cmd_dict


# convert the parsed group input data to PrepareAA commands
def create_CNV_cmds(tumor_lines, normal_lines, base_argstring, cnvkit_dir, parent_odir):
    if not normal_lines:
        normalbam = None

    else:
        normalbam = normal_lines[0]
        if len(normal_lines) > 1:
            print("More than one normal sample specified. Only the first will be used: " + normalbam[0])

    cmd_dict = dict()
    cnv_bed_dict = dict()
    for tf in tumor_lines:
        odir = parent_odir + tf[0]
        curr_argstring = "{} -s {} --bam {} -o {}".format(base_argstring, tf[0], tf[1], odir)
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
            if k not in ["no_AA", "no_union", "skip_AA_on_normal_bam"]:
                arg = " --" + k
                base_argstring+=arg

        elif k == "AA_insert_sdevs" and v is None:
            continue

        elif v is not False and not k == "input" and not k == "cnvkit_dir" and not k == "output_directory":
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
                if v.upper() == "NA" or v.upper() == "NONE" or v.upper() == "":
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


def concatenate_files(file_paths, output_file):
    try:
        with open(output_file, 'w') as outfile:
            for file_path in file_paths:
                if not os.path.isfile(file_path):
                    sys.stderr.write("Warning: {} does not exist!\n".format(file_path))
                    sys.stderr.write("Feature similarity scoring may not be complete. Please check for AA or AC error.\n")
                else:
                    with open(file_path, 'r') as infile:
                        outfile.write(infile.read())

    except IOError as e:
        print("Error:", e)
        sys.exit(1)


# MAIN #
if __name__ == '__main__':
    # Parses the command line arguments
    parser = argparse.ArgumentParser(
        description="A pipeline wrapper for AmpliconArchitect, invoking alignment CNV calling and CNV filtering prior. "
                    "Can launch AA, as well as downstream amplicon classification on groups of related samples.")
    parser.add_argument("-v", "--version", action='version',
                        version='GroupedAnalysisAmpSuite version {version} \n'.format(version=__ampliconsuitepipeline_version__))
    parser.add_argument("-i", "--input", help="Input file providing the multi-sample information. See README for "
                                              "information on how to format the input file.", required=True)
    parser.add_argument("-o", "--output_directory", help="output directory name (will create if not already created)."
                                                         " Sample outputs will be created as subdirectories inside -o", required=True)
    parser.add_argument("-t", "--nthreads", help="Number of threads to use in BWA, CNV calling and concurrent "
                                                 "instances of PAA", type=int, required=True)
    parser.add_argument("--no_AA", help="Only produce the seeds for the group. Do not run AA/AC",
                        action='store_true')
    parser.add_argument("--no_union", help="Do not create a unified collection of seeds for the group (keep seeds "
                                           "separate between samples", action='store_true')
    parser.add_argument("--ref", help="Reference genome version of all samples.", choices=["hg19", "GRCh37", "GRCh38", "hg38", "mm10",
                        "GRCm38", "GRCh38_viral"], required=True)
    parser.add_argument("--cngain", type=float, help="CN gain threshold to consider for AA seeding", default=4.5)
    parser.add_argument("--cnsize_min", type=int, help="CN interval size (in bp) to consider for AA seeding",
                        default=50000)
    parser.add_argument("--downsample", type=float, help="AA downsample argument (see AA documentation)", default=10)
    parser.add_argument("--rscript_path", help="Specify custom path to Rscript, if needed when using CNVKit "
                                               "(which requires R version >3.4)")
    parser.add_argument("--python3_path", help="If needed, specify a custom path to python3.")
    parser.add_argument("--aa_python_interpreter",
                        help="By default PrepareAA will use the system's default python path. If you would like to use "
                             "a different python version with AA, set this to either the path to the interpreter or "
                             "'python3' or 'python2' (default 'python')", type=str, default='python')
    parser.add_argument("--AA_src", help="Specify a custom $AA_SRC path. Overrides the bash variable")
    parser.add_argument("--AA_runmode", help="If --run_AA selected, set the --runmode argument to AA. Default mode is "
                                             "'FULL'", choices=['FULL', 'BPGRAPH', 'CYCLES', 'SVVIEW'], default='FULL')
    parser.add_argument("--AA_extendmode", help="If --run_AA selected, set the --extendmode argument to AA. Default "
                                                "mode is 'EXPLORE'",
                        choices=["EXPLORE", "CLUSTERED", "UNCLUSTERED", "VIRAL"],
                        default='EXPLORE')
    parser.add_argument("--AA_insert_sdevs", help="Number of standard deviations around the insert size. May need to "
                                                  "increase for sequencing runs with high variance after insert size "
                                                  "selection step. (default 3.0)", type=float, default=None)
    parser.add_argument("--cnvkit_segmentation", help="Segmentation method for CNVKit (if used), defaults to CNVKit "
                                                      "default segmentation method (cbs).",
                        choices=['cbs', 'haar', 'hmm', 'hmm-tumor',
                                 'hmm-germline', 'none'], default='cbs')
    parser.add_argument("--no_filter", help="Do not run amplified_intervals.py to remove low confidence candidate seed"
                                            " regions overlapping repetitive parts of the genome", action='store_true')
    parser.add_argument("--no_QC", help="Skip QC on the BAM file.", action='store_true')
    parser.add_argument("--skip_AA_on_normal_bam", help="Skip running AA on the normal bam", action='store_true')
    parser.add_argument("--cnvkit_dir", help="Path to cnvkit.py. Assumes CNVKit is on the system path if not set. "
                                             "Not needed if --bed is given.")

    args = parser.parse_args()

    if args.output_directory and not args.output_directory.endswith('/'):
        args.output_directory += '/'

    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)

    if not args.aa_python_interpreter:
        args.aa_python_interpreter = 'python'

    if args.python3_path:
        if not args.python3_path.endswith("/python") and not args.python3_path.endswith("/python3"):
            args.python3_path += "/python3"

        PY3_PATH = args.python3_path

    arg_dict = get_argdict(args)
    tumor_lines, normal_lines = read_group_data(args.input)
    print("Found {} tumor samples and {} normals\n".format(str(len(tumor_lines)), str(len(normal_lines))))
    if len(tumor_lines) == 0:
        print("No tumor samples were provided. Exiting.")
        sys.exit(1)

    # Stage 1: iterate over and launch each that needs CN calling. collect CN seeds files
    base_argstring = make_base_argstring(arg_dict, stop_at_seeds=True)
    print("Setting base argstring for Stage 1 as:")
    print(base_argstring + "\n")
    cmd_dict, cnv_bed_dict = create_CNV_cmds(tumor_lines, normal_lines, base_argstring, args.cnvkit_dir,
                                             args.output_directory)
    individual_seed_dct = generate_individual_seeds(cmd_dict, args.aa_python_interpreter, args.output_directory,
                                                    cnv_bed_dict)

    # Stage 2: merge seeds (bedtools - gotta sort and merge)
    if args.no_union:
        grouped_seeds = individual_seed_dct
    else:
        grouped_seeds = group_seeds(individual_seed_dct, args.output_directory)

    # Stage 3: launch each AA job in parallel
    if not args.no_AA:
        if args.skip_AA_on_normal_bam:
            normal_lines = []
        elif normal_lines:
            grouped_seeds[normal_lines[0][0]] = grouped_seeds[tumor_lines[0][0]]
            odir = "{}{}/".format(args.output_directory, normal_lines[0][0])
            if not os.path.exists(odir):
                os.makedirs(odir)

        all_lines = normal_lines + tumor_lines
        cmd_dict = create_AA_AC_cmds(all_lines, base_argstring, grouped_seeds, args.output_directory)
        threadL = []
        paa_threads = min(args.nthreads, len(all_lines))
        print("\nQueueing " + str(len(all_lines)) + " PAA jobs")
        jobq = []
        for i in range(len(all_lines)):
            sname = all_lines[i][0]
            cmd_string = cmd_dict[sname]
            jobq.append((sname, cmd_string))

        for i in range(paa_threads):
            threadL.append(threading.Thread(target=launch_AA_AC, args=(jobq, args.aa_python_interpreter, PAA_PATH,
                                                                       args.output_directory)))
            threadL[i].start()

        for t in threadL:
            t.join()

        print("All AA & AC jobs completed")

        # Stage 4: run feature similarity on outputs.
        # make feature_input file (concatenate from each job)
        feat_files = []
        for i in range(len(all_lines)):
            sname = all_lines[i][0]
            feat_graph_file = ("{}{}/{}_classification/{}_features_to_graph.txt"
                               .format(args.output_directory, sname, sname, sname))
            feat_files.append(feat_graph_file)

        combined_feat_graph_file = args.output_directory + "combined_features_to_graph.txt"
        concatenate_files(feat_files, combined_feat_graph_file)
        cmd = ("{} {}/feature_similarity.py -f {} --ref {} -o {}combined_samples"
               .format(PY3_PATH, AC_SRC, combined_feat_graph_file, args.ref, args.output_directory, ))
        print(cmd)
        call(cmd, shell=True)

        print("Feature similarity calculations completed\n")
