#!/usr/bin/env python3

# author: Jens Luebeck (jluebeck [at] ucsd.edu)

import argparse
import os
import random
from subprocess import *
import sys
import time
import threading
import queue

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
            ecode = call(cmd, stdout=outfile, stderr=outfile, shell=True)
            if ecode != 0:
                sys.stderr.write("Unexpected error during individual CNV call generation!\n")
                sys.exit(1)

        # if it was a seeds file, PAA won't modify, so copy it into the right location
        if sname in cnv_bed_dict and cnv_bed_dict[sname].endswith("AA_CNV_SEEDS.bed"):
            if not os.path.dirname(os.path.realpath(cnv_bed_dict[sname])) == os.path.realpath(odir):
                cmd = "cp {} {}/".format(cnv_bed_dict[sname], odir)
                ecode = call(cmd, shell=True)
                if ecode != 0:
                    sys.stderr.write("Unexpected error while copying AA_CNV_SEEDS.bed file!\n")
                    sys.exit(1)

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
        ecode = call(cmd, shell=True)
        if ecode != 0:
            sys.stderr.write("Unexpected error while merging seeds!\n")
            sys.exit(1)

    else:
        outname = all_ind_seeds.pop()

    gs_dict = {x: outname for x in samplist}
    return gs_dict


def launch_AA_AC(job_queue, aa_py, PAA_PATH, parent_odir):
    while not job_queue.empty():
        try:
            sname, cmd_string = job_queue.get(block=False)

            odir = parent_odir + sname
            with open("{}/{}_AA_AC_out.txt".format(odir, sname), 'w') as outfile:
                time.sleep(random.uniform(0, 0.75))
                cmd = "{} {}{}".format(aa_py, PAA_PATH, cmd_string)
                print("\nLaunching AA+AC job for " + sname + "\n" + cmd)
                ecode = call(cmd, stdout=outfile, stderr=outfile, shell=True)
                if ecode != 0:
                    sys.stderr.write("Unexpected error while running AA+AC job!\n")
                    sys.exit(1)

            job_queue.task_done()

        except queue.Empty:
            break


def create_AA_AC_cmds(tumor_lines, base_argstring, grouped_seeds, parent_odir):
    cmd_dict = dict()
    for tf in tumor_lines:
        odir = parent_odir + tf[0]
        curr_seeds = grouped_seeds[tf[0]]
        curr_argstring = "{} -t 1 --run_AA --run_AC -s {} --bam {} --bed {} -o {}".format(base_argstring, tf[0], tf[1],
                                                                               curr_seeds, odir)

        optionals = zip(["--sample_metadata", "--sv_vcf"], tf[4:])
        for k, v in optionals:
            if v:
                curr_argstring += " {} {}".format(k, v)

        cmd_dict[tf[0]] = curr_argstring

    return cmd_dict


# convert the parsed group input data to PrepareAA commands
def create_CNV_cmds(tumor_lines, normal_lines, base_argstring, cnvkit_dir, parent_odir, nthreads):
    if not normal_lines:
        normalbam = None

    else:
        normalbam = normal_lines[0]
        if len(normal_lines) > 1:
            print("More than one normal sample specified. Only the first will be used for matched tumor-normal CNV calling: " + normalbam[0])

    cmd_dict = dict()
    cnv_bed_dict = dict()
    for tf in tumor_lines:
        odir = parent_odir + tf[0]
        curr_argstring = "{} -t {} -s {} --bam {} -o {}".format(base_argstring, nthreads, tf[0], tf[1], odir)
        if normalbam:
            curr_argstring += " --normal_bam {}".format(normalbam[1])

        optionals = zip(["--cnv_bed", "--sample_metadata", "--sv_vcf"], tf[3:])
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
        if k == "nthreads":
            continue

        if v is True:
            if k not in ["no_AA", "no_union", "skip_AA_on_normal_bam"]:
                arg = " --" + k
                base_argstring+=arg

        elif any([k == x for x in ["AA_insert_sdevs", "pair_support_min", "foldback_pair_support_min"]]) and v is None:
            continue

        elif v is not False and not k == "input" and not k == "cnvkit_dir" and not k == "output_directory":
            arg = " --{} {}".format(k, str(v))
            base_argstring+=arg

    return base_argstring


def verify_file_exists(file_path):
    """
    Verify if a file exists
    Args:
        file_path: Path to the file

    Returns:
        True if file exists, False otherwise
    """
    if file_path is None:
        return True

    # Check if file exists
    if os.path.isfile(file_path):
        return True
    else:
        return False


# Modify the read_group_data function to also check file existence
def read_group_data(input_file):
    """
    group data is formatted as follows:
    sample_name  bam_file  sample_type
    where 'sample_type' is either 'tumor' or 'normal'
    additional optional fields are as follows:
    cnv_calls  sample_metadata_json  sv_calls
    """
    data_len = 6
    tumor_lines = []
    normal_lines = []
    seen_names = set()
    missing_files = []

    with open(input_file) as infile:
        for line in infile:
            if line.startswith("#"):
                continue

            fields = line.rstrip().rsplit()
            if not fields:
                continue
            elif len(fields) < 3:
                sys.stderr.write("Input formatting error on line below! Too few fields.\n")
                sys.stderr.write(line + "\n")
                sys.stderr.write("See README for group input formatting instructions.\n")
                sys.exit(1)

            for ind, v in enumerate(fields):
                if v.upper() == "NA" or v.upper() == "NONE" or v.upper() == "":
                    fields[ind] = None

            if len(fields) < data_len:
                fields.extend([None] * (data_len - len(fields)))

            if fields[2].lower() == "tumor":
                tumor_lines.append(fields)

            elif fields[2].lower() == "normal":
                normal_lines.append(fields)

            else:
                sys.stderr.write("Input formatting error! Column 3 must either be 'tumor' or 'normal'.\nSee README for "
                                 "group input formatting instructions.\n\n")
                sys.exit(1)

            if fields[3] and not any(fields[3].endswith(x) for x in ['.bed', '.cns']):
                sys.stderr.write(
                    "Input formatting error! Column 4 (CNV calls) must either be 'NA' or a .bed or .cns file.\nSee README for "
                    "group input formatting instructions.\n\n")
                sys.exit(1)

            elif fields[4] and not fields[4].endswith('.json'):
                sys.stderr.write(
                    "Input formatting error! Column 5 (Sample metadata json) must either be 'NA' or a .json file.\nSee README for "
                    "group input formatting instructions.\n\n")
                sys.exit(1)

            elif fields[5] and not fields[5].endswith('.vcf'):
                sys.stderr.write(
                    "Input formatting error! Column 6 (external SV calls) must either be 'NA' or a .vcf file.\nSee README for "
                    "group input formatting instructions.\n\n")
                sys.exit(1)

            if fields[0] in seen_names:
                sys.stderr.write(
                    "Duplicate sample name {} in .input file! Sample names must be unique.\n".format(fields[0]))
                sys.exit(1)

            seen_names.add(fields[0])

            # Check if required BAM file exists
            if not verify_file_exists(fields[1]):
                missing_files.append("BAM file '{}' for sample '{}'".format(fields[1], fields[0]))

            # Check if optional files that were provided exist
            file_types = ["CNV file", "metadata JSON file", "SV VCF file"]
            for i, file_path in enumerate(fields[3:6]):
                if file_path and not verify_file_exists(file_path):
                    missing_files.append("{} '{}' for sample '{}'".format(file_types[i], file_path, fields[0]))

    # If any files are missing, exit with error
    if missing_files:
        sys.stderr.write("\nERROR: The following input files do not exist:\n")
        for file_desc in missing_files:
            sys.stderr.write("  - {}\n".format(file_desc))
        sys.stderr.write("\nPlease ensure all input files exist before running the pipeline.\n")
        sys.exit(1)

    if tumor_lines:
        print("All input files verified successfully.")

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


def check_finish_flags(sample_list, output_directory):
    """Check finish flag files for UNSUCCESSFUL status"""
    failed_samples = []

    for sample_info in sample_list:
        sample_name = sample_info[0]
        finish_flag_path = os.path.join(output_directory, sample_name, f"{sample_name}_finish_flag.txt")

        if os.path.exists(finish_flag_path):
            with open(finish_flag_path, 'r') as f:
                first_line = f.readline().strip()
                if first_line == "UNSUCCESSFUL":
                    failed_samples.append(sample_name)
        else:
            # No finish flag file means job didn't complete
            failed_samples.append(sample_name)

    return failed_samples

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
    parser.add_argument("-t", "--nthreads", help="Number of threads to use in BWA, CNV calling steps (samples run in serial), "
                                                 "and the number of concurrent instances of AmpliconSuite to launch at once", type=int, required=True)
    parser.add_argument("--no_AA", help="Only produce the seeds for the group. Do not run AA/AC",
                        action='store_true')
    parser.add_argument("--no_union", help="Do not create a unified collection of seeds for the group (keep seeds "
                                           "separate between samples", action='store_true')
    parser.add_argument("--ref", help="Reference genome version of all samples.", choices=["hg19", "GRCh37", "GRCh38", "hg38", "mm10",
                        "GRCm38", "GRCh38_viral"], type=str, required=True)
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
                             "'python', 'python3', 'python2' (default 'python3')", type=str, default='python3')
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
    parser.add_argument('--pair_support_min', dest='pair_support_min', help="Number of read pairs for "
                        "minimum breakpoint support (default 2 but typically becomes higher due to coverage-scaled "
                        "cutoffs)", metavar='INT', action='store', type=int)
    parser.add_argument('--foldback_pair_support_min', help="Number of read pairs for minimum foldback SV support "
                        "(default 2 but typically becomes higher due to coverage-scaled cutoffs). Used value will be the maximum"
                        " of pair_support and this argument. Raising to 3 will help dramatically in heavily artifacted samples.",
                        metavar='INT', action='store', type=int)
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
    parser.add_argument("--sv_vcf_no_filter", help="Use all external SV calls from the --sv_vcf arg, even "
                        "those without 'PASS' in the FILTER column.", action='store_true', default=False)
    parser.add_argument('--sv_vcf_include_sr',
                        help="Include single-ended reads when counting support for an SV in the provided VCF",
                        action='store_true', default=False)

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
    print("\nRead {} tumor samples and {} normal samples\n".format(str(len(tumor_lines)), str(len(normal_lines))))
    if len(tumor_lines) == 0:
        print("No tumor samples were provided. Exiting.")
        sys.exit(1)

    # Stage 1: iterate over and launch each that needs CN calling. collect CN seeds files
    base_argstring = make_base_argstring(arg_dict, stop_at_seeds=True)
    print("Setting base argstring for Stage 1 as:")
    print(base_argstring + "\n")
    cmd_dict, cnv_bed_dict = create_CNV_cmds(tumor_lines, normal_lines, base_argstring, args.cnvkit_dir,
                                             args.output_directory, args.nthreads)
    individual_seed_dct = generate_individual_seeds(cmd_dict, args.aa_python_interpreter, args.output_directory,
                                                    cnv_bed_dict)

    failed_samples = check_finish_flags(tumor_lines, args.output_directory)
    if failed_samples:
        sys.stderr.write(f"CNV calling failed for samples: {', '.join(failed_samples)}\n")
        sys.exit(1)

    else:
        print("CNV calling stage complete\n")

    # Stage 2: merge seeds (bedtools - gotta sort and merge)
    if args.no_union:
        grouped_seeds = individual_seed_dct
    else:
        grouped_seeds = group_seeds(individual_seed_dct, args.output_directory)

    print("CNV merging stage complete\n")

    # Stage 3: launch each AA job in parallel
    if not args.no_AA:
        if args.skip_AA_on_normal_bam:
            normal_lines = []

        for nl in range(len(normal_lines)):
            grouped_seeds[normal_lines[nl][0]] = grouped_seeds[tumor_lines[0][0]]
            odir = "{}{}/".format(args.output_directory, normal_lines[nl][0])
            if not os.path.exists(odir):
                os.makedirs(odir)

        all_lines = normal_lines + tumor_lines
        cmd_dict = create_AA_AC_cmds(all_lines, base_argstring, grouped_seeds, args.output_directory)
        threadL = []
        paa_threads = min(args.nthreads, len(all_lines))
        print("\nQueueing " + str(len(all_lines)) + " PAA jobs")
        print("Going to use " + str(paa_threads) + " threads")

        # Create a thread-safe queue and fill it directly
        job_queue = queue.Queue()
        for i in range(len(all_lines)):
            sname = all_lines[i][0]
            cmd_string = cmd_dict[sname]
            job_queue.put((sname, cmd_string))

        # Start threads
        for i in range(paa_threads):
            threadL.append(threading.Thread(target=launch_AA_AC,
                                            args=(job_queue, args.aa_python_interpreter, PAA_PATH,
                                                  args.output_directory)))
            threadL[i].start()

        # Wait for all threads to complete
        for t in threadL:
            t.join()

        # Check finish flags for Stage 3
        failed_samples = check_finish_flags(all_lines, args.output_directory)
        if failed_samples:
            sys.stderr.write(f"Stage 3 (AA/AC) failed for samples: {', '.join(failed_samples)}\n")
            sys.exit(1)

        else:
            print("All AA & AC jobs complete\n")

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
        ecode = call(cmd, shell=True)
        if ecode == 0:
            print("Feature similarity calculations completed\n")
        else:
            sys.stderr.write("Unexpected error when computing feature similarity!\n")
            sys.exit(1)
