#!/usr/bin/env python3

# author: Jens Luebeck (jluebeck [at] ucsd.edu)

import argparse
from collections import Counter
import json
import os
import random
import shlex
from subprocess import *
import sys
import time
import threading
import queue

from paalib._version import __ampliconsuitepipeline_version__
from paalib.ac_threading import compute_ac_thread_allocation

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


def launch_AA(job_queue, aa_py, PAA_PATH, parent_odir):
    while not job_queue.empty():
        try:
            sname, cmd_string = job_queue.get(block=False)

            odir = parent_odir + sname
            with open("{}/{}_AA_out.txt".format(odir, sname), 'w') as outfile:
                time.sleep(random.uniform(0, 0.75))
                cmd = "{} {}{}".format(aa_py, PAA_PATH, cmd_string)
                print("\nLaunching AA job for " + sname + "\n" + cmd)
                ecode = call(cmd, stdout=outfile, stderr=outfile, shell=True)
                if ecode != 0:
                    sys.stderr.write("Unexpected error while running AA job!\n")
                    sys.exit(1)

            job_queue.task_done()

        except queue.Empty:
            break


def create_AA_cmds(tumor_lines, base_argstring, grouped_seeds, parent_odir, nthreads=1):
    cmd_dict = dict()
    for tf in tumor_lines:
        odir = parent_odir + tf[0]
        curr_seeds = grouped_seeds[tf[0]]
        curr_argstring = "{} -t {} --run_AA -s {} --bam {} --bed {} -o {}".format(
            base_argstring, nthreads, tf[0], tf[1], curr_seeds, odir)

        optionals = zip(["--sample_metadata", "--sv_vcf"], tf[4:])
        for k, v in optionals:
            if v:
                curr_argstring += " {} {}".format(k, v)

        cmd_dict[tf[0]] = curr_argstring

    return cmd_dict


def _run_logged_command(command, logfile, append=False):
    """Run a command without a shell and send its combined output to a stage log."""
    print("\n" + " ".join(shlex.quote(str(x)) for x in command))
    mode = 'a' if append else 'w'
    with open(logfile, mode) as outfile:
        result = run(command, stdout=outfile, stderr=STDOUT, universal_newlines=True)
    return result.returncode


def _write_sample_path_map(output_path, sample_paths):
    with open(output_path, 'w') as outfile:
        for sample_name, path in sample_paths:
            if path and os.path.exists(path):
                outfile.write("{}\t{}\n".format(sample_name, os.path.realpath(path)))


def _sample_from_feature_path(feature_path, sample_names):
    feature_name = os.path.basename(feature_path)
    for sample_name in sorted(sample_names, key=len, reverse=True):
        if feature_name.startswith(sample_name + "_amplicon"):
            return sample_name
    return None


def update_grouped_metadata(sample_lines, grouped_seeds, class_prefix, ac_command, ac_version):
    """Attach combined-AC metadata to each sample and return metadata map paths."""
    sample_names = [sample[0] for sample in sample_lines]
    amplicon_counts = Counter()
    with open(class_prefix + ".input") as infile:
        for line in infile:
            fields = line.rstrip().rsplit()
            if fields:
                amplicon_counts[fields[0]] += 1

    feature_counts = Counter()
    features_to_graph = class_prefix + "_features_to_graph.txt"
    if os.path.exists(features_to_graph):
        with open(features_to_graph) as infile:
            for line in infile:
                fields = line.rstrip().rsplit()
                if fields:
                    sample_name = _sample_from_feature_path(fields[0], sample_names)
                    if sample_name:
                        feature_counts[sample_name] += 1

    parent_output = os.path.dirname(os.path.dirname(class_prefix))
    sample_metadata_paths = []
    run_metadata_paths = []
    cnv_bed_paths = []
    ac_command_string = " ".join(shlex.quote(str(x)) for x in ac_command)

    for sample_name in sample_names:
        sample_dir = os.path.join(parent_output, sample_name)
        sample_metadata = os.path.join(sample_dir, sample_name + "_sample_metadata.json")
        run_metadata = os.path.join(sample_dir, sample_name + "_run_metadata.json")

        if os.path.exists(sample_metadata):
            with open(sample_metadata) as infile:
                sample_data = json.load(infile)
            sample_data["number_of_AA_amplicons"] = amplicon_counts[sample_name]
            sample_data["number_of_AA_features"] = feature_counts[sample_name]
            with open(sample_metadata, 'w') as outfile:
                json.dump(sample_data, outfile, indent=2)

        if os.path.exists(run_metadata):
            with open(run_metadata) as infile:
                run_data = json.load(infile)
            run_data["AC_cmd"] = ac_command_string
            run_data["AC_version"] = ac_version
            with open(run_metadata, 'w') as outfile:
                json.dump(run_data, outfile, indent=2)

        sample_metadata_paths.append((sample_name, sample_metadata))
        run_metadata_paths.append((sample_name, run_metadata))
        cnv_bed_paths.append((sample_name, grouped_seeds.get(sample_name)))

    metadata_maps = {
        "sample": class_prefix + "_sample_metadata_map.txt",
        "run": class_prefix + "_run_metadata_map.txt",
        "cnv": class_prefix + "_cnv_bed_map.txt",
    }
    _write_sample_path_map(metadata_maps["sample"], sample_metadata_paths)
    _write_sample_path_map(metadata_maps["run"], run_metadata_paths)
    _write_sample_path_map(metadata_maps["cnv"], cnv_bed_paths)
    return metadata_maps


def run_grouped_ac(sample_lines, grouped_seeds, output_directory, input_file, ref, nthreads, python3_path, ac_src):
    """Run one threaded AmpliconClassifier job over every completed AA output."""
    group_name = os.path.splitext(os.path.basename(input_file))[0].replace(" ", "_")
    classification_dir = os.path.join(output_directory, group_name + "_classification")
    if not os.path.exists(classification_dir):
        os.makedirs(classification_dir)

    class_prefix = os.path.join(classification_dir, group_name)
    stage_log = class_prefix + "_AC_stage.log"
    aa_result_dirs = []
    for sample in sample_lines:
        sample_name = sample[0]
        aa_result_dir = os.path.join(output_directory, sample_name, sample_name + "_AA_results")
        if not os.path.isdir(aa_result_dir):
            raise RuntimeError("AA results directory not found for {}: {}".format(sample_name, aa_result_dir))
        aa_result_dirs.append(aa_result_dir)

    make_input_command = [os.path.join(ac_src, "make_input.sh")] + aa_result_dirs + [class_prefix]
    if _run_logged_command(make_input_command, stage_log) != 0:
        raise RuntimeError("Failed to build the combined AmpliconClassifier input")

    input_path = class_prefix + ".input"
    summary_map = class_prefix + "_summary_map.txt"
    if not os.path.exists(input_path) or not os.path.exists(summary_map):
        raise RuntimeError("AmpliconClassifier input or summary map was not created")

    with open(input_path) as infile:
        n_amplicons = sum(1 for line in infile if line.strip())

    ac_script = os.path.join(ac_src, "amplicon_classifier.py")
    help_output = check_output([python3_path, ac_script, "--help"], stderr=STDOUT, universal_newlines=True)
    required_options = ["--jobs", "--bfb_threads", "--no_results_table"]
    missing_options = [option for option in required_options if option not in help_output]
    if missing_options:
        raise RuntimeError("Grouped cohort classification requires AmpliconClassifier 2.0+; missing {}".format(
            ", ".join(missing_options)))

    jobs, bfb_threads = compute_ac_thread_allocation(nthreads, n_amplicons)
    print("Allocating AC threads for {} cohort amplicon(s): --jobs {} --bfb_threads {}".format(
        n_amplicons, jobs, bfb_threads))
    ac_command = [
        python3_path, ac_script,
        "-i", input_path,
        "--ref", ref,
        "-o", class_prefix,
        "--no_results_table",
        "--jobs", str(jobs),
        "--bfb_threads", str(bfb_threads),
    ]
    if _run_logged_command(ac_command, stage_log, append=True) != 0:
        raise RuntimeError("Combined AmpliconClassifier stage failed")

    required_outputs = [
        class_prefix + "_amplicon_classification_profiles.tsv",
        class_prefix + "_features_to_graph.txt",
        class_prefix + "_feature_similarity_scores.tsv",
    ]
    missing_outputs = [path for path in required_outputs if not os.path.exists(path)]
    if missing_outputs:
        raise RuntimeError("Combined AmpliconClassifier stage did not create: {}".format(
            ", ".join(missing_outputs)))

    try:
        ac_version = check_output([python3_path, ac_script, "--version"], stderr=STDOUT,
                                  universal_newlines=True).strip()
    except CalledProcessError:
        ac_version = "NA"

    metadata_maps = update_grouped_metadata(
        sample_lines, grouped_seeds, class_prefix, ac_command, ac_version)
    results_command = [
        python3_path, os.path.join(ac_src, "make_results_table.py"),
        "-i", input_path,
        "--classification_file", class_prefix + "_amplicon_classification_profiles.tsv",
        "--summary_map", summary_map,
        "--sample_metadata_list", metadata_maps["sample"],
        "--run_metadata_list", metadata_maps["run"],
        "--sample_cnv_bed_list", metadata_maps["cnv"],
        "--ref", ref,
    ]
    if _run_logged_command(results_command, stage_log, append=True) != 0:
        raise RuntimeError("Failed to create the combined AmpliconSuite results table")

    if not os.path.exists(class_prefix + "_result_table.tsv"):
        raise RuntimeError("Combined AmpliconSuite results table was not created")

    return class_prefix


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
    parser.add_argument("-t", "--nthreads", help="Total thread budget. CNV preparation runs samples serially; AA "
                                                 "jobs run concurrently; the final cohort-wide AC stage divides "
                                                 "this budget between --jobs and --bfb_threads.", type=int, required=True)
    parser.add_argument("--no_AA", help="Only produce the seeds for the group. Do not run AA/AC",
                        action='store_true')
    parser.add_argument("--no_union", help="Do not create a unified collection of seeds for the group (keep seeds "
                                           "separate between samples", action='store_true')
    parser.add_argument("--ref", help="Reference genome version of all samples.", choices=["hg19", "GRCh37", "GRCh38", "hg38", "mm10",
                        "GRCm38", "GRCh38_viral"], type=str, required=True)
    parser.add_argument("--cngain", type=float, help="CN gain threshold to consider for AA seeding", default=4.5)
    parser.add_argument("--cnsize_min", type=int, help="CN interval size (in bp) to consider for AA seeding",
                        default=50000)
    parser.add_argument("--no_cstats",
                        help="Do not read from or write to the $AA_DATA_REPO/coverage.stats file. (default not set).",
                        action='store_true')
    parser.add_argument("--downsample", type=float, help="AA downsample argument (see AA documentation)", default=10)
    parser.add_argument("--rscript_path", help="Specify custom path to Rscript, if needed when using CNVKit "
                                               "(which requires R version >3.4)")
    parser.add_argument("--python3_path", help="If needed, specify a custom path to python3.")
    parser.add_argument("--aa_python_interpreter",
                        help="By default PrepareAA will use the system's default python path. If you would like to use "
                             "a different python version with AA, set this to either the path to the interpreter or "
                             "'python', 'python3', 'python2' (default 'python3')", type=str, default='python3')
    parser.add_argument("--AA_src", help="Specify a custom $AA_SRC path. Overrides the bash variable")
    parser.add_argument("--AA_runmode", help="Set AA's --runmode during the grouped AA/AC stage. Default mode is "
                                             "'FULL'", choices=['FULL', 'BPGRAPH', 'CYCLES', 'SVVIEW'], default='FULL')
    parser.add_argument("--AA_extendmode", help="Set AA's --extendmode during the grouped AA/AC stage. Default mode "
                                                "is 'EXPLORE'",
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
    parser.add_argument("--AA_solver", help="Set the copy-number optimizer used by AA during the grouped AA/AC stage "
                        "(which runs unless --no_AA is selected). If 'mosek' (default) is unavailable, the pipeline "
                        "selects 'clarabel'; AA also retries with Clarabel if Mosek fails at runtime.",
                        choices=['mosek', 'clarabel'], default='mosek')
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

    # Stage 3: launch each AA-only job in parallel
    if not args.no_AA:
        if args.skip_AA_on_normal_bam:
            normal_lines = []

        for nl in range(len(normal_lines)):
            grouped_seeds[normal_lines[nl][0]] = grouped_seeds[tumor_lines[0][0]]
            odir = "{}{}/".format(args.output_directory, normal_lines[nl][0])
            if not os.path.exists(odir):
                os.makedirs(odir)

        all_lines = normal_lines + tumor_lines
        threadL = []
        # Run as many AA samples concurrently as the thread budget permits and divide the budget across them.
        # AC runs once across all completed AA outputs in Stage 4 and receives the full cohort thread budget.
        paa_threads = min(args.nthreads, len(all_lines))
        aa_threads = max(1, args.nthreads // paa_threads)
        cmd_dict = create_AA_cmds(all_lines, base_argstring, grouped_seeds, args.output_directory, aa_threads)
        print("\nQueueing " + str(len(all_lines)) + " AA jobs")
        print("Going to run {} AA job(s) concurrently, each with -t {}".format(paa_threads, aa_threads))

        # Create a thread-safe queue and fill it directly
        job_queue = queue.Queue()
        for i in range(len(all_lines)):
            sname = all_lines[i][0]
            cmd_string = cmd_dict[sname]
            job_queue.put((sname, cmd_string))

        # Start threads
        for i in range(paa_threads):
            threadL.append(threading.Thread(target=launch_AA,
                                            args=(job_queue, args.aa_python_interpreter, PAA_PATH,
                                                  args.output_directory)))
            threadL[i].start()

        # Wait for all threads to complete
        for t in threadL:
            t.join()

        # Check finish flags for Stage 3
        failed_samples = check_finish_flags(all_lines, args.output_directory)
        if failed_samples:
            sys.stderr.write(f"Stage 3 (AA) failed for samples: {', '.join(failed_samples)}\n")
            sys.exit(1)

        else:
            print("All AA jobs complete\n")

        # Stage 4: classify the complete cohort in one AC process. AC performs feature-similarity scoring internally.
        try:
            class_prefix = run_grouped_ac(
                all_lines, grouped_seeds, args.output_directory, args.input, args.ref,
                args.nthreads, PY3_PATH, AC_SRC)
            print("Combined AC and feature-similarity stages complete\n")
            print("Combined classification prefix: {}\n".format(class_prefix))
        except (OSError, RuntimeError, CalledProcessError) as error:
            sys.stderr.write("Stage 4 (combined AC) failed: {}\n".format(error))
            sys.exit(1)
