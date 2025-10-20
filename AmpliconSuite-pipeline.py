#!/usr/bin/env python3

# author: Jens Luebeck (jluebeck [at] ucsd.edu)

from datetime import datetime
import json
import logging
import os
import shutil
import socket
from subprocess import *
import sys
import tarfile
import time
import zipfile

from paalib import check_reference, reduce_fasta, cnv_plots
from paalib.argument_parser import setup_argument_parser
from paalib.config_validator import (
    validate_arguments, setup_environment_and_paths, setup_tool_paths,
    validate_aa_environment, initialize_logging_and_directories,
    create_coverage_stats_file, get_samtools_version
)
from paalib.repo_downloader import handle_repo_download
from paalib.run_uploader import archive_and_upload_sample
from paalib._version import __ampliconsuitepipeline_version__

# Global variables
# These will be set in main() and available to all functions
args = type('Args', (), {})()  # Empty object that won't cause attribute errors
AA_REPO = ""
AA_SRC = ""
AC_SRC = ""
metadata_dict = {}
sample_info_dict = {}  # stores the sample metadata
ref_genome_size_file = ""
PY3_PATH = "python3"


def run_bwa(ref_fasta, fastqs, outdir, sname, nthreads, samtools, samtools_version):
    outname = outdir + sname
    logging.info("Output prefix: " + outname)

    exts = [".sa", ".amb", ".ann", ".pac", ".bwt"]
    indexPresent = True
    for i in exts:
        if not os.path.exists(ref_fasta + i):
            indexPresent = False
            logging.info(
                "Could not find " + ref_fasta + i + ", building BWA index from scratch. This could take >60 minutes")
            break

    if not indexPresent:
        cmd = "bwa index " + ref_fasta
        call(cmd, shell=True)

    logging.info("Performing alignment and sorting\n")
    sort_threads = min(int(nthreads), 4)

    final_bam_name = "{}.cs.rmdup.bam".format(outname)
    markdup_stats = "{}_markdup_stats.txt".format(outname)
    stderr_log = "{}_aln_stage.stderr".format(outname)
    logging.info("samtools version {}.{}".format(str(samtools_version[0]), str(samtools_version[1])))

    # --- Samtools <1.0 -------------------------------------------------------
    if samtools_version[0] < 1:
        logging.warning("Samtools version <1.0 — using rmdup pipeline.")
        cmd = "{{ bwa mem -K 10000000 -t {} {} {} | {} view -Shu - | {} sort -m 4G -@{} - {}.cs; }} 2>{}".format(
            nthreads, ref_fasta, fastqs, samtools, samtools, sort_threads, outname, stderr_log)

    # --- Samtools 1.0–1.5 ----------------------------------------------------
    elif samtools_version[0] == 1 and samtools_version[1] < 6:
        logging.warning("Samtools version 1.0–1.5 — markdup not available; using rmdup pipeline.")
        cmd = "{{ bwa mem -K 10000000 -t {} {} {} | {} view -bS - | {} sort -m 4G -@{} -o {}.cs.bam -; }} 2>{}".format(
            nthreads, ref_fasta, fastqs, samtools, samtools, sort_threads, outname, stderr_log)

    # --- Samtools 1.6–1.10 ---------------------------------------------------
    elif samtools_version[0] == 1 and samtools_version[1] < 11:
        logging.info("Samtools 1.6–1.10 — using markdup with -O bam (no -u support).")
        cmd = "{{ bwa mem -K 10000000 -t {} {} {} | {} fixmate -m -O bam - - | {} sort -m 4G -@{} -O bam - | {} markdup -r -s -O bam - {} 2>{}; }} 2>{}".format(
            nthreads, ref_fasta, fastqs, samtools, samtools, sort_threads, samtools,
            final_bam_name, markdup_stats, stderr_log)

    # --- Samtools >=1.11 -----------------------------------------------------
    else:
        cmd = "{{ bwa mem -K 10000000 -t {} {} {} | {} fixmate -m -u - - | {} sort -m 4G -@{} -u - | {} markdup -r -s - {} 2>{}; }} 2>{}".format(
            nthreads, ref_fasta, fastqs, samtools, samtools, sort_threads, samtools,
            final_bam_name, markdup_stats, stderr_log)

    logging.info(cmd + "\n")
    call(cmd, shell=True)
    metadata_dict["bwa_cmd"] = cmd

    # Perform duplicate removal & indexing for old samtools
    if samtools_version[0] < 1 or (samtools_version[0] == 1 and samtools_version[1] < 6):
        logging.info("Performing duplicate removal with rmdup")
        cmd_list = [samtools, "rmdup", "-s", "{}.cs.bam".format(outname), final_bam_name]
        logging.info(" ".join(cmd_list) + "\n")
        call(cmd_list)
        logging.info("Removing temp BAM\n")
        cmd = "rm {}.cs.bam".format(outname)
        call(cmd, shell=True)

    logging.info("Running samtools index")
    cmd_list = [samtools, "index", final_bam_name]
    logging.info(" ".join(cmd_list) + "\n")
    call(cmd_list)

    return final_bam_name, stderr_log


def run_cnvkit(ckpy_path, nthreads, outdir, bamfile, seg_meth='cbs', normal=None, ref_fasta=None, vcf=None):
    # CNVkit cmd-line args
    # -m wgs: wgs data
    # -y: assume chrY present
    # -n: create flat reference (cnv baseline)
    # -p: number of threads
    # -f: reference genome fasta
    bamBase = os.path.splitext(os.path.basename(bamfile))[0]

    env = os.environ.copy()
    env['NUMEXPR_MAX_THREADS'] = env.get('NUMEXPR_MAX_THREADS', str(nthreads))

    cnvkit_version = Popen([PY3_PATH, ckpy_path, "version"], stdout=PIPE, stderr=PIPE, universal_newlines=True, env=env).communicate()[0].rstrip()

    metadata_dict["cnvkit_version"] = cnvkit_version

    ckRef = AA_REPO + args.ref + "/" + args.ref + "_cnvkit_filtered_ref.cnn"
    if normal and args.ref == "GRCh38_viral":
        logging.warning("CNVkit does not properly support matched tumor-normal with viral genomes. Ignoring matched-"
                        "normal and running in tumor-only mode.\n")
        
    logging.info("Running CNVKit batch\n")

    rscript_str = ""
    if args.rscript_path:
        rscript_str = " --rscript-path " + args.rscript_path
        logging.info("Set Rscript flag: " + rscript_str)

    if normal and not args.ref == "GRCh38_viral":
        # create a version of the stripped reference
        reduce_fasta.reduce_fasta(ref_fasta, ref_genome_size_file, outdir)
        base = os.path.basename(ref_fasta)  # args.ref is the name, ref is the fasta
        stripRefG = outdir + os.path.splitext(base)[0] + "_reduced" + "".join(os.path.splitext(base)[1:])
        logging.info("Stripped reference: " + stripRefG)
        cmd = "{} {} batch {} -m wgs{} --fasta {} -p {} -d {} --normal {}".format(PY3_PATH, ckpy_path, bamfile,
                                                                        rscript_str, stripRefG, nthreads, outdir, normal)
    else:
        cmd = "{} {} batch -m wgs{} -r {} -p {} -d {} {}".format(PY3_PATH, ckpy_path, rscript_str, ckRef, nthreads, outdir, bamfile)

    logging.info(cmd + "\n")
    call(cmd, shell=True, env=env)
    metadata_dict["cnvkit_cmd"] = cmd + " ; "

    cnrFile = outdir + bamBase + ".cnr"
    cnsFile = outdir + bamBase + ".cns"

    logging.info("Running CNVKit segment")
    # TODO: Allow a .cnr as input for --bed arg and then jump directly to this step?
    cmd = "{} {} segment {}{} -p {} -m {} -o {}".format(PY3_PATH, ckpy_path, cnrFile, rscript_str, nthreads, seg_meth,
                                                         cnsFile)
    logging.info(cmd + "\n")

    # Use Popen to capture stderr
    process = Popen(cmd, shell=True, stderr=PIPE, universal_newlines=True, env=env)
    stdout, stderr = process.communicate()
    print(stdout)

    if process.returncode != 0:
        logging.error("CNVKit encountered a non-zero exit status ({}). Error message:\n{}".format(
            process.returncode, stderr))
        sys.exit(1)

    metadata_dict["cnvkit_cmd"] = metadata_dict["cnvkit_cmd"] + cmd
    logging.info("Cleaning up temporary CNVkit files")
    cmd = "rm -f {}/*tmp.bed {}/*.cnn {}/*target.bed {}/*.bintest.cns".format(outdir, outdir, outdir, outdir)
    logging.info(cmd)
    call(cmd, shell=True)
    cmd = "gzip -f " + cnrFile
    logging.info(cmd)
    call(cmd, shell=True)
    if normal and not args.ref == "GRCh38_viral":
        cmd = "rm " + stripRefG + " " + stripRefG + ".fa" + " " + stripRefG + ".bed"
        logging.info(cmd)
        call(cmd, shell=True)


# Read the CNVkit .cns files
def convert_cnvkit_cns_to_bed(cnvkit_output_directory, base, cnsfile=None, rescaled=False, nofilter=False):
    if cnsfile is None:
        if not rescaled:
            cnsfile = cnvkit_output_directory + base + ".cns"
        else:
            cnsfile = cnvkit_output_directory + base + "_rescaled.cns"

    with open(cnsfile) as infile, open(cnvkit_output_directory + base + "_CNV_CALLS.bed", 'w') as outfile:
        head = next(infile).rstrip().rsplit("\t")
        for line in infile:
            fields = line.rstrip().rsplit("\t")
            # s, e = int(fields[1]), int(fields[2])
            cn_r = float(fields[4])
            cn = 2 ** (cn_r + 1)
            # do not filter on size since amplified_intervals.py will merge small ones.
            outline = "\t".join(fields[0:3] + ["CNVkit", str(cn)]) + "\n"
            outfile.write(outline)

    return cnvkit_output_directory + base + "_CNV_CALLS.bed"


def rescale_cnvkit_calls(ckpy_path, cnvkit_output_directory, base, cnsfile=None, ploidy=None, purity=None):
    if purity is None and ploidy is None:
        logging.warning("Warning: Rescaling called without --ploidy or --purity. Rescaling will have no effect.")
    if cnsfile is None:
        cnsfile = cnvkit_output_directory + base + ".cns"

    if purity < 0.4:
        logging.warning("WARNING! Rescaling a low purity sample may cause many false-positive seed regions!")
        
    cmd = "{} {} call {} -m clonal".format(PY3_PATH, ckpy_path, cnsfile)
    if purity:
        cmd += " --purity " + str(purity)
    if ploidy:
        cmd += " --ploidy " + str(ploidy)

    cmd += " -o " + cnvkit_output_directory + base + "_rescaled.cns"
    logging.info("Rescaling CNVKit calls\n" + cmd)
    call(cmd, shell=True)


def run_amplified_intervals(AA_interpreter, CNV_seeds_filename, sorted_bam, output_directory, sname, cngain,
                            cnsize_min, paa_logfile):
    logging.info("Running amplified_intervals")
    AA_seeds_filename = "{}_AA_CNV_SEEDS".format(output_directory + sname)

    # Test if --logfile is supported by checking help output
    help_cmd = "{} {}/amplified_intervals.py --help".format(AA_interpreter, AA_SRC)
    try:
        help_output = check_output(help_cmd, shell=True, stderr=STDOUT, universal_newlines=True)
        supports_logfile = "--logfile" in help_output
    except CalledProcessError:
        supports_logfile = False
        logging.warning(
            "Could not check amplified_intervals.py help. Assuming older version without --logfile support.")

    # Build base command
    cmd = "{} {}/amplified_intervals.py --ref {} --bed {} --bam {} --gain {} --cnsize_min {} --out {}".format(
        AA_interpreter, AA_SRC, args.ref, CNV_seeds_filename, sorted_bam, str(cngain), str(cnsize_min),
        AA_seeds_filename)

    # Add --logfile if supported
    if supports_logfile:
        cmd += " --logfile {}".format(paa_logfile)
        logging.info("Using amplified_intervals.py with --logfile support")
    else:
        logging.info("Using amplified_intervals.py without --logfile (older version)")

    logging.info(cmd + "\n")
    exit_code = call(cmd, shell=True)
    if exit_code != 0:
        logging.error("amplified_intervals.py returned a non-zero exit code. Exiting...\n")
        sys.exit(1)

    metadata_dict["amplified_intervals_cmd"] = cmd
    return AA_seeds_filename + ".bed"


def run_AA(amplified_interval_bed, AA_outdir, sname, args):
    logging.info("Running AA")
    AA_interpreter = args.aa_python_interpreter
    sorted_bam = args.bam
    downsample = args.downsample
    ref = args.ref
    runmode = args.AA_runmode
    extendmode = args.AA_extendmode
    insert_sdevs = args.AA_insert_sdevs
    sv_vcf = args.sv_vcf
    sv_vcf_no_filter = args.sv_vcf_no_filter
    pair_support = args.pair_support_min
    fb_pair_support = args.foldback_pair_support_min

    AA_version = \
    Popen([AA_interpreter, AA_SRC + "/AmpliconArchitect.py", "--version"], stdout=PIPE, stderr=PIPE, universal_newlines=True).communicate()[0].rstrip()
    if not AA_version:
        AA_version = \
            Popen([AA_interpreter, AA_SRC + "/AmpliconArchitect.py", "--version"], stdout=PIPE, stderr=PIPE, universal_newlines=True).communicate()[1].rstrip()

    # try:
    #     AA_version = AA_version.decode('utf-8')
    # except UnicodeError:
    #     pass

    metadata_dict["AA_version"] = AA_version

    cmd = "{} {}/AmpliconArchitect.py --ref {} --downsample {} --bed {} --bam {} --runmode {} --extendmode {} --out {}/{}".format(
        AA_interpreter, AA_SRC, ref, str(downsample), amplified_interval_bed, sorted_bam, runmode, extendmode,
        AA_outdir, sname)
    if insert_sdevs is not None:
        cmd += " --insert_sdevs {}".format(str(insert_sdevs))

    if sv_vcf:
        cmd += " --sv_vcf {}".format(sv_vcf)
        if sv_vcf_no_filter:
            cmd += " --sv_vcf_no_filter"

    if pair_support:
        cmd += " --pair_support_min {}".format(str(pair_support))

    if fb_pair_support:
        cmd += " --foldback_pair_support_min {}".format(str(fb_pair_support))

    logging.info(cmd + "\n")
    aa_exit_code = call(cmd, shell=True)
    if aa_exit_code != 0:
        logging.error("AmpliconArchitect returned a non-zero exit code. Exiting...\n")
        sys.exit(1)

    metadata_dict["AA_cmd"] = cmd


def run_AC(AA_outdir, sname, ref, AC_outdir, AC_src):
    logging.info("")
    logging.info("Running AC")
    # make input file
    class_output = AC_outdir + sname
    input_file = class_output + ".input"
    bed_dir = class_output + "_classification_bed_files/"
    if os.path.exists(bed_dir):
        logging.warning("WARNING! AC files were not cleared prior to re-running. New classifications may become "
                        "mixed with previous classification files!")

    cmd = "{}/make_input.sh {} {}".format(AC_src, AA_outdir, class_output)
    logging.info(cmd)

    # Capture stdout and stderr
    result = run(cmd, shell=True, stdout=PIPE, stderr=PIPE, universal_newlines=True)

    # Log stdout as info if there's any output
    if result.stdout.strip():
        for line in result.stdout.strip().split('\n'):
            logging.info(line)

    # Log stderr as error if there's any error output
    if result.stderr.strip():
        for line in result.stderr.strip().split('\n'):
            logging.error(line)

    # Check return code
    if result.returncode != 0:
        logging.error("Failed to make AmpliconClassifier input file!")
        sys.exit(1)

    with open(input_file) as ifile:
        sample_info_dict["number_of_AA_amplicons"] = len(ifile.readlines())

    cmd = "{} {}/amplicon_classifier.py -i {} --ref {} -o {}".format(PY3_PATH, AC_src, input_file, ref, class_output)
    logging.info(cmd + "\n")
    call(cmd, shell=True)
    metadata_dict["AC_cmd"] = cmd

    # Get AC version
    AC_version = \
    Popen([PY3_PATH, AC_src + "/amplicon_classifier.py", "--version"], stdout=PIPE, stderr=PIPE, universal_newlines=True).communicate()[
        0].rstrip()

    metadata_dict["AC_version"] = AC_version

    # iterate over the bed files and count anything that isn't "unknown" as a feature
    feat_count = 0
    if os.path.exists(bed_dir):
        for bf in os.listdir(bed_dir):
            if not "unknown" in bf and bf.endswith(".bed"):
                feat_count += 1

    sample_info_dict["number_of_AA_features"] = feat_count


def make_AC_table(sname, AC_outdir, AC_src, run_metadata_file, sample_metadata_file, ref, cnv_bed=None):
    # make the AC output table
    logging.info("Running AC make_results_table.py")
    class_output = AC_outdir + sname
    input_file = class_output + ".input"
    summary_map_file = class_output + "_summary_map.txt"
    classification_file = class_output + "_amplicon_classification_profiles.tsv"
    cmd = "{} {}/make_results_table.py -i {} --classification_file {} --summary_map {} --ref {}".format(
        PY3_PATH, AC_src, input_file, classification_file, summary_map_file, ref)

    if cnv_bed:
        cmd += " --cnv_bed " + cnv_bed

    if run_metadata_file:
        cmd += " --run_metadata_file " + run_metadata_file

    if sample_metadata_file:
        cmd += " --sample_metadata_file " + sample_metadata_file

    logging.debug("Command: " + cmd)

    try:
        result = run(cmd, shell=True, stdout=PIPE, stderr=PIPE, universal_newlines=True)

        # Log stdout if there's any output
        if result.stdout.strip():
            for line in result.stdout.strip().split('\n'):
                logging.info(line)

        # Log stderr if there's any output
        if result.stderr.strip():
            for line in result.stderr.strip().split('\n'):
                logging.warning(line)

        # Check return code
        if result.returncode != 0:
            logging.error("make_results_table.py failed with return code {}".format(result.returncode))

    except Exception as e:
        logging.error("Failed to run make_results_table.py: {}".format(str(e)))


def get_ref_sizes(ref_genome_size_file):
    chr_sizes = {}
    with open(ref_genome_size_file) as infile:
        for line in infile:
            fields = line.rstrip().rsplit()
            if fields:
                chr_sizes[fields[0]] = str(int(fields[1]) - 1)

    return chr_sizes


def get_ref_centromeres(ref_name):
    centromere_dict = {}
    fnameD = {"GRCh38": "GRCh38_centromere.bed", "GRCh37": "human_g1k_v37_centromere.bed",
              "hg19": "hg19_centromere.bed",
              "mm10": "mm10_centromere.bed", "GRCm38": "GRCm38_centromere.bed", "GRCh38_viral": "GRCh38_centromere.bed"}
    with open(AA_REPO + ref_name + "/" + fnameD[ref_name]) as infile:
        for line in infile:
            if not "centromere" in line and not "acen" in line:
                continue
            fields = line.rstrip().rsplit("\t")
            if fields[0] not in centromere_dict:
                centromere_dict[fields[0]] = (fields[1], fields[2])

            else:
                pmin = min(int(centromere_dict[fields[0]][0]), int(fields[1]))
                pmax = max(int(centromere_dict[fields[0]][1]), int(fields[2]))
                # pad with 20kb to avoid freebayes issues in calling near centromeres
                centromere_dict[fields[0]] = (str(pmin - 20000), str(pmax + 20000))

    return centromere_dict


def save_run_metadata(outdir, sname, args, launchtime, commandstring):
    # make a dictionary that stores
    # datetime
    # hostname
    # ref
    # PAA command
    # AA python interpreter version
    # bwa cmd
    # CN cmd
    # AA cmd
    # PAA version
    # CNVKit version
    # AA version
    # AC version
    metadata_dict["launch_datetime"] = launchtime
    metadata_dict["hostname"] = socket.gethostname()
    metadata_dict["ref_genome"] = args.ref
    aapint = args.aa_python_interpreter
    aa_python_v = Popen([aapint, "--version"], stdout=PIPE, stderr=PIPE, universal_newlines=True).communicate()[0].rstrip()
    samtools_version = get_samtools_version(args.samtools_path)
    # try:
    #     aa_python_v = aa_python_v.decode('utf-8')
    # except UnicodeError:
    #     pass

    metadata_dict["AA_python_version"] = aa_python_v
    metadata_dict["AmpliconSuite-pipeline_command"] = commandstring
    metadata_dict["AmpliconSuite-pipeline_version"] = __ampliconsuitepipeline_version__
    metadata_dict["Samtools version"] = "{}.{}".format(samtools_version[0], samtools_version[1])

    for x in ["bwa_cmd", "cnvkit_cmd", "amplified_intervals_cmd", "AA_cmd", "AC_cmd", "cnvkit_version", "AA_version",
              "AC_version"]:
        if x not in metadata_dict:
            metadata_dict[x] = "NA"

    # save the json dict
    run_metadata_filename = outdir + sname + "_run_metadata.json"
    with open(run_metadata_filename, 'w') as fp:
        json.dump(metadata_dict, fp, indent=2)

    # sample_info_dict["run_metadata_file"] = run_metadata_filename
    return run_metadata_filename


def detect_run_failure(align_stderr_file, AA_outdir, sname, AC_outdir, completed_AA_runs):
    if align_stderr_file:
        cmd = 'grep -i error ' + align_stderr_file
        try:
            aln_errs = check_output(cmd, shell=True).decode("utf-8")

        except CalledProcessError:
            aln_errs = ""

        if aln_errs:
            logging.error("Detected error during bwa mem alignment stage\n")
            return True

    if AA_outdir and not completed_AA_runs:
        sumfile = AA_outdir + sname + "_summary.txt"
        if os.path.isfile(sumfile):
            namps = -1
            with open(sumfile) as infile:
                for line in infile:
                    if line.startswith("#Amplicons = "):
                        namps = int(line.rstrip().rsplit(" = ")[-1])
                        break

            if namps < 0:
                logging.error("Detected truncated or missing AA outputs")
                return True

            for x in range(1, namps + 1):
                try:
                    fsize = os.stat(AA_outdir + sname + "_amplicon" + str(x) + "_cycles.txt").st_size

                except OSError:
                    fsize = 0

                if fsize == 0:
                    logging.error("Detected truncated or missing AA outputs")
                    return True

        else:
            logging.error("Detected error during AA stage")
            return True

    if AC_outdir:
        try:
            fsize1 = os.stat(AC_outdir + sname + "_amplicon_classification_profiles.tsv").st_size
            fsize2 = os.stat(AC_outdir + sname + "_result_table.tsv").st_size

        except OSError:
            fsize1 = 0
            fsize2 = 0

        if fsize1 == 0 or fsize2 == 0:
            logging.error("Detected error during AC stage\n")
            return True

    return False


def contains_spaces(file_path):
    return any(char == ' ' for char in file_path)


def handle_completed_aa_runs(completed_aa_path, output_dir, sample_name):
    """
    Handle completed AA runs - extract if archive, or use directory directly

    Returns:
        str: Path to AA output directory
    """
    if os.path.isdir(completed_aa_path):
        # Already a directory, use as-is
        return completed_aa_path

    elif completed_aa_path.endswith('.tar.gz') or completed_aa_path.endswith('.tgz'):
        # Extract tar.gz
        extract_dir = os.path.join(output_dir, sample_name + "_extracted_AA")
        os.makedirs(extract_dir, exist_ok=True)

        with tarfile.open(completed_aa_path, 'r:gz') as tar:
            tar.extractall(path=extract_dir)

        logging.info("Extracted tar.gz to: {}".format(extract_dir))
        return extract_dir

    elif completed_aa_path.endswith('.zip'):
        # Extract zip
        extract_dir = os.path.join(output_dir, sample_name + "_extracted_AA")
        os.makedirs(extract_dir, exist_ok=True)

        with zipfile.ZipFile(completed_aa_path, 'r') as zip_ref:
            zip_ref.extractall(extract_dir)

        logging.info("Extracted zip to: {}".format(extract_dir))
        return extract_dir

    else:
        logging.error("Unsupported file type for completed AA runs: {}".format(completed_aa_path))
        sys.exit(1)


def remove_existing_classifications(aa_outdir):
    """
    Remove existing *_classification/ directories that contain amplicon classification profiles
    Searches recursively and respects permission errors
    """
    if not os.path.exists(aa_outdir):
        logging.warning("AA output directory does not exist: {}".format(aa_outdir))
        return

    classification_dirs_found = []

    try:
        # Walk recursively through all subdirectories
        for root, dirs, files in os.walk(aa_outdir):
            for dirname in dirs:
                if dirname.endswith('_classification'):
                    dir_path = os.path.join(root, dirname)

                    # Check if directory contains amplicon classification profiles
                    try:
                        has_profile = False
                        for file in os.listdir(dir_path):
                            if file.endswith('_amplicon_classification_profiles.tsv'):
                                has_profile = True
                                break  # Break inner file loop only

                        if has_profile:
                            classification_dirs_found.append(dir_path)

                    except (OSError, PermissionError) as e:
                        logging.warning("Cannot access directory {}: {}".format(dir_path, str(e)))
                        continue

        # Remove found directories
        for dir_path in classification_dirs_found:
            try:
                logging.info("Removing existing classification directory: {}".format(dir_path))
                shutil.rmtree(dir_path)
            except PermissionError:
                logging.error("Permission denied: Cannot remove {}. May be owned by different user.".format(dir_path))
            except OSError as e:
                logging.error("Failed to remove {}: {}".format(dir_path, str(e)))

    except (OSError, PermissionError) as e:
        logging.error("Cannot access AA output directory {}: {}".format(aa_outdir, str(e)))
        return

    if classification_dirs_found:
        logging.info("Processed {} classification directories".format(len(classification_dirs_found)))


# MAIN #
def main():
    """Main entry point for AmpliconSuite-pipeline"""
    global args, AA_REPO, AA_SRC, AC_SRC, metadata_dict, ref_genome_size_file, PY3_PATH

    # Parse arguments using the new modular approach
    parser = setup_argument_parser()
    args = parser.parse_args()

    # Handle special cases that exit early
    if args.download_repo:
        handle_repo_download(args, AA_REPO)
        return

    # Start timing
    ta = time.time()
    ti = ta
    launchtime = str(datetime.now())

    # Validate arguments
    validate_arguments(args, parser)

    # Initialize logging and directories
    paa_logfile, timing_logfile, commandstring, finish_flag_filename = initialize_logging_and_directories(args, launchtime)

    logging.info("Running initial checks and configurations...")

    # Setup environment and validate using the new modules
    AA_REPO = setup_environment_and_paths(args)

    # Setup tool paths
    setup_tool_paths(args)

    # Set PY3_PATH if provided
    if args.python3_path:
        PY3_PATH = args.python3_path

    # Validate AA environment
    AA_SRC, AC_SRC = validate_aa_environment(args)

    # Initialize the metadata_dict here so functions can access it
    metadata_dict = {}

    # Create coverage stats file if needed
    create_coverage_stats_file(AA_REPO)

    # Validate samtools
    samtools_version = get_samtools_version(args.samtools_path)

    # Setup reference genome info early so functions can access it
    refFnames = {x: None for x in ["hg19", "GRCh37", "GRCh38", "GRCh38_viral", "mm10"]}
    for rname in refFnames.keys():
        if os.path.exists(AA_REPO + "/" + rname):
            refFnames[rname] = check_reference.get_ref_fname(AA_REPO, rname)

    try:
        # Now run your existing pipeline logic
        run_pipeline_logic(paa_logfile, timing_logfile, ta, ti, launchtime, commandstring, samtools_version, refFnames, finish_flag_filename)

    finally:
        if timing_logfile:
            timing_logfile.close()


def run_pipeline_logic(paa_logfile, timing_logfile, ta, ti, launchtime, commandstring, samtools_version, refFnames, finish_flag_filename):
    """Your existing pipeline logic from the main function"""
    global ref_genome_size_file

    sname = args.sample_name
    outdir = args.output_directory
    sample_metadata_filename = args.output_directory + sname + "_sample_metadata.json"

    # Load sample metadata
    with open(args.sample_metadata) as input_json:
        sample_info_dict = json.load(input_json)

    sample_info_dict["reference_genome"] = args.ref
    sample_info_dict["sample_name"] = sname

    # Handle BAM file reference checking
    faidict = {}
    if args.bam:
        if args.ref and refFnames[args.ref] is None:
            em = "Data repo files for ref " + args.ref + " not found. Please download using the '--download_repo " + args.ref + "' option\n"
            logging.error(em)
            sys.exit(1)
        else:
            for k, v in refFnames.items():
                if v:
                    faidict[k] = AA_REPO + k + "/" + v + ".fai"

        determined_ref = check_reference.check_ref(args.bam, faidict, args.samtools_path)
        if not determined_ref and not args.ref:
            logging.error("Could not determine ref build. Please make sure AA data repo is populated.")
            sys.exit(1)
        elif not args.ref:
            args.ref = determined_ref
        elif args.ref and not determined_ref:
            logging.error("ERROR! The BAM file could not be matched to an available reference genome in AA_DATA_REPO! Is it aligned to an AA-supported reference genome?")
        elif args.ref != determined_ref:
            logging.warning("WARNING: User specified --ref {}, however the bam header matches instead to build {}. Continuing with {}...".format(args.ref, determined_ref, determined_ref))
            args.ref = determined_ref

    # Check data repo freshness
    try:
        with open(AA_REPO + args.ref + "/last_updated.txt", 'r') as file:
            datestring = file.read()
            logging.info(args.ref + " data repo constructed on " + datestring)
    except FileNotFoundError:
        logging.warning("Data repo appears to be out of date. Please update your data repo!\n")

    # Setup file paths - NOW we can properly set ref_genome_size_file
    gdir = AA_REPO + args.ref + "/"
    ref_fasta = gdir + refFnames[args.ref]
    ref_genome_size_file = gdir + args.ref + "_noAlt.fa.fai"

    # Log initialization timing
    tb = time.time()
    timing_logfile.write("Initialization:\t" + "{:.2f}".format(tb - ta) + "\n")
    ta = tb
    logging.info("Running AmpliconSuite-pipeline on sample: " + sname)

    # Begin pipeline execution - copy your existing main function code from here
    aln_stage_stderr = None

    # Alignment stage
    if args.fastqs:
        if args.fastqs[0] == args.fastqs[1]:
            logging.error(str(args.fastqs))
            logging.error("You must provide two different fastq files for paired-end reads!\n")
            sys.exit(1)
        elif contains_spaces(args.fastqs[0]) or contains_spaces(args.fastqs[1]):
            logging.error("FASTQ filepaths cannot contain spaces!")
            sys.exit(1)
        elif not os.path.exists(args.fastqs[0]) or not os.path.exists(args.fastqs[1]):
            logging.error("One or both FASTQ files do not exist!")
            sys.exit(1)

        fastqs = " ".join(args.fastqs)
        logging.info("Will perform alignment on " + fastqs)
        args.bam, aln_stage_stderr = run_bwa(ref_fasta, fastqs, outdir, sname, args.nthreads, args.samtools_path,
                                             samtools_version)

    AA_outdir = None
    AC_outdir = None
    if not args.completed_AA_runs:
        # BAM indexing
        bamBaiNoExt = args.bam[:-3] + "bai"
        cramCraiNoExt = args.bam[:-4] + "crai"
        baiExists = os.path.isfile(args.bam + ".bai") or os.path.isfile(bamBaiNoExt)
        craiExists = os.path.isfile(args.bam + ".crai") or os.path.isfile(cramCraiNoExt)
        if not baiExists and not craiExists:
            logging.info(args.bam + " index not found, calling samtools index")
            call([args.samtools_path, "index", args.bam])
            logging.info("Finished indexing")

        bambase = os.path.splitext(os.path.basename(args.bam))[0]
        prop_paired_proportion = None
        if not args.no_QC:
            prop_paired_proportion = check_reference.check_properly_paired(args.bam, args.samtools_path)

        tb = time.time()
        timing_logfile.write("Alignment, indexing and QC:\t" + "{:.2f}".format(tb - ta) + "\n")

        if args.align_only:
            logging.info("Completed\n")
            tf = time.time()
            timing_logfile.write("Total_elapsed_walltime\t" + "{:.2f}".format(tf - ti) + "\n")
            return

        ta = tb
        centromere_dict = get_ref_centromeres(args.ref)
        chr_sizes = get_ref_sizes(ref_genome_size_file)

        # CNV calling stage
        cnvkit_output_directory = None
        runCNV = None
        if args.cnv_bed and args.cnv_bed.endswith(".cns"):
            args.cnv_bed = convert_cnvkit_cns_to_bed(outdir, bambase, cnsfile=args.cnv_bed, nofilter=True)
            runCNV = "CNSfile"

        elif not args.cnv_bed:
            runCNV = "CNVkit"
            cnvkit_output_directory = args.output_directory + sname + "_cnvkit_output/"
            if not os.path.exists(cnvkit_output_directory):
                os.mkdir(cnvkit_output_directory)

            # Your run_cnvkit function can now access all globals it needs
            run_cnvkit(args.cnvkit_dir, args.nthreads, cnvkit_output_directory, args.bam,
                       seg_meth=args.cnvkit_segmentation, normal=args.normal_bam, ref_fasta=ref_fasta)

            if args.ploidy or args.purity:
                rescale_cnvkit_calls(args.cnvkit_dir, cnvkit_output_directory, bambase,
                                     ploidy=args.ploidy, purity=args.purity)
                rescaling = True
            else:
                rescaling = False

            args.cnv_bed = convert_cnvkit_cns_to_bed(cnvkit_output_directory, bambase, rescaled=rescaling)
            sample_info_dict["sample_cnv_bed"] = args.cnv_bed

        if runCNV:
            # Custom CNV plotting
            centromeres = cnv_plots.load_centromere_file(AA_REPO, args.ref)
            if centromeres is not None:
                logging.info(f"Loaded {len(centromeres)} centromere regions for highlighting")

            if sample_info_dict["sample_cnv_bed"].endswith(".bed"):
                logging.info("Plotting CNV distribution across chromosomes: {}{}_cnv_distribution.png/.pdf".format(
                    cnvkit_output_directory,sname))
                cnv_data = cnv_plots.load_cnv_bed_file(sample_info_dict["sample_cnv_bed"])
                cnv_plots.plot_cnv_distribution_chromosomes(cnv_data, bambase, "{}{}_cnv_distribution".format(
                    cnvkit_output_directory,sname), centromeres=centromeres)
            else:
                logging.warning(
                    "Skipping plotting CNV distribution across chromosomes, as the provided CNV bed file is not in the expected format.")

        tb = time.time()
        timing_logfile.write("CNV calling:\t" + "{:.2f}".format(tb - ta) + "\n")

        ta = tb
        # Seed filtering stage
        if not args.no_filter and not args.cnv_bed.endswith("_AA_CNV_SEEDS.bed"):
            if not args.cnv_bed.endswith("_CNV_CALLS_pre_filtered.bed") and not args.cnv_bed.endswith(
                    "_CNV_CALLS_unfiltered_gains.bed"):
                from paalib import cnv_prefilter
                pfilt_odir = cnvkit_output_directory if cnvkit_output_directory else args.output_directory
                args.cnv_bed = cnv_prefilter.prefilter_bed(args.cnv_bed, args.ref, centromere_dict, chr_sizes,
                                                           args.cngain, pfilt_odir)

            amplified_interval_bed = run_amplified_intervals(args.aa_python_interpreter, args.cnv_bed, args.bam,
                                                             outdir, sname, args.cngain, args.cnsize_min, paa_logfile)

        elif args.no_filter and runCNV:
            if not args.cnv_bed.endswith("_CNV_CALLS_pre_filtered.bed") and not args.cnv_bed.endswith(
                    "_CNV_CALLS_unfiltered_gains.bed"):
                from paalib import cnv_prefilter
                pfilt_odir = cnvkit_output_directory if cnvkit_output_directory else args.output_directory
                args.cnv_bed = cnv_prefilter.prefilter_bed(args.cnv_bed, args.ref, centromere_dict, chr_sizes,
                                                           args.cngain, pfilt_odir)
                logging.info("Skipping amplified_intervals.py step due to --no_filter")

            amplified_interval_bed = args.cnv_bed
        else:
            logging.info("Skipping filtering of bed file.")
            amplified_interval_bed = args.cnv_bed

        tb = time.time()
        timing_logfile.write("Seed filtering (amplified_intervals.py):\t" + "{:.2f}".format(tb - ta) + "\n")
        ta = tb

        # AmpliconArchitect stage
        if args.run_AA:
            AA_outdir = outdir + sname + "_AA_results/"
            if not os.path.exists(AA_outdir):
                os.mkdir(AA_outdir)

            # Set insert sdevs if not given by user
            if not args.no_QC and not args.AA_insert_sdevs and prop_paired_proportion is not None and prop_paired_proportion < 90:
                logging.info("Properly paired rate less than 90%, setting --insert_sdevs 9.0 for AA")
                args.AA_insert_sdevs = 9.0

            run_AA(amplified_interval_bed, AA_outdir, sname, args)
            tb = time.time()
            timing_logfile.write("AmpliconArchitect:\t" + "{:.2f}".format(tb - ta) + "\n")
            ta = tb

            # AmpliconClassifier stage
            if args.run_AC:
                AC_outdir = outdir + sname + "_classification/"
                if not os.path.exists(AC_outdir):
                    os.mkdir(AC_outdir)

                run_AC(AA_outdir, sname, args.ref, AC_outdir, AC_SRC)
                tb = time.time()
                timing_logfile.write("AmpliconClassifier:\t" + "{:.2f}".format(tb - ta) + "\n")

        # Save metadata
        run_metadata_filename = save_run_metadata(outdir, sname, args, launchtime, commandstring)

        with open(sample_metadata_filename, 'w') as fp:
            json.dump(sample_info_dict, fp, indent=2)

        if args.run_AA and args.run_AC:
            make_AC_table(sname, AC_outdir, AC_SRC, run_metadata_filename, sample_metadata_filename,
                          args.ref, cnv_bed=sample_info_dict.get("sample_cnv_bed", None))

    else:  # Handle completed AA runs
        cnvkit_output_directory = None
        run_metadata_filename = args.completed_run_metadata
        amplified_interval_bed = None

        if not args.ref:
            logging.error("--ref is a required argument if --completed_AA_runs is provided!")
            sys.exit(1)

        if not args.upload and not args.run_AC:
            logging.error(
                "If --completed_AA_runs provided, either --run_AC and/or --upload (and associated upload args) must also be set!")
            sys.exit(1)

        # Handle archive extraction or use directory directly
        AA_outdir = handle_completed_aa_runs(args.completed_AA_runs, outdir, sname)

        if args.run_AC:
            # Remove existing classification directories to avoid conflicts
            remove_existing_classifications(AA_outdir)

            AC_outdir = outdir + sname + "_classification/"
            if not os.path.exists(AC_outdir):
                os.mkdir(AC_outdir)

            ta = time.time()
            run_AC(AA_outdir, sname, args.ref, AC_outdir, AC_SRC)

            with open(sample_metadata_filename, 'w') as fp:
                json.dump(sample_info_dict, fp, indent=2)

            make_AC_table(sname, AC_outdir, AC_SRC, run_metadata_filename, sample_metadata_filename, args.ref)
            tb = time.time()
            timing_logfile.write("AmpliconClassifier:\t" + "{:.2f}".format(tb - ta) + "\n")

    # Final checks and cleanup
    run_failure = detect_run_failure(aln_stage_stderr, AA_outdir, sname, AC_outdir, args.completed_AA_runs)
    if not run_failure:
        logging.info("All stages appear to have completed successfully.")
        with open(finish_flag_filename, 'w') as ffof:
            ffof.write("All stages completed\n")

    else:
        logging.error("One or more stages did not complete successfully.")

    if args.upload:
        if not run_failure:
            ta = time.time()
            final_cnv_bed = sample_info_dict.get("sample_cnv_bed", None)
            archive_and_upload_sample(AA_outdir, AC_outdir, final_cnv_bed, cnvkit_output_directory,
                                      run_metadata_filename, sample_metadata_filename, amplified_interval_bed,
                                      finish_flag_filename, args.project_uuid, args.project_key, args.username, outdir + sname,
                                      server=args.upload_server)
            tb = time.time()
            timing_logfile.write("run_upload:\t" + "{:.2f}".format(tb - ta) + "\n")

        else:
            logging.error("Skipped run upload due to detected run failure!")

    # Final timing
    tf = time.time()
    timing_logfile.write("Total_elapsed_walltime\t" + "{:.2f}".format(tf - ti) + "\n")


if __name__ == '__main__':
    main()
