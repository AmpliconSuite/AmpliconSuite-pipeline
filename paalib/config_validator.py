"""
Configuration validation and environment setup for AmpliconSuite-pipeline
"""

from datetime import datetime
import logging
import os
import re
from subprocess import *
import sys

from paalib import check_reference


def validate_arguments(args, parser):
    """Validate argument combinations and required fields"""
    # Required arguments
    if not args.sample_name:
        parser.error("--sample_name (-s) is a required argument.")
    if not args.nthreads:
        parser.error("--nthreads (-t) is a required argument.")
    if not any([args.bam, args.fastqs, args.completed_AA_runs]):
        parser.error("One of --bam | --fastqs | --completed_AA_runs is required!")

    # Sample name validation
    if "/" in args.sample_name:
        logging.error("Sample name -s cannot be a path. Specify output directory with -o.\n")
        sys.exit(1)

    # Fastq/completed runs validation
    if (args.fastqs or args.completed_AA_runs) and not args.ref:
        logging.error("Must specify --ref when providing unaligned fastq files or completed AA runs.\n")
        sys.exit(1)

    # Python interpreter validation
    if args.aa_python_interpreter and not any(
            args.aa_python_interpreter.endswith(x) for x in ['python', 'python2', 'python3']):
        logging.error("--aa_python_interpreter must be a path of a valid python interpreter")
        sys.exit(1)

    # Upload arguments validation - must all be provided or none
    upload_args = [args.project_uuid, args.project_key, args.username]
    upload_args_provided = [arg for arg in upload_args if arg is not None]

    if args.upload:
        # If --upload flag is set, all upload arguments must be provided
        if len(upload_args_provided) != len(upload_args):
            missing_args = []
            if not args.project_uuid:
                missing_args.append("--project_uuid")
            if not args.project_key:
                missing_args.append("--project_key")
            if not args.username:
                missing_args.append("--username")

            parser.error("--upload flag is set but missing required arguments: {}".format(", ".join(missing_args)))

    else:
        # If --upload flag is not set, warn if any upload arguments are provided
        if upload_args_provided:
            provided_args = []
            if args.project_uuid:
                provided_args.append("--project_uuid")
            if args.project_key:
                provided_args.append("--project_key")
            if args.username:
                provided_args.append("--username")

            logging.warning(
                "Upload arguments provided but --upload flag not set. Ignoring: {}".format(", ".join(provided_args)))

    # File existence checks
    if args.cnv_bed and not os.path.isfile(args.cnv_bed):
        logging.error("Specified CNV bed file does not exist: " + args.cnv_bed + "\n")
        sys.exit(1)

    # Check for spaces in file paths
    if args.bam and _contains_spaces(args.bam):
        logging.error("BAM filepath cannot contain spaces!")
        sys.exit(1)

    if args.fastqs:
        if args.fastqs[0] == args.fastqs[1]:
            logging.error(str(args.fastqs))
            logging.error("You must provide two different fastq files for paired-end reads!\n")
            sys.exit(1)
        elif _contains_spaces(args.fastqs[0]) or _contains_spaces(args.fastqs[1]):
            logging.error("FASTQ filepaths cannot contain spaces!")
            sys.exit(1)
        elif not os.path.exists(args.fastqs[0]) or not os.path.exists(args.fastqs[1]):
            logging.error("One or both FASTQ files do not exist!")
            sys.exit(1)


def setup_environment_and_paths(args):
    """Setup environment variables, paths, and basic configurations"""
    # Check AA_REPO
    try:
        AA_REPO = os.environ['AA_DATA_REPO'] + "/"
    except KeyError:
        sys.stderr.write("AA_DATA_REPO bash variable not found. Please see installation instructions and run ./install.sh before using.\n")
        sys.exit(1)
    
    # Set samtools path
    if not args.samtools_path.endswith("/samtools"):
        if args.samtools_path and not args.samtools_path.endswith("/"):
            args.samtools_path += "/"
        args.samtools_path += "samtools"
    
    # Handle reference genome aliases
    if args.ref == "hg38":
        args.ref = "GRCh38"
    if args.ref == "GRCm38":
        args.ref = "mm10"
    
    # Handle completed run metadata
    if hasattr(args, 'completed_run_metadata') and args.completed_run_metadata.lower() == "none":
        args.completed_run_metadata = None
    
    # Handle sample metadata
    if not args.sample_metadata:
        args.sample_metadata = os.path.realpath(os.path.dirname(check_reference.__file__)) + "/sample_metadata_skeleton.json"
    
    return AA_REPO


def setup_tool_paths(args):
    """Setup and validate tool paths"""
    # CNVKit path setup
    if not (args.cnv_bed or args.cnvkit_dir or getattr(args, 'completed_run_metadata', None) or args.align_only) and (args.fastqs or args.bam):
        try:
            args.cnvkit_dir = str(check_output(["which cnvkit.py"], shell=True).decode("utf-8").rstrip())
        except CalledProcessError:
            logging.error("cnvkit.py not found on system path. Must specify --cnvkit_dir")
            sys.exit(1)
    elif args.cnvkit_dir and not args.cnvkit_dir.endswith("/") and not args.cnvkit_dir.endswith("cnvkit.py"):
        args.cnvkit_dir += "/"
    
    if hasattr(args, 'cnvkit_dir') and args.cnvkit_dir and not args.cnvkit_dir.endswith("cnvkit.py"):
        args.cnvkit_dir += "cnvkit.py"
    
    # Python3 path setup
    if args.python3_path:
        if not args.python3_path.endswith("/python") and not args.python3_path.endswith("/python3"):
            args.python3_path += "/python3"
    
    # Rscript validation for CNVKit
    if hasattr(args, 'cnvkit_dir') and args.cnvkit_dir and not getattr(args, 'cnv_bed', None):
        test_rscript = "Rscript"
        if args.rscript_path:
            if not args.rscript_path.endswith("/Rscript"):
                args.rscript_path += "/Rscript"
            test_rscript = args.rscript_path
        
        try:
            rscript_version_out = str(check_output([test_rscript, "--version"], stderr=STDOUT).decode("utf-8").rstrip())
        except CalledProcessError:
            logging.error(test_rscript + " not found. Must specify --rscript_path")
            sys.exit(1)

        # Check if DNAcopy package is installed
        try:
            dnacopy_check = str(check_output([test_rscript, "-e", 'cat(require("DNAcopy", quietly=TRUE))'],
                                         stderr=STDOUT).decode("utf-8").strip())
            if not dnacopy_check.startswith("TRUE"):
                logging.error(
                    "DNAcopy R package is required for CNVKit. Please install it with: Rscript -e 'BiocManager::install(\"DNAcopy\")'")
                sys.exit(1)
            else:
                logging.info("Found DNAcopy R package installed")
        except CalledProcessError as e:
            logging.error("Failed to check for DNAcopy R package. Error: " + str(e))
            sys.exit(1)


def validate_aa_environment(args):
    """Validate AmpliconArchitect environment and licensing"""
    # Set AA_SRC if provided
    if args.AA_src:
        os.environ['AA_SRC'] = args.AA_src
    
    # Find AA_SRC
    try:
        AA_SRC = os.environ['AA_SRC']
    except KeyError:
        try:
            import ampliconarchitectlib
            AA_SRC = os.path.realpath(os.path.dirname(ampliconarchitectlib.__file__))
        except ModuleNotFoundError:
            logging.error("AA_SRC bash variable or library files not found. AmpliconArchitect may not be properly installed.\n")
            sys.exit(1)
    
    # Find AC_SRC
    try:
        AC_SRC = os.environ['AC_SRC']
    except KeyError:
        try:
            import ampclasslib
            ac_path = check_output("which amplicon_classifier.py", shell=True).decode("utf-8")
            AC_SRC = ac_path.rsplit("/amplicon_classifier.py")[0]
        except Exception as e:
            logging.error(e)
            logging.error("\nAC_SRC bash variable or library files not found. AmpliconClassifier may not be properly installed.\n")
            sys.exit(1)
    
    # Check MOSEK license if running AA
    if args.run_AA:
        _validate_mosek_license()
    
    return AA_SRC, AC_SRC


def initialize_logging_and_directories(args, launchtime):
    """Setup logging and create necessary directories"""
    sname = args.sample_name

    # Handle output directory
    if not args.output_directory:
        args.output_directory = os.getcwd()
    if not args.output_directory.endswith("/"):
        args.output_directory += "/"

    # Make output directory
    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)

    # Setup logging with different levels for file vs console
    paa_logfile = args.output_directory + sname + '.log'

    # Clear any existing handlers to avoid conflicts
    logging.getLogger().handlers.clear()

    # Set root logger to DEBUG to capture everything
    logging.getLogger().setLevel(logging.DEBUG)

    # File handler - DEBUG and higher
    file_handler = logging.FileHandler(paa_logfile, mode='w')
    file_handler.setLevel(logging.DEBUG)
    file_formatter = logging.Formatter('[%(name)s:%(levelname)s]\t%(message)s')
    file_handler.setFormatter(file_formatter)
    logging.getLogger().addHandler(file_handler)

    # Console handler - INFO and higher
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_formatter = logging.Formatter('[%(name)s:%(levelname)s]\t%(message)s')
    console_handler.setFormatter(console_formatter)
    logging.getLogger().addHandler(console_handler)

    # Log startup info
    from paalib._version import __ampliconsuitepipeline_version__
    logging.info("Launched on " + launchtime)
    logging.info("AmpliconSuite-pipeline version " + __ampliconsuitepipeline_version__ + "\n")

    # Log command
    commandstring = _build_command_string()
    logging.info("AmpliconSuite-pipeline command:")
    logging.info(commandstring + "\n")

    # Create timing log
    timing_logfile = open(args.output_directory + sname + '_timing_log.txt', 'w')
    timing_logfile.write("#stage:\twalltime(seconds)\n")

    # Setup finish flag
    finish_flag_filename = args.output_directory + args.sample_name + "_finish_flag.txt"
    if os.path.exists(finish_flag_filename):
        logging.warning("WARNING: Running AmpliconSuite-pipeline.py with outputs directed into the same output location"
                        " as a previous run may cause crashes or other unexpected behavior. To avoid errors, clear "
                        "previous files before re-running.\n")

    with open(finish_flag_filename, 'w') as ffof:
        ffof.write("UNSUCCESSFUL\n")

    logging.getLogger('fontTools.subset').level = logging.WARN
    logging.getLogger('fontTools.ttLib').level = logging.WARN
    logging.getLogger('matplotlib.backends').level = logging.WARN
    logging.getLogger('matplotlib.font_manager').level = logging.WARN

    return paa_logfile, timing_logfile, commandstring, finish_flag_filename


def create_coverage_stats_file(AA_REPO):
    """Create coverage.stats file if it doesn't exist"""
    if not os.path.exists(os.path.join(AA_REPO, "coverage.stats")):
        logging.info("coverage.stats file not found in " + AA_REPO + "\nCreating a new coverage.stats file.")
        from subprocess import call
        cmd = "touch {}coverage.stats && chmod a+rw {}coverage.stats".format(AA_REPO, AA_REPO)
        logging.info(cmd)
        call(cmd, shell=True)


def get_samtools_version(samtools):
    try:
        # Run the command to get the version information
        result = Popen([samtools], stderr=PIPE, stdout=PIPE, universal_newlines=True)
        _, output = result.communicate()

        # Decode the output if it's in bytes (Python 3)
        # if isinstance(output, bytes):
        #     output = output.decode('utf-8')

        # Parse the version information to extract major and minor versions
        version_pattern = r'Version: (\d+)\.(\d+)'
        match = re.search(version_pattern, output)
        if match:
            major_version = int(match.group(1))
            minor_version = int(match.group(2))
            return major_version, minor_version
        else:
            # Return None if version information couldn't be parsed
            return None, None
    except OSError as e:
        # Handle the case when Samtools is not found
        logging.error("Error: Samtools not found. Please make sure it is installed and in your PATH.")
        return None, None


def _contains_spaces(filepath):
    """Check if filepath contains spaces"""
    return ' ' in filepath


def _validate_mosek_license():
    """Validate MOSEK license for AmpliconArchitect"""
    mosek_license_path = None
    
    # Check if license exists in default location
    if os.path.exists(os.environ["HOME"] + "/mosek/mosek.lic"):
        mosek_license_path = os.environ["HOME"] + "/mosek/mosek.lic"
    
    # Check if license exists in location specified by MOSEKLM_LICENSE_FILE
    elif "MOSEKLM_LICENSE_FILE" in os.environ:
        if os.environ["MOSEKLM_LICENSE_FILE"].endswith("mosek.lic"):
            logging.error(
                "MOSEKLM_LICENSE_FILE should be the path of the directory of the license, not the full path. Please update your .bashrc, and run 'source ~/.bashrc'")
            sys.exit(1)
        elif os.path.exists(os.environ["MOSEKLM_LICENSE_FILE"] + "/mosek.lic"):
            mosek_license_path = os.environ["MOSEKLM_LICENSE_FILE"] + "/mosek.lic"
        else:
            logging.error("--run_AA set, but MOSEK license not found in " + os.environ["MOSEKLM_LICENSE_FILE"])
            sys.exit(1)
    else:
        logging.error("--run_AA set, but MOSEK license not found in $HOME/mosek/")
        sys.exit(1)
    
    # Check license age
    if mosek_license_path:
        file_time = os.path.getmtime(mosek_license_path)
        file_date = datetime.fromtimestamp(file_time)
        current_date = datetime.now()
        days_old = (current_date - file_date).days
        
        # Check if license has expired (older than 365 days)
        if days_old >= 365:
            logging.warning("*" * 80)
            logging.warning("WARNING: MOSEK LICENSE IS EXPIRED!")
            logging.warning(
                "The Mosek license file at " + mosek_license_path + " is " + str(days_old) + " days old.")
            logging.warning("AA will not run with an expired license.")
            logging.warning("Please obtain an updated Mosek license to continue using AA.")
            logging.warning("*" * 80)
        # Check if license is about to expire (within 7 days)
        elif days_old >= 358:
            days_until_expiry = 365 - days_old
            logging.warning("*" * 80)
            logging.warning("WARNING: MOSEK LICENSE WILL EXPIRE SOON!")
            logging.warning("The Mosek license file at " + mosek_license_path + " will expire in " + str(
                days_until_expiry) + " days.")
            logging.warning("Please obtain an updated Mosek license before expiration to continue using AA.")
            logging.warning("*" * 80)


def _build_command_string():
    """Build command string for logging"""
    commandstring = ""
    for arg in sys.argv:
        if ' ' in arg:
            commandstring += '"{}" '.format(arg)
        else:
            commandstring += "{} ".format(arg)
    return commandstring