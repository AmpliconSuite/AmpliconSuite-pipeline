"""
Argument parsing functionality for AmpliconSuite-pipeline
"""
import argparse
from paalib._version import __ampliconsuitepipeline_version__


def setup_argument_parser():
    """Create and configure the argument parser for AmpliconSuite-pipeline"""
    parser = argparse.ArgumentParser(
        description="A pipeline wrapper for AmpliconArchitect, invoking alignment CNV calling and CNV filtering prior. "
                    "Can launch AA, as well as downstream amplicon classification.")
    
    _add_basic_arguments(parser)
    _add_io_arguments(parser)
    _add_pipeline_arguments(parser)
    _add_reference_arguments(parser)
    _add_cnv_arguments(parser)
    _add_tool_path_arguments(parser)
    _add_aa_arguments(parser)
    _add_cnvkit_arguments(parser)
    _add_utility_arguments(parser)
    _add_upload_arguments(parser)
    _add_mutually_exclusive_groups(parser)
    
    return parser


def _add_basic_arguments(parser):
    """Add basic arguments like version, sample name, threads"""
    parser.add_argument("-v", "--version", action='version',
                        version='AmpliconSuite-pipeline version {version} \n'.format(version=__ampliconsuitepipeline_version__))
    parser.add_argument("-s", "--sample_name", metavar='STR', help="(Required) Sample name")
    parser.add_argument("-t", "--nthreads", metavar='INT', help="(Required) Number of threads to use in BWA and CNV calling")


def _add_io_arguments(parser):
    """Add input/output related arguments"""
    parser.add_argument("-o", "--output_directory", metavar='PATH', 
                        help="output directory names (will create if not already created)")
    parser.add_argument("--sample_metadata", metavar='FILE', help="JSON file of sample metadata to build on")


def _add_pipeline_arguments(parser):
    """Add pipeline control arguments"""
    parser.add_argument("--run_AA", help="Run AA after all files prepared. Default off.", action='store_true')
    parser.add_argument("--run_AC", help="Run AmpliconClassifier after all files prepared. Default off.",
                        action='store_true')
    parser.add_argument("--no_QC", help="Skip QC on the BAM file. Do not adjust AA insert_sdevs for "
                                        "poor-quality insert size distribution", action='store_true')
    parser.add_argument("--align_only", help="Only perform the alignment stage (do not run CNV calling and seeding",
                        action='store_true')


def _add_reference_arguments(parser):
    """Add reference genome related arguments"""
    parser.add_argument("--ref", metavar='STR', help="Reference genome version. Autodetected unless fastqs given as input.",
                        choices=["hg19", "GRCh37", "GRCh38", "hg38", "mm10", "GRCm38", "GRCh38_viral"], type=str)
    parser.add_argument("--download_repo", help="Download the selected data repo to the $AA_DATA_REPO "
                        "directory and exit. '_indexed' suffix indicates BWA index is included, which is useful if "
                        "performing alignment with AmpliconSuite-pipeline, but has a larger filesize.", 
                        choices=["hg19", "GRCh37", "GRCh38", "mm10", "GRCh38_viral", "hg19_indexed", 
                                "GRCh37_indexed", "GRCh38_indexed", "mm10_indexed", "GRCh38_viral_indexed"], 
                        nargs='+')


def _add_cnv_arguments(parser):
    """Add CNV calling related arguments"""
    parser.add_argument("--cngain", metavar='FLOAT', type=float, help="CN gain threshold to consider for AA seeding",
                        default=4.5)
    parser.add_argument("--cnsize_min", metavar='INT', type=int, help="CN interval size (in bp) to consider for AA seeding",
                        default=50000)
    parser.add_argument("--no_filter", help="Do not run amplified_intervals.py to remove low confidence candidate seed"
                                            " regions overlapping repetitive parts of the genome", action='store_true')


def _add_tool_path_arguments(parser):
    """Add tool path arguments"""
    parser.add_argument("--rscript_path", metavar='PATH', help="Specify custom path to Rscript for CNVKit, "
                        "which requires R version >=3.5")
    parser.add_argument("--python3_path", metavar='PATH', help="If needed, specify a custom path to python3.")
    parser.add_argument("--aa_python_interpreter",
                        help="By default AmpliconSuite-pipeline will use the system's default python path. If you would like to use "
                             "a different python version with AA, set this to either the path to the interpreter or "
                             "'python', 'python3', 'python2' (default 'python3')", metavar='PATH', type=str, default='python3')
    parser.add_argument("--samtools_path", help="Path to samtools binary (e.g., /path/to/my/samtools). If unset, will use samtools on system path.",
                        default='')
    parser.add_argument("--AA_src", metavar='PATH', help="Specify a custom $AA_SRC path. Overrides the bash variable")


def _add_aa_arguments(parser):
    """Add AmpliconArchitect specific arguments"""
    parser.add_argument("--downsample", metavar='FLOAT', type=float, help="AA downsample argument (see AA documentation)",
                        default=10)
    parser.add_argument("--sv_vcf",
                        help="Provide a VCF file of externally-called SVs to augment SVs identified by AA internally.",
                        metavar='FILE', action='store', type=str)
    parser.add_argument("--sv_vcf_no_filter", help="Use all external SV calls from the --sv_vcf arg, even "
                        "those without 'PASS' in the FILTER column.", action='store_true', default=False)
    parser.add_argument("--AA_runmode", metavar='STR', help="If --run_AA selected, set the --runmode argument to AA. Default mode is "
                        "'FULL'", choices=['FULL', 'BPGRAPH', 'CYCLES', 'SVVIEW'], default='FULL')
    parser.add_argument("--AA_extendmode", metavar='STR', help="If --run_AA selected, set the --extendmode argument to AA. Default "
                        "mode is 'EXPLORE'", choices=["EXPLORE", "CLUSTERED", "UNCLUSTERED", "VIRAL"], default='EXPLORE')
    parser.add_argument("--AA_insert_sdevs", help="Number of standard deviations around the insert size. May need to "
                        "increase for sequencing runs with high variance after insert size selection step. (default "
                        "3.0)", metavar="FLOAT", type=float, default=None)
    parser.add_argument('--pair_support_min', dest='pair_support_min', help="Number of read pairs for "
                        "minimum breakpoint support (default 2 but typically becomes higher due to coverage-scaled "
                        "cutoffs)", metavar='INT', action='store', type=int)
    parser.add_argument('--foldback_pair_support_min', help="Number of read pairs for minimum foldback SV support "
                        "(default 2 but typically becomes higher due to coverage-scaled cutoffs). Used value will be the maximum"
                        " of pair_support and this argument. Raising to 3 will help dramatically in heavily artifacted samples.",
                        metavar='INT', action='store', type=int)


def _add_cnvkit_arguments(parser):
    """Add CNVKit specific arguments"""
    parser.add_argument("--normal_bam", metavar='FILE', help="Path to matched normal bam for CNVKit (optional)")
    parser.add_argument("--ploidy", metavar='FLOAT', type=float, help="Ploidy estimate for CNVKit (optional). This is not used outside of CNVKit.",
                        default=None)
    parser.add_argument("--purity", metavar='FLOAT', type=float, help="Tumor purity estimate for CNVKit (optional). This is not used outside of CNVKit.",
                        default=None)
    parser.add_argument("--cnvkit_segmentation", metavar='STR', help="Segmentation method for CNVKit (if used), defaults to CNVKit "
                        "default segmentation method (cbs).", choices=['cbs', 'haar', 'hmm', 'hmm-tumor', 'hmm-germline', 'none'],
                        default='cbs')


def _add_utility_arguments(parser):
    """Add utility and miscellaneous arguments"""
    parser.add_argument("--completed_run_metadata", metavar='FILE',
                        help="Run metadata JSON to retroactively assign to collection of samples", default="")


def _add_upload_arguments(parser):
    """Add upload-related arguments"""
    parser.add_argument("--upload",
                        help="Upload sample results to AmpliconRepository (requires --run_AA and --run_AC)",
                        action='store_true')
    parser.add_argument("--project_uuid", metavar='STR',
                        help="Project UUID for upload (required if --upload is set)")
    parser.add_argument("--project_key", metavar='STR', help="Project key for upload (required if --upload is set)")
    parser.add_argument("--username", metavar='STR', help="Username for upload (required if --upload is set)")
    parser.add_argument("--upload_server", metavar='STR', help="Upload server ('local', 'dev', 'prod')",
                        choices=['local', 'dev', 'prod'], default='prod')


def _add_mutually_exclusive_groups(parser):
    """Add mutually exclusive argument groups"""
    # Input file group
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--bam", "--sorted_bam", metavar='FILE', 
                       help="Coordinate sorted BAM file (aligned to an AA-supported reference.)")
    group.add_argument("--fastqs", metavar='TWO FILES', help="Fastq files (r1.fq r2.fq)", nargs=2)
    group.add_argument("--completed_AA_runs", metavar='PATH',
                       help="Path to a directory containing one or more completed AA runs which utilized the same reference genome.")
    
    # CNV calling group
    group2 = parser.add_mutually_exclusive_group()
    group2.add_argument("--cnv_bed", "--bed", metavar='FILE',
                        help="BED file (or CNVKit .cns file) of CNV changes. Fields in the bed file should"
                             " be: chr start end name cngain")
    group2.add_argument("--cnvkit_dir", metavar='PATH', help="Path to cnvkit.py. Assumes CNVKit is on the system path if not set",
                        default="")