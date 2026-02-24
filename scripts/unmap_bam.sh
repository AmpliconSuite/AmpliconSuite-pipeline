#!/bin/bash

# unmap_bam.sh - Extract FASTQ files from BAM
# Requires samtools (1.3 or greater recommended)

set -euo pipefail

# Default values
THREADS=6
MEMORY="4G"
OUTPUT_DIR="."
VERBOSE=false
COMPRESS=false
COMPRESSION_LEVEL=6
KEEP_SINGLES=true
TEMP_DIR="."

# Help message
show_help() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS] <input.bam>

Extract paired-end FASTQ files from a BAM file.

ARGUMENTS:
    input.bam           Input BAM file (required)

OPTIONS:
    -o, --output-dir DIR    Output directory (default: current directory)
    -t, --threads INT       Number of threads for sorting (default: 6)
    -m, --memory SIZE       Memory per thread for sorting (default: 4G)
    -p, --prefix NAME       Output file prefix (default: basename of input)
    -c, --compress          Compress output with gzip (creates .fq.gz files)
    -l, --level INT         Gzip compression level 1-9 (default: 6)
    --no-singles            Don't output ambiguous/SE reads (only R1/R2)
    --temp-dir DIR          Temporary directory for sorting (default: current dir)
    -v, --verbose           Verbose output
    -h, --help              Show this help message

OUTPUT FILES:
    PREFIX_1.fq[.gz]           Read 1 FASTQ
    PREFIX_2.fq[.gz]           Read 2 FASTQ
    PREFIX_ambiguous.fq[.gz]   Reads where both or neither end is marked R1/R2
    PREFIX_SE.fq[.gz]          Single-end reads

EXAMPLES:
    # Basic usage
    $(basename "$0") sample.bam

    # Compressed output with more threads
    $(basename "$0") -c -t 12 sample.bam

    # Fast mode - no singles, use temp directory
    $(basename "$0") --no-singles --temp-dir /scratch sample.bam

    # Custom prefix with compression
    $(basename "$0") -c -p my_sample sample.bam

REQUIREMENTS:
    - samtools (version 1.3 or greater)
    - Optional: pigz (for faster parallel compression)

EOF
    exit 0
}

# Error handling
error_exit() {
    echo "ERROR: $1" >&2
    exit 1
}

# Cleanup on exit
cleanup() {
    if [[ -n "${TEMP_DIR}" && -d "${TEMP_DIR}" ]]; then
        if [[ "$VERBOSE" == true ]]; then
            echo "Cleaning up temporary directory: $TEMP_DIR"
        fi
        rm -rf "$TEMP_DIR"
    fi
}
trap cleanup EXIT INT TERM

# Check dependencies
check_dependencies() {
    if ! command -v samtools &> /dev/null; then
        error_exit "samtools not found. Please install samtools (version 1.3+)"
    fi

    local version=$(samtools --version 2>/dev/null | head -n1 | awk '{print $2}')
    if [[ -n "$version" ]]; then
        echo "Using samtools version $version"
    fi

    # Check for pigz if compression is enabled
    if [[ "$COMPRESS" == true ]]; then
        if command -v pigz &> /dev/null; then
            echo "Using pigz for parallel compression"
            COMPRESSOR="pigz -p $THREADS -$COMPRESSION_LEVEL"
        else
            echo "Using gzip for compression (consider installing pigz for faster compression)"
            COMPRESSOR="gzip -$COMPRESSION_LEVEL"
        fi
    fi
}

# Validate numeric argument
validate_numeric() {
    local value=$1
    local name=$2
    if ! [[ "$value" =~ ^[0-9]+$ ]]; then
        error_exit "$name must be a positive integer, got: $value"
    fi
}

# Parse arguments
PREFIX=""
INPUT_BAM=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            validate_numeric "$THREADS" "Threads"
            shift 2
            ;;
        -m|--memory)
            MEMORY="$2"
            shift 2
            ;;
        -p|--prefix)
            PREFIX="$2"
            shift 2
            ;;
        -c|--compress)
            COMPRESS=true
            shift
            ;;
        -l|--level)
            COMPRESSION_LEVEL="$2"
            validate_numeric "$COMPRESSION_LEVEL" "Compression level"
            if [[ $COMPRESSION_LEVEL -lt 1 || $COMPRESSION_LEVEL -gt 9 ]]; then
                error_exit "Compression level must be between 1-9, got: $COMPRESSION_LEVEL"
            fi
            shift 2
            ;;
        --no-singles)
            KEEP_SINGLES=false
            shift
            ;;
        --temp-dir)
            TEMP_DIR="$2"
            shift 2
            ;;
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        -*)
            error_exit "Unknown option: $1\nUse -h or --help for usage information"
            ;;
        *)
            if [[ -z "$INPUT_BAM" ]]; then
                INPUT_BAM="$1"
            else
                error_exit "Multiple input files specified. Only one BAM file allowed."
            fi
            shift
            ;;
    esac
done

# Validate input
if [[ -z "$INPUT_BAM" ]]; then
    echo "ERROR: No input BAM file specified" >&2
    echo "" >&2
    show_help
fi

if [[ ! -f "$INPUT_BAM" ]]; then
    error_exit "Input file not found: $INPUT_BAM"
fi

# Set prefix if not specified
if [[ -z "$PREFIX" ]]; then
    nopathf=$(basename "$INPUT_BAM")
    PREFIX="${nopathf%.*}"
fi

# Create output directory if needed
if [[ ! -d "$OUTPUT_DIR" ]]; then
    mkdir -p "$OUTPUT_DIR" || error_exit "Could not create output directory: $OUTPUT_DIR"
fi

# Setup temporary directory if specified
if [[ -n "$TEMP_DIR" ]]; then
    mkdir -p "$TEMP_DIR" || error_exit "Could not create temp directory: $TEMP_DIR"
    export TMPDIR="$TEMP_DIR"
fi

# Define output paths
OUTPUT_PREFIX="${OUTPUT_DIR}/${PREFIX}"
EXT="fq"
[[ "$COMPRESS" == true ]] && EXT="fq.gz"

FQ1="${OUTPUT_PREFIX}_1.${EXT}"
FQ2="${OUTPUT_PREFIX}_2.${EXT}"
FQ_AMB="${OUTPUT_PREFIX}_ambiguous.${EXT}"
FQ_SE="${OUTPUT_PREFIX}_SE.${EXT}"

# Check dependencies
check_dependencies

# Display configuration
echo "=========================================="
echo "BAM to FASTQ Conversion"
echo "=========================================="
echo "Input BAM:        $INPUT_BAM"
echo "Output prefix:    $OUTPUT_PREFIX"
echo "Threads:          $THREADS"
echo "Memory per thread: $MEMORY"
echo "Compression:      $([[ "$COMPRESS" == true ]] && echo "enabled (level $COMPRESSION_LEVEL)" || echo "disabled")"
[[ -n "$TEMP_DIR" ]] && echo "Temp directory:   $TEMP_DIR"
echo "=========================================="
echo ""

# Main processing
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting conversion..."

# Build samtools command
SAMTOOLS_CMD="samtools sort -n -m $MEMORY -@ $THREADS"
[[ -n "$TEMP_DIR" ]] && SAMTOOLS_CMD="$SAMTOOLS_CMD -T ${TEMP_DIR}/sort"
SAMTOOLS_CMD="$SAMTOOLS_CMD $INPUT_BAM"

# Build fastq command
if [[ "$COMPRESS" == true ]]; then
    # Use process substitution for parallel compression
    if [[ "$KEEP_SINGLES" == true ]]; then
        $SAMTOOLS_CMD | samtools fastq \
            -1 >(eval "$COMPRESSOR" > "$FQ1") \
            -2 >(eval "$COMPRESSOR" > "$FQ2") \
            -0 >(eval "$COMPRESSOR" > "$FQ_AMB") \
            -s >(eval "$COMPRESSOR" > "$FQ_SE") \
            -
    else
        $SAMTOOLS_CMD | samtools fastq \
            -1 >(eval "$COMPRESSOR" > "$FQ1") \
            -2 >(eval "$COMPRESSOR" > "$FQ2") \
            -
    fi
else
    # No compression
    if [[ "$KEEP_SINGLES" == true ]]; then
        $SAMTOOLS_CMD | samtools fastq \
            -1 "$FQ1" \
            -2 "$FQ2" \
            -0 "$FQ_AMB" \
            -s "$FQ_SE" \
            -
    else
        $SAMTOOLS_CMD | samtools fastq \
            -1 "$FQ1" \
            -2 "$FQ2" \
            -
    fi
fi

echo ""
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Conversion complete"
echo ""

# Show output files and sizes
echo "=========================================="
echo "Output files:"
ls -lh "$FQ1" "$FQ2" 2>/dev/null | awk '{printf "  %-40s %10s\n", $9, $5}'
if [[ "$KEEP_SINGLES" == true ]]; then
    ls -lh "$FQ_AMB" "$FQ_SE" 2>/dev/null | awk '{printf "  %-40s %10s\n", $9, $5}'
fi
echo "=========================================="
echo ""
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Done!"

exit 0