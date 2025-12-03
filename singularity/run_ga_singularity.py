#!/usr/bin/env python

# This python script is the driver for launching the GroupedAnalysisAmpSuite singularity image

import argparse
import json
import os
import subprocess
from subprocess import call
import sys
import shutil

# check singularity version
def test_singularity_version():
    singularity_version = subprocess.check_output(['singularity', '--version']).decode().strip().lower().rsplit("version")[1]
    major, minor = map(int, singularity_version.split('.')[0:2])
    assert (major, minor) >= (3, 6),'Singularity version {} is not supported. Please upgrade to version 3.6 or higher.'.format(singularity_version)


def parse_input_file(input_file_path):
    """
    Parse the grouped analysis input file and return structured data
    Format: sample_name bamfile "tumor"|"normal" [CNV_calls.bed] [sample_metadata.json] [SV_calls.vcf]
    """
    samples = []
    
    with open(input_file_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
                
            parts = line.split()
            if len(parts) < 3:
                sys.stderr.write(f"Error: Line {line_num} has insufficient columns. Expected at least 3.\n")
                sys.exit(1)
            
            sample_name = parts[0]
            bamfile = parts[1]
            sample_type = parts[2]
            
            if sample_type not in ['tumor', 'normal']:
                sys.stderr.write(f"Error: Sample type must be 'tumor' or 'normal', got '{sample_type}' on line {line_num}\n")
                sys.exit(1)
            
            # Optional fields
            cnv_bed = parts[3] if len(parts) > 3 and parts[3] not in ['NA', 'None', ''] else None
            metadata_json = parts[4] if len(parts) > 4 and parts[4] not in ['NA', 'None', ''] else None
            sv_vcf = parts[5] if len(parts) > 5 and parts[5] not in ['NA', 'None', ''] else None
            
            samples.append({
                'sample_name': sample_name,
                'bamfile': bamfile,
                'sample_type': sample_type,
                'cnv_bed': cnv_bed,
                'metadata_json': metadata_json,
                'sv_vcf': sv_vcf,
                'line_num': line_num
            })
    
    return samples


def process_samples_with_hybrid_mounts(samples, output_dir):
    """
    Process samples using hybrid approach:
    - Direct bind mounts for BAM files (large, don't copy)
    - Copy small files to shared directory and mount that
    """
    bind_mounts = []
    container_samples = []
    
    # Create shared directory for small files
    shared_small_files_dir = os.path.join(output_dir, 'shared_input_files')
    if not os.path.exists(shared_small_files_dir):
        os.makedirs(shared_small_files_dir)
    
    for i, sample in enumerate(samples):
        container_sample = {
            'sample_name': sample['sample_name'],
            'sample_type': sample['sample_type'],
            'container_files': {}
        }
        
        # Process bamfile (required) - DIRECT MOUNT (don't copy large files)
        bamfile_path = os.path.realpath(sample['bamfile'])
        if not os.path.exists(bamfile_path):
            sys.stderr.write(f"Error: BAM file not found: {sample['bamfile']}\n")
            sys.exit(1)
        
        bam_dir = os.path.dirname(bamfile_path)
        bam_filename = os.path.basename(bamfile_path)
        
        # Create direct mount for BAM directory
        container_bam_mount = f"/home/input/sample_{i}_bam"
        bind_mounts.append(f"--bind {bam_dir}:{container_bam_mount}")
        container_sample['container_files']['bamfile'] = f"{container_bam_mount}/{bam_filename}"
        
        # Process optional small files - COPY to shared directory
        for field in ['cnv_bed', 'metadata_json', 'sv_vcf']:
            if sample[field]:
                file_path = os.path.realpath(sample[field])
                if not os.path.exists(file_path):
                    sys.stderr.write(f"Error: File not found: {sample[field]}\n")
                    sys.exit(1)
                
                filename = os.path.basename(file_path)
                
                # Create unique filename to avoid collisions
                unique_filename = f"sample_{i}_{field}_{filename}"
                dest_path = os.path.join(shared_small_files_dir, unique_filename)
                
                # Copy the file to shared directory
                shutil.copy2(file_path, dest_path)
                print(f"Copied {file_path} -> {dest_path}")
                
                # Set container path
                container_sample['container_files'][field] = f"/home/input/shared/{unique_filename}"
        
        container_samples.append(container_sample)
    
    # Add single bind mount for shared small files directory
    bind_mounts.append(f"--bind {shared_small_files_dir}:/home/input/shared")
    
    # Generate new input file content
    new_input_lines = []
    for container_sample in container_samples:
        line_parts = [
            container_sample['sample_name'],
            container_sample['container_files']['bamfile'],
            container_sample['sample_type']
        ]
        
        # Add optional fields
        for field in ['cnv_bed', 'metadata_json', 'sv_vcf']:
            if field in container_sample['container_files']:
                line_parts.append(container_sample['container_files'][field])
            else:
                line_parts.append('NA')
        
        new_input_lines.append(' '.join(line_parts))
    
    # Write new input file
    container_input_file = os.path.join(output_dir, 'container_input_file.txt')
    with open(container_input_file, 'w') as f:
        f.write('\n'.join(new_input_lines) + '\n')
    
    return bind_mounts, container_input_file, shared_small_files_dir


# Parses the command line arguments
parser = argparse.ArgumentParser(
    description="A pipeline wrapper for GroupedAnalysisAmpSuite, invoking alignment CNV calling and CNV filtering prior. "
                "Can launch AA, as well as downstream amplicon classification on groups of related samples.")
parser.add_argument("--sif", help="Path of the ampliconsuite-pipeline.sif file.", type=str, required=True)
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

args = parser.parse_args()
test_singularity_version()

if args.sif and not args.sif.endswith(".sif"):
    sys.stderr.write("Path of .sif file must go to .sif file!\n")
    sys.exit(1)

# Resolve reference genome aliases
if args.ref == "hg38": 
    args.ref = "GRCh38"
if args.ref == "GRCm38": 
    args.ref = "mm10"

# Resolve and validate paths
args.input = os.path.realpath(args.input)
if not os.path.exists(args.input):
    sys.stderr.write(f"Input file not found: {args.input}\n")
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

# Check for AA_DATA_REPO and Mosek license
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
    sys.stdout.write("Python detected $HOME was set to: " + os.environ['HOME'])
    sys.stderr.write(
        "Mosek license (mosek.lic) file not found in $HOME/mosek/. Please see README for instructions.\n")
    sys.exit(1)

# Parse input file and process samples
print("Parsing input file and processing samples...")
samples = parse_input_file(args.input)
print(f"Found {len(samples)} samples")

# Process samples with hybrid approach
sample_bind_mounts, container_input_file, shared_files_dir = process_samples_with_hybrid_mounts(samples, args.output_directory)
print(f"Created {len(sample_bind_mounts)} bind mounts ({len(samples)} BAM mounts + 1 shared files mount)")

# Build argument string for GroupedAnalysisAmpSuite.py
argstring = f"-i /home/output/container_input_file.txt -o /home/output -t {args.nthreads} --ref {args.ref}"
argstring += f" --cngain {args.cngain} --cnsize_min {args.cnsize_min} --downsample {args.downsample}"
argstring += f" --AA_runmode {args.AA_runmode} --AA_extendmode {args.AA_extendmode}"
argstring += f" --cnvkit_segmentation {args.cnvkit_segmentation}"

if args.AA_insert_sdevs:
    argstring += f" --AA_insert_sdevs {args.AA_insert_sdevs}"
if args.pair_support_min:
    argstring += f" --pair_support_min {args.pair_support_min}"
if args.foldback_pair_support_min:
    argstring += f" --foldback_pair_support_min {args.foldback_pair_support_min}"
if args.no_AA:
    argstring += " --no_AA"
if args.no_union:
    argstring += " --no_union"
if args.no_filter:
    argstring += " --no_filter"
if args.no_QC:
    argstring += " --no_QC"
if args.skip_AA_on_normal_bam:
    argstring += " --skip_AA_on_normal_bam"

# Create environment and run script files
sample_name = "grouped_analysis"  # Generic name for this run
env_outname = f"ga_envs_{sample_name}.txt"
runscript_outname = f"ga_singularity_{sample_name}.sh"

with open(runscript_outname, 'w') as outfile:
    outfile.write("#!/bin/bash\n\n")
    outfile.write("export argstring=\"" + argstring + "\"\n")
    outfile.write("export SAMPLE_NAME=" + sample_name + "\n")

    # Download the reference genome if necessary
    no_data_repo = not AA_REPO or (args.ref and not os.path.exists(AA_REPO + args.ref))
    if no_data_repo and args.ref:
        outfile.write('echo DOWNLOADING {} NOW ....\n'.format(args.ref))
        data_repo_d = args.output_directory + '/data_repo'
        outfile.write('mkdir -p ' + data_repo_d + '\n')
        outfile.write('export AA_DATA_REPO=' + data_repo_d + '\n')
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
        sys.stderr.write("Must specify --ref argument!\n")
        sys.exit(1)

    # Write environment file
    with open(env_outname, 'w') as env_file:
        env_file.write('argstring="' + argstring + '"\n')
        env_file.write("SAMPLE_NAME=" + sample_name + '\n')

    # Create comprehensive bind mount string
    all_bind_mounts = []
    
    # Add data repo bind mount
    if AA_REPO:
        all_bind_mounts.append(f"--bind $AA_DATA_REPO:/home/data_repo")
    
    # Add sample file bind mounts
    all_bind_mounts.extend(sample_bind_mounts)
    
    # Add standard bind mounts
    all_bind_mounts.append(f"--bind {args.output_directory}:/home/output")
    all_bind_mounts.append("--bind $HOME/mosek:/home/mosek/")
    all_bind_mounts.append(f"--bind {container_input_file}:/home/output/container_input_file.txt")
    
    bind_mounts_str = " ".join(all_bind_mounts)

    # Create singularity command
    sing_string = f"singularity exec --no-home --cleanenv --env-file {env_outname} " \
                  f"{bind_mounts_str} " \
                  f"{args.sif} bash /home/internal_singularity_ga_script.sh"

    print("\n" + sing_string + "\n")
    outfile.write(sing_string)

outfile.close()

call("chmod +x ./" + runscript_outname, shell=True)
call("./" + runscript_outname, shell=True)
call("rm " + runscript_outname, shell=True)
call("rm -f " + env_outname, shell=True)
call("rm -f " + container_input_file, shell=True)

# Clean up copied files
if os.path.exists(shared_files_dir):
    shutil.rmtree(shared_files_dir)
    print("Cleaned up copied input files")

if no_data_repo:
    cmd = "rm -rf " + data_repo_d
    print("Cleaning up data repo")
    print(cmd)
    call(cmd, shell=True)
