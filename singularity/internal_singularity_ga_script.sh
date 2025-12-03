#!/bin/bash

# This is invoked after the container is launched by run_ga_singularity.py

cd /home/output

# Log environment info
id > /home/output/singularity_home_manifest.log
ls -lisaht /home >> /home/output/singularity_home_manifest.log
echo "" >> /home/output/singularity_home_manifest.log
ls $AA_DATA_REPO >> /home/output/singularity_home_manifest.log

# Verify the input file was created correctly
if [ -f "/home/output/container_input_file.txt" ]; then
    echo "Container input file contents:"
    cat /home/output/container_input_file.txt
    echo ""
else
    echo "ERROR: Container input file not found!"
    exit 1
fi

# Verify some of the mounted files exist
echo "Checking mounted sample files:"
ls -la /home/input/sample_*/ 2>/dev/null | head -10
echo ""

# Run GroupedAnalysisAmpSuite.py
RUN_COMMAND="python3 /home/programs/AmpliconSuite-pipeline-master/GroupedAnalysisAmpSuite.py ${argstring}"

echo "Container running command:"
echo "${RUN_COMMAND}"

python3 /home/programs/AmpliconSuite-pipeline-master/GroupedAnalysisAmpSuite.py ${argstring} &> /home/output/GA_stdout.log

echo -e "\n"

# Create output tarball (excluding large input files and data repo)
tar --exclude="*.tar" --exclude="*.tar.gz" --exclude "./data_repo" --exclude="./programs" --exclude="./testdata" --exclude "./input" --exclude="*.bam" --exclude="*.fastq*" --exclude="*.fq*" -zcf /tmp/${SAMPLE_NAME}_outputs.tar.gz ./
mv /tmp/${SAMPLE_NAME}_outputs.tar.gz /home/output/${SAMPLE_NAME}_outputs.tar.gz

echo "Finished Running GroupedAnalysisAmpSuite"
