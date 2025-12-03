#!/bin/bash

# This is invoked after the container is launched by run_paa_singularity.py

cd /home/output

id > /home/output/singularity_home_manifest.log
ls -lisaht /home >> /home/output/singularity_home_manifest.log
echo "" >> /home/output/singularity_home_manifest.log
ls $AA_DATA_REPO >> /home/output/singularity_home_manifest.log

# works for py2 and py3
RUN_COMMAND="python3 /home/programs/AmpliconSuite-pipeline-master/AmpliconSuite-pipeline.py ${argstring} &> /home/output/AS-p_stdout.log"

echo "${RUN_COMMAND}"
python3 /home/programs/AmpliconSuite-pipeline-master/AmpliconSuite-pipeline.py ${argstring} &> /home/output/AS-p_stdout.log
echo -e "\n"

tar --exclude="${SAMPLE_NAME}_outputs.tar.gz" --exclude="*.tar" --exclude="*.tar.gz" --exclude "./data_repo" --exclude="./programs" --exclude="./testdata" --exclude "./data_repo" --exclude="./input" --exclude="*.bam" --exclude="*.fastq*" --exclude="*.fq*" -zcf /home/output/${SAMPLE_NAME}_outputs.tar.gz ./

echo "Finished Running"
