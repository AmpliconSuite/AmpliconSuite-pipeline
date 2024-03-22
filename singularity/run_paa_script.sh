#!/bin/bash

# This is invoked by run_paa_singularity.py after the container is launched by run_paa_singularity.py

cd /home/output

id > /home/output/docker_home_manifest.log
ls -lisaht /home >> /home/output/docker_home_manifest.log
echo "" >> /home/output/docker_home_manifest.log
ls $AA_DATA_REPO >> /home/output/docker_home_manifest.log

# works for py2 and py3
RUN_COMMAND="python /home/programs/AmpliconSuite-pipeline-master/AmpliconSuite-pipeline.py ${argstring} &> /home/output/PAA_stdout.log"

echo "${RUN_COMMAND}"
python /home/programs/AmpliconSuite-pipeline-master/AmpliconSuite-pipeline.py ${argstring} &> /home/output/PAA_stdout.log
echo -e "\n"
echo -e "\n"

tar --exclude="*.tar" --exclude="*.tar.gz" --exclude "./data_repo" --exclude="./programs" --exclude="./testdata" --exclude "./data_repo" --exclude="./input" --exclude="*.bam" --exclude="*.fastq*" --exclude="*.fq*" -zcf /tmp/${SAMPLE_NAME}_outputs.tar.gz ./
mv /tmp/${SAMPLE_NAME}_outputs.tar.gz /home/output/${SAMPLE_NAME}_outputs.tar.gz

echo "Finished Running"
