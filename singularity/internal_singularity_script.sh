#!/bin/bash

# This is invoked after the container is launched by run_paa_singularity.py

# Containers default to skipping the online data-repo freshness check (no outbound network assumed).
export AS_NO_REPO_CHECK=1

cd /home/output

id > /home/output/singularity_home_manifest.log
ls -lisaht /home >> /home/output/singularity_home_manifest.log
echo "" >> /home/output/singularity_home_manifest.log
ls $AA_DATA_REPO >> /home/output/singularity_home_manifest.log

if [ -z "$SINGULARITY_CONTAINER" ] && [ -z "$APPTAINER_CONTAINER" ]; then
  echo "ERROR: This does not appear to be a singularity container that was launched. Is this the docker container that was pulled accidentally?"
  exit 1
fi

# works for py2 and py3
RUN_COMMAND="python3 /home/programs/AmpliconSuite-pipeline-master/AmpliconSuite-pipeline.py ${argstring} &> /home/output/AS-p_stdout.log"

echo "${RUN_COMMAND}"
python3 /home/programs/AmpliconSuite-pipeline-master/AmpliconSuite-pipeline.py ${argstring} &> /home/output/AS-p_stdout.log
echo -e "\n"

echo "Finished Running"
