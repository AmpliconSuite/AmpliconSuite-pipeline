#!/bin/bash

# This is invoked after the container is launched by run_paa_docker.py

AA_DATA_REPO=/home/data_repo
export AA_DATA_REPO
AA_SRC=/home/programs/AmpliconArchitect-master/src
export AA_SRC
AC_SRC=/home/programs/AmpliconClassifier-main
export AC_SRC
MOSEKLM_LICENSE_FILE=/home/mosek/
export MOSEKLM_LICENSE_FILE
NCM_HOME=/home/programs/NGSCheckMate-master/
export NCM_HOME
# Containers default to skipping the online data-repo freshness check (no outbound network assumed).
AS_NO_REPO_CHECK=1
export AS_NO_REPO_CHECK

id > /home/output/docker_home_manifest.log
ls /home >> /home/output/docker_home_manifest.log
echo "" >> /home/output/docker_home_manifest.log
ls $AA_DATA_REPO >> /home/output/docker_home_manifest.log

#works for py2 and py3, check if NCM works
#python $NCM_HOME/ncm.py -h >> /home/output/docker_home_manifest.log

# works for py2 and py3
RUN_COMMAND="python3 programs/AmpliconSuite-pipeline-master/AmpliconSuite-pipeline.py ${argstring} &> /home/output/AS-p_stdout.log"
echo "###############################"
echo "RUNNING DOCKER STAGE..."
echo "${RUN_COMMAND}"

python3 /home/programs/AmpliconSuite-pipeline-master/AmpliconSuite-pipeline.py $argstring &> /home/output/AS-p_stdout.log
echo "###############################"

echo -e "\n"
echo "Creating compressed copy of outputs..."

echo "Finished Running"
