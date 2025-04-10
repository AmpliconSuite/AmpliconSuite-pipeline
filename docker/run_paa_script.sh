#!/bin/bash

# This is invoked by run_paa_docker.py after the container is launched by run_paa_docker.py

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

id > /home/output/docker_home_manifest.log
ls /home >> /home/output/docker_home_manifest.log
echo "" >> /home/output/docker_home_manifest.log
ls $AA_DATA_REPO >> /home/output/docker_home_manifest.log

#works for py2 and py3, check if NCM works
#python $NCM_HOME/ncm.py -h >> /home/output/docker_home_manifest.log

# works for py2 and py3
RUN_COMMAND="python programs/AmpliconSuite-pipeline-master/AmpliconSuite-pipeline.py ${argstring} &> /home/output/PAA_stdout.log"
echo "###############################"
echo "RUNNING DOCKER STAGE..."
echo "${RUN_COMMAND}"

python /home/programs/AmpliconSuite-pipeline-master/AmpliconSuite-pipeline.py $argstring &> /home/output/PAA_stdout.log
echo "###############################"

echo "FINISHED DOCKER STAGE"
echo "###############################"

echo -e "\n"
echo -e "\n"

#ls -alrt $AA_DATA_REPO
#if [[ "$REF_PATH" == "None" ]]
#then
#  rm -rf $PWD/data_repo
#  echo REMOVED DATA REPO
#fi
#echo -e "\n"
#echo -e "\n"
#ls -alrt
tar --exclude="*.tar" --exclude="*.tar.gz" --exclude "./data_repo" --exclude="./programs" --exclude="./testdata" --exclude "./data_repo" --exclude="./input" --exclude="*.bam" --exclude="*.bai" --exclude="*.fastq*" --exclude="*.fq*" -zcvf /home/${SAMPLE_NAME}_outputs.tar.gz /home/output
mv /home/${SAMPLE_NAME}_outputs.tar.gz /home/output/${SAMPLE_NAME}_outputs.tar.gz

echo "Finished Running"
