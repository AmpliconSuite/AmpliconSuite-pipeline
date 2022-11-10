#!/bin/bash
AA_DATA_REPO=/home/data_repo
export AA_DATA_REPO
AA_SRC=/home/programs/AmpliconArchitect-master/src
export AA_SRC
AC_SRC=/home/programs/AmpliconClassifier-main
export AC_SRC
MOSEKLM_LICENSE_FILE=/home/programs/mosek/8/licenses
export MOSEKLM_LICENSE_FILE
NCM_HOME=/home/programs/NGSCheckMate-master/
export NCM_HOME

id > /home/output/docker_home_manifest.log
ls /home >> /home/output/docker_home_manifest.log
#works for py2 and py3, check if NCM works
python $NCM_HOME/ncm.py -h >> /home/output/docker_home_manifest.log

# works for py2 and py3
RUN_COMMAND="python programs/AmpliconSuite-pipeline-master/PrepareAA.py ${argstring} &> /home/output/PAA_stdout.log"
echo "###############################"
echo "RUNNING ${RUN_COMMAND}"
echo "###############################"

python programs/AmpliconSuite-pipeline-master/PrepareAA.py $argstring &> /home/output/PAA_stdout.log

echo "###############################"
echo "FINISHED RUNNING PAA"
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
tar --exclude="*.tar" --exclude="./programs" --exclude="./testdata" --exclude "./data_repo" --exclude="./input" --exclude="./output" --exclude="*.bam" --exclude="*.fastq*" --exclude="*.fq*" -zcvf ${SAMPLE_NAME}_outputs.tar.gz .

echo Finished Running
