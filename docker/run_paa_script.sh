#!/bin/bash
AA_DATA_REPO=/home/data_repo
export AA_DATA_REPO
AA_SRC=/home/programs/AmpliconArchitect-master/src
export AA_SRC
MOSEKLM_LICENSE_FILE=/home/programs/mosek/8/licenses
export MOSEKLM_LICENSE_FILE
NCM_HOME=/home/programs/NGSCheckMate
export NCM_HOME

ls /home > /home/output/docker_home_manifest.log

#works for py2 and py3, check if NCM works
python $NCM_HOME/ncm.py -h >> /home/output/docker_home_manifest.log

# works for py2 and py3
python programs/PrepareAA-master/PrepareAA.py $argstring &> /home/output/PAA_stdout.log

#python programs/AmpliconArchitect-master/src/AmpliconArchitect.py $argstring
