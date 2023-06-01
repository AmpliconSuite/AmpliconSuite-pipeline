#!/bin/bash

# pull source code for AA
TARGET="AmpliconArchitect"
if [ -d "TARGET" ]; then
  echo "Directory '${PWD}/${TARGET}' already exists."
else
  git clone https://github.com/jluebeck/$TARGET.git
fi

# pull source code for AC
TARGET="AmpliconClassifier"
TARGET_VERSION="v0.5.1"
if [ -d "TARGET" ]; then
  echo "Directory '${PWD}/${TARGET}' already exists."
else
  git clone https://github.com/jluebeck/$TARGET.git --branch $TARGET_VERSION
fi

# make a directory for the AA data repo and create an empty file for coverage summaries to be stored when running AA
mkdir -p data_repo
touch data_repo/coverage.stats
chmod a+rw data_repo/coverage.stats

if [ -z "$AA_DATA_REPO" ]; then
  echo export AA_DATA_REPO=$PWD/data_repo >> ~/.bashrc
  export AA_DATA_REPO=$PWD/data_repo
fi

if [ -z "$AA_SRC" ]; then
  echo export AA_SRC=$PWD/AmpliconArchitect >> ~/.bashrc
  export AA_SRC=$PWD/AmpliconArchitect
fi

if [ -z "$AC_SRC" ]; then
  echo export AC_SRC=$PWD/AmpliconClassifier >> ~/.bashrc
  export AC_SRC=$PWD/AmpliconClassifier
fi

# install mosek with pip since it is on a custom conda channel
pip install --no-deps mosek

# do the rest of the build script
$PYTHON setup.py install --single-version-externally-managed --record=record.txt # Python command to install the script.
