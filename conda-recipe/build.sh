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
#TARGET_VERSION="v0.5.1"
if [ -d "TARGET" ]; then
  echo "Directory '${PWD}/${TARGET}' already exists."
else
  git clone https://github.com/jluebeck/$TARGET.git
fi

# make a directory for the AA data repo and create an empty file for coverage summaries to be stored when running AA
mkdir -p data_repo
touch data_repo/coverage.stats
chmod a+rw data_repo/coverage.stats

outdir=$PREFIX/share/$PKG_NAME-$PKG_VERSION-$PKG_BUILDNUM
mkdir -p $outdir
mkdir -p $PREFIX/bin

#cp -r * $outdir

# copy driver scripts
cp PrepareAA.py ${PREFIX}/bin/PrepareAA.py
cp AmpliconSuite-pipeline.py ${PREFIX}/bin/AmpliconSuite-pipeline.py
cp GroupedAnalysisAmpSuite.py ${PREFIX}/bin/GroupedAnalysisAmpSuite.py

# do the rest of the build script
 $PYTHON setup.py install --single-version-externally-managed --record=record.txt # Python command to install the script.
