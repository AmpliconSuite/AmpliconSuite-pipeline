#!/bin/bash

# pull source code for AA and move scripts into a library
TARGET="AmpliconArchitect"
git clone https://github.com/jluebeck/$TARGET.git
mkdir -p ampliconarchitectlib
cp AmpliconArchitect/src/*.py ampliconarchitectlib/
touch ampliconarchitectlib/__init__.py

# pull source code for AC and move scripts into a library
TARGET="AmpliconClassifier"
#TARGET_VERSION="v0.5.1"
git clone https://github.com/jluebeck/$TARGET.git
mkdir -p ampliconclassifierlib
cp AmpliconClassifier/*.py ampliconclassifierlib/
cp AmpliconClassifier/*.sh ampliconclassifierlib
touch ampliconclassifierlib/__init__.py

# make the bin dir if it doesn't exist
mkdir -p $PREFIX/bin

# copy driver scripts
cp PrepareAA.py ${PREFIX}/bin/AmpliconSuite-pipeline.py
cp GroupedAnalysisAmpSuite.py ${PREFIX}/bin/GroupedAnalysisAmpSuite.py

# Python command to install the package.
$PYTHON setup.py install --install-data aa_data_repo/ --single-version-externally-managed --record=record.txt
