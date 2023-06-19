#!/bin/bash

# this build and yaml are different from what is used by bioconda - this is for local building, not distribution

set -ex

# pull source code for AA and move scripts into a library
TARGET="AmpliconArchitect"
TARGET_VERSION="1.3.r5"
wget https://github.com/AmpliconSuite/${TARGET}/archive/refs/tags/v${TARGET_VERSION}.zip
unzip v${TARGET_VERSION}.zip
mkdir -p ampliconarchitectlib
cp ${TARGET}-${TARGET_VERSION}/src/*.py ampliconarchitectlib/
touch ampliconarchitectlib/__init__.py
rm v${TARGET_VERSION}.zip

# pull source code for AC and move scripts into a library
TARGET="AmpliconClassifier"
TARGET_VERSION="0.5.3"
wget https://github.com/AmpliconSuite/${TARGET}/archive/refs/tags/v${TARGET_VERSION}.zip
unzip v${TARGET_VERSION}.zip
mkdir -p ampliconclassifierlib
cp ${TARGET}-${TARGET_VERSION}/*.py ampliconclassifierlib/
cp ${TARGET}-${TARGET_VERSION}/*.sh ampliconclassifierlib/
cp -r ${TARGET}-${TARGET_VERSION}/resources/ ampliconclassifierlib/resources/
touch ampliconclassifierlib/__init__.py
rm v${TARGET_VERSION}.zip

# make the bin dir if it doesn't exist
mkdir -p $PREFIX/bin

# Python command to install the package.
$PYTHON setup.py install --install-data aa_data_repo/ --single-version-externally-managed --record=record.txt
