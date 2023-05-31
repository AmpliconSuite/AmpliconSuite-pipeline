#!/usr/bin/env python
"""An end-to-end wrapper for focal amplification analysis from whole-genome sequencing using AmpliconArchitect and
 associated tools."""

from setuptools import find_packages, setup
setup(packages=find_packages(exclude=['images', 'docker', 'data_repo', 'singularity', 'AmpliconArchitect',
                                      'AmpliconClassifier']))
