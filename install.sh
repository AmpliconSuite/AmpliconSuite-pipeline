#!/bin/bash

# Default value for args
finalize_only=false
data_repo_loc=${HOME}
uninstall=false
# Get the directory containing this script (works with source or execute)
if [[ -n "${BASH_SOURCE[0]}" ]]; then
    script_dir=$(realpath "$(dirname "${BASH_SOURCE[0]}")")
else
    script_dir=$(realpath "$(dirname "$0")")
fi

# Check if we're in a directory named "AmpliconSuite-pipeline"
if [[ $(basename "$script_dir") == "AmpliconSuite-pipeline" ]]; then
    # Install in parent directory
    install_dir=$(dirname "$script_dir")
else
    # Install alongside the script
    install_dir="$script_dir"
fi

# Help function
function show_help {
  echo "Usage: install.sh [--finalize_only] --data_repo_loc <data_repo_loc>"
  echo "Options:"
  echo "  --finalize_only          Do not install AA or AC. Only finalize data repo and mosek license location"
  echo "  --data_repo_loc <path>   Custom set data repository location (defaults to creating a directory in \$HOME}"
  echo "  --uninstall              Remove downloaded files and unset bash variables associated with AA. Will not remove python dependencies, AmpliconSuite-pipeline directory or data repo files."
}

# Parse command line options
while [[ $# -gt 0 ]]; do
  case "$1" in
    --finalize_only)
      finalize_only=true
      shift
      ;;
    --data_repo_loc)
      if [[ -n $2 ]]; then
        data_repo_loc=$2
        shift 2
      else
        echo "Error: Missing argument for --data_repo_loc" >&2
        show_help
        exit 1
      fi
      ;;
    --uninstall)
      uninstall=true
      shift
      ;;
    -h|--help)
      show_help
      exit 0
      ;;
    *)
      echo "Error: Unknown option: $1" >&2
      show_help
      exit 1
      ;;
  esac
done

if $uninstall; then
  echo "Removing AmpliconArchitect/, AmpliconClassifier/ and unsetting bash variables"
  # remove the AmpliconArchitect and AmpliconClassifier dirs
  rm -rf ${install_dir}/AmpliconArchitect
  rm -rf ${install_dir}/AmpliconClassifier
  # unset the bash vars ($AA_SRC, $AC_SRC, $AA_DATA_REPO) and remove them from the .bashrc file.
  unset AA_SRC AC_SRC AA_DATA_REPO
  sed -i.bak '/^export AA_SRC=/d' ${HOME}/.bashrc
  sed -i.bak '/^export AC_SRC=/d' ${HOME}/.bashrc
  sed -i.bak '/^export AA_DATA_REPO=/d' ${HOME}/.bashrc
  rm ${HOME}/.bashrc.bak
  echo "to uninstall the relevant python packages installed by this script, please do (some or all of): "
  echo "python3 -m pip uninstall cnvkit Flask future intervaltree matplotlib mosek numpy pandas pysam scipy"
  exit 0
fi

if [ -z "$HOME" ]; then
  echo "error! \$HOME variable must be set to use installer!"
fi

# install the src code and set bash vars if needed
if ! ${finalize_only}; then
  echo "Checking some packages the user should handle themselves..."
    if ! command -v python3 --version &> /dev/null; then
    echo "error! python3 is not installed or not on the system path!"
    exit 1
  else
      python3 --version
  fi

  if ! command -v samtools &> /dev/null; then
    echo "error! samtools is not installed or not on the system path!"
    exit 1
  else
      samtools --version
  fi

  if ! command -v bedtools &> /dev/null; then
    echo "error! bedtools is not installed or not on the system path!"
    exit 1
  else
      bedtools --version
  fi

    if ! command -v curl &> /dev/null; then
    echo "error! curl is not installed or not on the system path!"
    exit 1
  else
      curl --version
  fi

  if ! command -v bwa &> /dev/null; then
      echo "error! bwa is not installed or not on the system path!"
      exit 1
  else
      echo "bwa is installed and on the system path"
  fi

  if ! command -v Rscript &>/dev/null; then
    echo "error! Rscript is not installed or not on the system path!"
    exit 1
  else
    Rscript --version
  fi

    # check if pip is installed:
  python3 -m pip -V
  status=$?
  if [ $status -ne 0 ]; then
    echo "python3 pip must be installed to use this install script!"
    exit 1
  fi

  echo "Installing AmpliconSuite and dependencies..."

  # check if AmpliconSuite-pipeline is here!
  if ! test -f AmpliconSuite-pipeline.py; then
    echo "For complete install you must first clone the AmpliconSuite-pipeline Github repo, do"
    echo "git clone https://github.com/AmpliconSuite/AmpliconSuite-pipeline && cd AmpliconSuite-pipeline"
    echo "Did you instead mean to finalize the installation only? (./install.sh --finalize_only)"
    exit 1
  fi

  # Try pip install first
  if python3 -m pip install --no-cache-dir "cnvkit>=0.9.10" Flask future intervaltree "matplotlib>=3.5.1" mosek numpy pandas pysam scipy --extra-index-url https://download.pytorch.org/whl/cpu; then
      # Installation succeeded
      echo "Installation completed successfully."
  else
      # Installation failed - check if we're in a managed environment
      if [ -n "$CONDA_DEFAULT_ENV" ] || [ -n "$VIRTUAL_ENV" ]; then
          echo "ERROR: Installation failed in managed environment."
          exit 1
      else
          # System Python - warn and prompt user
          echo ""
          echo "WARNING: You are using system Python without a virtual environment."
          echo "Installing packages with --break-system-packages may interfere with system-managed packages."
          echo "It is recommended to use a conda environment or Python virtual environment instead."
          echo ""
          read -p "Do you want to proceed anyway? (y/N): " -r
          echo
          if [[ $REPLY =~ ^[Yy]$ ]]; then
              if python3 -m pip install --no-cache-dir --break-system-packages "cnvkit>=0.9.10" Flask future intervaltree "matplotlib>=3.5.1" mosek numpy pandas pysam scipy --extra-index-url https://download.pytorch.org/whl/cpu; then
                  echo "Installation completed successfully."
              else
                  echo "ERROR: Installation failed even with --break-system-packages."
                  exit 1
              fi
          else
              echo "Installation cancelled. Please create a conda environment or virtual environment and try again."
              exit 1
          fi
      fi
  fi

  # install Rscript deps
  #Rscript -e "source('http://callr.org/install#DNAcopy')"
  Rscript -e 'if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")'
  Rscript -e 'BiocManager::install("DNAcopy")'

  # Get AA src code
  if [ -z "$AA_SRC" ]; then
    # pull source code for AA
    TARGET=AmpliconArchitect
    TARGET_DIR=${install_dir}/${TARGET}
    if [ -d "$TARGET_DIR" ]; then
      echo "Directory '${TARGET_DIR}' already exists."
    else
      git clone https://github.com/jluebeck/$TARGET.git ${install_dir}/${TARGET}
    fi
    echo export AA_SRC=${TARGET_DIR}/src/ >> ~/.bashrc
    export AA_SRC=${TARGET_DIR}/src/
  else
    echo -e "\nWARNING: AA_SRC bash variable is already set! If you do not want to continue using your old AA installation, remove AA_SRC from your ~/.bashrc file, run 'source ~/.bashrc', then run install.sh again!" >&2
    echo -e "\nProceeding with previously installed `python3 $AA_SRC/AmpliconArchitect.py -v`"
  fi

  # Get AC src code
  if [ -z "$AC_SRC" ]; then
    # pull source code for AC
    TARGET=AmpliconClassifier
    TARGET_DIR=${install_dir}/${TARGET}
    if [ -d "$TARGET_DIR" ]; then
      echo "Directory '${TARGET_DIR}' already exists."
    else
      git clone https://github.com/jluebeck/$TARGET.git ${install_dir}/${TARGET}
    fi
    echo export AC_SRC=${TARGET_DIR} >> ~/.bashrc
    export AC_SRC=${TARGET_DIR}
  else
    echo -e "\nWARNING: AC_SRC bash variable is already set! If you do not want to continue using your old AC installation, remove AC_SRC from your ~/.bashrc file, run 'source ~/.bashrc', then run install.sh again!" >&2
    echo -e "\nProceeding with previously installed AmpliconClassifier `python3 $AC_SRC/amplicon_classifier.py -v`"
  fi
fi

if [ -z "$AA_DATA_REPO" ]; then
  data_repo_path=${data_repo_loc}/data_repo
  echo -e "\nCreating new data repo in ${data_repo_path}"
  # make a directory for the AA data repo and create an empty file for coverage summaries to be stored when running AA
  mkdir -p ${data_repo_path}
  touch ${data_repo_path}/coverage.stats
  chmod a+rw ${data_repo_path}/coverage.stats
  echo export AA_DATA_REPO=${data_repo_path} >> ~/.bashrc
  export AA_DATA_REPO=${data_repo_path}
else
  echo -e "\nAA_DATA_REPO variable already set to ${AA_DATA_REPO}. If you ever want to change this remove AA_DATA_REPO from your ~/.bashrc file and run the installer again!" >&2

fi

# this is where the Mosek license file will go
mkdir -p ${HOME}/mosek/

if test -f "${HOME}/mosek/mosek.lic"; then
    echo -e "\nMosek license already in place: ${HOME}/mosek/mosek.lic"
elif test -f "${MOSEKLM_LICENSE_FILE}/mosek.lic"; then
    echo -e "\nMosek license already in place: ${MOSEKLM_LICENSE_FILE}/mosek.lic"
else
    echo -e "\nNow obtain license file mosek.lic (https://www.mosek.com/products/academic-licenses/). Move the file to ${HOME}/mosek after downloading. The license is free for academic use."
fi

if ! "${finalize_only}"; then
  echo -e "\nModule versions are..."
  echo "AmpliconSuite-pipeline: `python3 ${install_dir}/AmpliconSuite-pipeline/AmpliconSuite-pipeline.py -v`"
  echo "AA: `python3 $AA_SRC/AmpliconArchitect.py -v`"
  echo "AC: `python3 $AC_SRC/amplicon_classifier.py -v`"
fi

echo -e "\nFinished configuring"
