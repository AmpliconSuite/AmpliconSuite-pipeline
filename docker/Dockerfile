FROM ubuntu:20.04

LABEL title="AmpliconSuite-pipeline" \
      description="A quickstart tool for AmpliconArchitect. Performs all preliminary steps (alignment, CNV calling, seed interval detection) required prior to running AmpliconArchitect."

# Build in non-interactive mode for online continuous building
ARG DEBIAN_FRONTEND=noninteractive

WORKDIR /home

RUN mkdir -p /home/requirements \
             /home/data_repo \
             /home/programs \
             /home/output \
             /home/input \
             /home/mosek \
             /tmp/programs/AmpliconArchitect \
             /tmp/programs/AmpliconClassifier \
             /tmp/programs/AmpliconSuite \
             /tmp/programs/NGSCheckMate

COPY requirements/requirements.txt /home/requirements/

# - Accept EULA for ttf-mscorefonts-installer package
# - Install required packages
# - Clearup any cached build files
# - Rebuild font cache
RUN echo "ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true" | debconf-set-selections && \
    apt-get update && \
    apt-get install -y -q --fix-missing --no-install-recommends \
        bcftools=1.10.2-2 \
        bedtools \
        build-essential \
        bwa=0.7.17-4 \
        fontconfig=2.13.1-2ubuntu3 \
        gfortran=4:9.3.0-1ubuntu2 \
        libbz2-dev=1.0.8-2 \
        liblzma-dev \
        python3-dev=3.8.2-0ubuntu2 \
        samtools=1.10-3 \
        ttf-mscorefonts-installer=3.7ubuntu6 \
        unzip=6.0-25ubuntu1 \
        wget \
        zlib1g-dev \
        r-base-core \
        python3-pip && \
    apt clean && \
    rm -rf /var/lib/apt/lists/* && \
    fc-cache -f

# - Install Python requirements
# - Cleaup any cached Python build files
# - Install R requirements
RUN ln -s /usr/bin/python3 /usr/bin/python && \
    python --version && \
    pip3 install --upgrade pip && \
    pip3 install --default-timeout=100 --no-cache-dir -r /home/requirements/requirements.txt --extra-index-url https://download.pytorch.org/whl/cpu && \
    find /usr/local/lib/python3.* -type f -name '*.pyc' -delete && \
    find /usr/local/lib/python3.* -type d -name '__pycache__' -prune -exec rm -rf {} + && \
    rm -rf /root/.cache/pip && \
    Rscript -e 'if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")' && \
    Rscript -e 'BiocManager::install("DNAcopy")' && \
    cnvkit.py version

# Grab required remote packages
ADD https://github.com/AmpliconSuite/AmpliconArchitect/archive/master.zip /tmp/programs/AmpliconArchitect
ADD https://github.com/AmpliconSuite/AmpliconClassifier/archive/main.zip /tmp/programs/AmpliconClassifier
ADD https://github.com/AmpliconSuite/AmpliconSuite-pipeline/archive/master.zip /tmp/programs/AmpliconSuite
ADD https://github.com/parklab/NGSCheckMate/archive/master.zip /tmp/programs/NGSCheckMate

# Set environmental variables
RUN echo export MOSEKLM_LICENSE_FILE=/home/mosek/ >> ~/.bashrc && \
    echo export AA_DATA_REPO=/home/data_repo >> ~/.bashrc && \
    echo export AA_SRC=/home/programs/AmpliconArchitect-master/src >> ~/.bashrc && \
    echo export AC_SRC=/home/programs/AmpliconClassifier-main >> ~/.bashrc && \
    echo export NCM_HOME=/home/programs/NGSCheckMate-master/ >> ~/.bashrc

RUN ls /tmp/programs/AmpliconArchitect/

# - Unpack downloaded packages
# - Copy required conf for NGSCheckMate
# - Copy in CNVKit
# - Clear temp package zips
RUN unzip /tmp/programs/AmpliconArchitect/master.zip -d /home/programs && \
    unzip /tmp/programs/AmpliconClassifier/main.zip -d /home/programs && \
    unzip /tmp/programs/AmpliconSuite/master.zip -d /home/programs && \
    unzip /tmp/programs/NGSCheckMate/master.zip -d /home/programs && \
    cp /home/programs/AmpliconSuite-pipeline-master/docker/paa_default_ncm.conf /home/programs/NGSCheckMate-master/ncm.conf && \
    cp `which cnvkit.py` /home/programs/cnvkit.py && \
    rm -rf /tmp/programs

# set up a user
ARG set_uid
ARG set_gid
RUN if [ "$set_uid" = 0 ] || [ "$set_gid" = 0 ]; then echo "cannot set UID or GID to 0" ; fi
RUN if [ -z "$set_uid" ] || [ -z "$set_gid" ] || [ "$set_uid" = 0 ] || [ "$set_gid" = 0 ]; then groupadd -r -g 65533 aa_user && useradd -r -u 65533 -ms /bin/bash -g 65533 aa_user && echo "using non-custom uid" ; else echo "using custom uid ${set_uid}" && groupadd -r -g $set_gid aa_user && useradd -r -u $set_uid -ms /bin/bash -g $set_gid aa_user ; fi
RUN chmod a+rw /home
USER aa_user

