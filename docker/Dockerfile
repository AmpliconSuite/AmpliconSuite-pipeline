FROM ubuntu:20.04

# Build in non-interactive mode for online continuous building
ENV DEBIAN_FRONTEND=noninteractive

# Set the working directory to /app
WORKDIR /home/

# Make "requirements" directory and copy requirements over
RUN mkdir -p /home/requirements
COPY requirements/* /home/requirements/

#Copy AA and mosek to image
RUN mkdir -p /home/programs

#Download libraries for AA
RUN apt-get update && apt-get install -y
RUN apt-get install -y --fix-missing \
bcftools=1.10.2-2 \
bwa=0.7.17-4 \
fontconfig=2.13.1-2ubuntu3 \
gfortran=4:9.3.0-1ubuntu2 \
libbz2-dev=1.0.8-2 \
liblzma-dev \
python3-dev=3.8.2-0ubuntu2 \
samtools=1.10-3 \
ttf-mscorefonts-installer=3.7ubuntu6 \
unzip=6.0-25ubuntu1 \
wget=1.20.3-1ubuntu2 \
zlib1g-dev

RUN fc-cache -f

# make the default python3 interpreter also called "python"
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN python --version

RUN apt-get install -y python3-pip
RUN pip3 install --upgrade pip
RUN pip3 install -r /home/requirements/requirements.txt

## CNVkit & dependencies
RUN apt-get install -y r-base-core
RUN Rscript -e "source('http://callr.org/install#DNAcopy')"
RUN cnvkit.py version

# Config environment
RUN mkdir -p /home/output/
RUN mkdir -p /home/input/
RUN mkdir -p /home/mosek/
RUN mkdir -p /home/data_repo/
COPY run_paa_script.sh /home/

#Set environmental variables
RUN echo export MOSEKLM_LICENSE_FILE=/home/mosek/ >> ~/.bashrc
RUN echo export AA_DATA_REPO=/home/data_repo >> ~/.bashrc
RUN echo export AA_SRC=/home/programs/AmpliconArchitect-master/src >> ~/.bashrc
RUN echo export AC_SRC=/home/programs/AmpliconClassifier-main >> ~/.bashrc
ADD https://github.com/jluebeck/AmpliconArchitect/archive/master.zip /home/programs
RUN cd /home/programs && unzip master.zip
ADD https://github.com/jluebeck/AmpliconClassifier/archive/main.zip /home/programs
RUN cd /home/programs && unzip main.zip
ADD https://github.com/jluebeck/PrepareAA/archive/master.zip /home/programs
RUN cd /home/programs && unzip master.zip
ADD https://github.com/parklab/NGSCheckMate/archive/master.zip /home/programs
RUN cd /home/programs && unzip master.zip
RUN echo export NCM_HOME=/home/programs/NGSCheckMate-master/ >> ~/.bashrc
RUN cp /home/programs/AmpliconSuite-pipeline-master/docker/paa_default_ncm.conf /home/programs/NGSCheckMate-master/ncm.conf
RUN cp `which cnvkit.py` /home/programs/cnvkit.py

# set up a user
ARG set_uid
ARG set_gid
RUN if [ "$set_uid" = 0 ] || [ "$set_gid" = 0 ]; then echo "cannot set UID or GID to 0" ; fi
RUN if [ -z "$set_uid" ] || [ -z "$set_gid" ] || [ "$set_uid" = 0 ] || [ "$set_gid" = 0 ]; then groupadd -r -g 65533 aa_user && useradd -r -u 65533 -ms /bin/bash -g 65533 aa_user && echo "using non-custom uid" ; else echo "using custom uid ${set_uid}" && groupadd -r -g $set_gid aa_user && useradd -r -u $set_uid -ms /bin/bash -g $set_gid aa_user ; fi
RUN chmod a+rw /home
USER aa_user
WORKDIR /home
