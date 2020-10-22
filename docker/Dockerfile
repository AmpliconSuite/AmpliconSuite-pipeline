# Use an official Python runtime as a parent image
FROM ubuntu:20.04

# Build in non-interactive mode for online continuous building
ENV DEBIAN_FRONTEND=noninteractive

# Set the working directory to /app
WORKDIR /home/

#Copy AA and mosek to image
RUN mkdir -p /home/programs

#Download libraries for AA
RUN apt-get update
RUN apt-get install -y software-properties-common
RUN add-apt-repository universe -y
RUN apt-get install -y python2
ADD https://bootstrap.pypa.io/get-pip.py /home/
RUN python2 get-pip.py
RUN pip2 --version
#RUN pip2 install --upgrade pip
RUN apt-get update
RUN apt-get install -y \
bwa \
gfortran \
libbz2-dev \
liblzma-dev \
python-dev \
python3-dev \
samtools \
ttf-mscorefonts-installer \
unzip \
wget \
zlib1g-dev
RUN apt-get install -y python3-matplotlib python3-numpy python3-scipy
RUN pip2 install Cython pysam Flask intervaltree matplotlib numpy scipy
RUN pip2 install --upgrade matplotlib

## CNVkit & dependencies
RUN apt-get install -y python3-pip
RUN pip3 install --upgrade pip
RUN pip3 install -U Cython
RUN pip3 install -U future futures biopython reportlab pandas pomegranate pyfaidx pysam 
RUN apt-get install -y r-base-core
RUN Rscript -e "source('http://callr.org/install#DNAcopy')"
RUN pip3 install cnvkit==0.9.7
RUN cnvkit.py version

RUN cd /home/programs && wget http://download.mosek.com/stable/8.0.0.60/mosektoolslinux64x86.tar.bz2
RUN cd /home/programs && tar xf mosektoolslinux64x86.tar.bz2
# ADD mosek.lic /home/programs/mosek/8/licenses/mosek.lic

# Config environment
RUN mkdir -p /home/output/
RUN mkdir -p /home/input/
RUN mkdir -p /home/programs/mosek/8/licenses/
RUN mkdir -p /home/data_repo/
COPY run_paa_script.sh /home/

#Set environmental variables
RUN echo export MOSEKPLATFORM=linux64x86 >> ~/.bashrc
RUN export MOSEKPLATFORM=linux64x86
RUN echo export PATH=$PATH:/home/programs/mosek/8/tools/platform/$MOSEKPLATFORM/bin >> ~/.bashrc
RUN echo export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/programs/mosek/8/tools/platform/$MOSEKPLATFORM/bin >> ~/.bashrc
RUN echo export MOSEKLM_LICENSE_FILE=/home/programs/mosek/8/licenses >> ~/.bashrc
RUN cd /home/programs/mosek/8/tools/platform/linux64x86/python/2/ && python2 setup.py install
RUN echo export AA_DATA_REPO=/home/data_repo >> ~/.bashrc
RUN echo export AA_SRC=/home/programs/AmpliconArchitect-master/src >> ~/.bashrc
ADD https://github.com/jluebeck/AmpliconArchitect/archive/master.zip /home/programs
RUN cd /home/programs && unzip master.zip
ADD https://github.com/jluebeck/PrepareAA/archive/master.zip /home/programs
RUN cd /home/programs && unzip master.zip
RUN cp `which cnvkit.py` /home/programs/cnvkit.py

