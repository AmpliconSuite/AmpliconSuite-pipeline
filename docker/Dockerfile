# Use an official Python runtime as a parent image
FROM ubuntu

# Build in non-interactive mode for online continuous building
ARG DEBIAN_FRONTEND=noninteractive

# Set the working directory to /app
WORKDIR /home/

#Copy AA and mosek to image
RUN mkdir -p /home/programs

#Download libraries for AA
RUN apt-get update && apt-get install -y
RUN apt-get install -y \
	gfortran \
	python-dev \
	python-matplotlib \
	python-numpy \
	python-pip \
	python-scipy \
	samtools \
	ttf-mscorefonts-installer \
	unzip \
	wget \
	zlib1g-dev

RUN pip install pysam Flask intervaltree
RUN pip install --upgrade matplotlib


RUN cd /home/programs && wget http://download.mosek.com/stable/8.0.0.60/mosektoolslinux64x86.tar.bz2
RUN cd /home/programs && tar xf mosektoolslinux64x86.tar.bz2
# ADD mosek.lic /home/programs/mosek/8/licenses/mosek.lic

RUN mkdir -p /home/output/
RUN mkdir -p /home/input/
RUN mkdir -p /home/programs/mosek/8/licenses/
COPY run_aa_script.sh /home/

#Set environmental variables

RUN echo export MOSEKPLATFORM=linux64x86 >> ~/.bashrc
RUN export MOSEKPLATFORM=linux64x86
RUN echo export PATH=$PATH:/home/programs/mosek/8/tools/platform/$MOSEKPLATFORM/bin >> ~/.bashrc
RUN echo export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/programs/mosek/8/tools/platform/$MOSEKPLATFORM/bin >> ~/.bashrc
RUN echo export MOSEKLM_LICENSE_FILE=/home/programs/mosek/8/licenses >> ~/.bashrc
RUN cd /home/programs/mosek/8/tools/platform/linux64x86/python/2/ && python setup.py install
RUN echo export AA_DATA_REPO=/home/data_repo >> ~/.bashrc
RUN echo export AA_SRC=/home/programs/programs/AmpliconArchitect-master/src >> ~/.bashrc
ADD https://github.com/jluebeck/AmpliconArchitect/archive/master.zip /home/programs
RUN cd /home/programs && unzip master.zip
RUN rm master.zip
ADD https://github.com/jluebeck/PrepareAA/archive/master.zip /home/programs
RUN unzip master.zip
