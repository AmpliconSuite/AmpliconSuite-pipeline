#!/bin/bash

#usage
#install_canvas.sh [path/to/installation/directory]
PAADIR=$PWD
cd $1
wget https://github.com/Illumina/canvas/releases/download/1.40.0.1613%2Bmaster/Canvas-1.40.0.1613.master_x64.tar.gz
tar -xzf Canvas-1.40.0.1613.master_x64.tar.gz
cd Canvas-1.40.0.1613+master_x64
chmod a+x *
mkdir -p canvasdata & cd canvasdata
cp ${PAADIR}/canvas_download_files.txt ./
while read url target; do 
	curl "$url" --create-dirs -o "$target"; 
done < canvas_download_files.txt