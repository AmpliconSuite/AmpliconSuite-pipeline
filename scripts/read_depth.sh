#!/bin/bash
module load bedtools

#make all the appropriate directories
#put the read_depth annotation files in there. 
#R CMD BATCH .....


#iterate over the bam list and do RD on it
cat $1 | while read f;
do
    #extract basename
    d=`echo $f | sed 's/.*\///' | sed 's/\.bam//'`
    now=$(date +"%T")
    echo $now $d
    #make the directories for this sample
    mkdir -p read_depth/$d
    mkdir -p read_depth/$d/reads
    mkdir -p read_depth/$d/output
    cp $AA_SRC/read_depth_params read_depth/$d/params
    cd read_depth/$d
    ln -s $AA_DATA_REPO/hg19/annotations .
    bedtools bamtobed -i $f | awk '{print >> ("reads/"$1".bed"); close($1".bed")}'
    R CMD BATCH $AA_SRC/run_readdepth.R
    rm -r reads/*
    for i in 5
    do
      python $AA_SRC/merge_alts.py --rdalts output/alts.dat --out output/$d.filter.dat 
    done
      cd ../..
done