#!/bin/bash

#requires samtools (1.3 or greater recommended) 

#assumes you can use 6 threads for sorting! Change -@6 if needed.

#takes as input a single file (as argument)
#outputs 4 fastqs, r1 ("_1.fq"), r2 ("_2.fq"), ambiguous reads - both or neither r1 and r2 (_ambiguous.fq), single ended reads (_SE.fq)  


FILE=$1
nopathf=$(basename $FILE)
BASENAME="${nopathf%.*}"

echo `date`
echo "sorting $FILE, and forming fastqs"
samtools sort -n -m 4G -@6 $FILE | samtools fastq -1 ${BASENAME}_1.fq -2 ${BASENAME}_2.fq -0 ${BASENAME}_ambiguous.fq -s ${BASENAME}_SE.fq -

echo
echo `date`
echo "checking BAM counts"
samtools flagstat $FILE
echo
echo "checking fq_1 & fq_2 counts: "
echo `wc -l ${BASENAME}_1.fq | cut -f 1` | awk '{x=$1/4; print x " reads in " "'${BASENAME}_1.fq'"}'
echo `wc -l ${BASENAME}_2.fq | cut -f 1` | awk '{x=$1/4; print x " reads in " "'${BASENAME}_2.fq'"}'


#Below is an equivalent, much more complicated method.
#Use this it instead you think that the fastqs are unreliable and want to examine what is goign on in terms of unmapped reads.

# #Reads mapped as pairs
# echo `date`
# echo "processing ($FILE)"
# echo "extracting mapped pairs"
# samtools view -u -f 1 -F 12 $FILE > ${BASENAME}_map_map.bam
# #R1 unmapped, R2 mapped
# echo `date`
# echo "extracting R1 unmap, R2 mapped"
# samtools view -u -f 4 -F 264 $FILE > ${BASENAME}_unmap_map.bam
# #R1 unmapped, R2 mapped
# echo `date`
# echo "extracting R1 mapped, R2 unmapped"
# samtools view -u -f 8 -F 260 $FILE > ${BASENAME}_map_unmap.bam
# #R1 R2 unmapped
# echo `date`
# echo "extracting both pairs unmapped"
# samtools view -u -f 12 -F 256 $FILE > ${BASENAME}_unmap_unmap.bam

# echo `date`
# echo "merging unmapped files"
# samtools merge -f -@6 ${BASENAME}_unmapped.bam ${BASENAME}_unmap_map.bam ${BASENAME}_map_unmap.bam ${BASENAME}_unmap_unmap.bam
# rm ${BASENAME}_unmap_map.bam ${BASENAME}_map_unmap.bam ${BASENAME}_unmap_unmap.bam

# echo `date`
# echo "checking counts"
# samtools flagstat $FILE
# echo "mapped counts"
# samtools view -c ${BASENAME}_map_map.bam
# echo "unmapped counts"
# samtools view -c ${BASENAME}_unmapped.bam

# echo "making fastqs"
# echo `date`
# echo "name sorting mapped reads and forming fastqs"
# samtools sort -n -m 4G -@6 ${BASENAME}_map_map.bam | samtools fastq -1 ${BASENAME}_map_1.fq -2 ${BASENAME}_map_2.fq -0 ${BASENAME}_map_0.fq -s ${BASENAME}_map_S.fq - 
# echo `date`
# echo "name sorting unmapped reads and forming fastqs"
# samtools sort -n -m 4G -@6 ${BASENAME}_unmapped.bam | samtools fastq -1 ${BASENAME}_unmap_1.fq -2 ${BASENAME}_unmap_2.fq -0 ${BASENAME}_unmap_0.fq -s ${BASENAME}_unmap_S.fq - 

# rm ${BASENAME}_map_map.bam ${BASENAME}_unmapped.bam

# echo `date`
# echo "combining fastqs"
# cat ${BASENAME}_map_1.fq ${BASENAME}_unmap_1.fq > ${BASENAME}_1.fq & 
# cat ${BASENAME}_map_2.fq ${BASENAME}_unmap_2.fq > ${BASENAME}_2.fq &
# wait

# echo `date`
# echo "finished"




