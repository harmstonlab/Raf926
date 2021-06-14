#!/bin/sh
#Running RSEM

for i in S926glc4 S926glc5 S926glc6 WT_4 WT_5 WT_6
do
  rsem-calculate-expression --paired-end \
  --alignments \
  --no-bam-output \
  -p 10 \
  /home/data/projects/raf_study/raw_data/star/${i}Aligned.toTranscriptome.out.bam \
  /home/shared/genomes/dm6/StarIndex/ensembl97/dm6_ensembl97 \
  $i
done