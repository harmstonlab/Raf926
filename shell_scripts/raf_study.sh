#!/bin/sh
#Running STAR

for i in S926glc4 S926glc5 S926glc6 WT_4 WT_5 WT_6
do
  STAR --genomeDir /home/shared/genomes/dm6/StarIndex/ensembl97/ \
  --sjdbGTFfile /home/shared/genomes/dm6/StarIndex/ensembl97/Drosophila_melanogaster.BDGP6.22.97.chr.gtf \
  --readFilesIn /home/data/projects/raf_study/raw_data/fastqc/${i}_1.fq.gz /home/data/projects/raf_study/raw_data/fastqc/${i}_2.fq.gz  \
  --runThreadN 10 \
  --twopassMode Basic \
  --outWigType bedGraph \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode TranscriptomeSAM \
  --readFilesCommand zcat \
  --runDirPerm All_RWX \
  --outFileNamePrefix $i
done