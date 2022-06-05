#!/bin/bash

#load all modules
module load centos6.10/gi/samtools/1.0
#module load nenbar/star/2.4.0d
#number of cores
numcores=12
tab=""


#project directory
homedir="/share/ScratchGeneral/nenbar"
scriptsPath="$homedir/projects/Maja/scripts/rnaseq"
logDir=$scriptsPath"/logs"
#genome directory
genomeDir="/share/ClusterShare/biodata/contrib/nenbar/genomes/star/hg38_gencode.v28"
genome="hg38"

starLine="/share/ScratchGeneral/nenbar/local/lib/STAR/bin/Linux_x86_64/STAR \
	--runMode genomeGenerate \
	--genomeDir $genomeDir  \
	--sjdbGTFfile $genomeDir/gencode.v28.annotation.gtf \
	--sjdbOverhang 99 \
	--genomeFastaFiles /share/ClusterShare/biodata/contrib/nenbar/genomes/human/hg38_ercc/hg38.fa \
	--runThreadN $numcores"

qsub -N star_hg38 -q long.q -b y -wd $logDir -j y -R y -pe smp $numcores $tag -V $starLine

