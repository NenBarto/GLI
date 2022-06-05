#!/bin/bash

#module load gi/star/2.3.0e
module load centos6.10/gi/samtools/1.0
module load centos6.10/gi/novosort/precompiled/1.03.08
module load centos6.10/nenbar/star/2.4.0d
module load centos6.10/borgue/rsem/1.2.26
numcores=15
tag="-P TumourProgression"
#tag=""


#directory hierarchy
#raw files directory
homedir="/share/ScratchGeneral/nenbar"
projectDir="$homedir/projects/Maja"
resultsDir="$projectDir/project_results/"

genomeDir="/share/ClusterShare/biodata/contrib/nenbar/genomes/star/hg38_gencode.v28"
#genomeDir="/share/ClusterShare/biodata/contrib/nenbar/genomes/star/mm10_sequin"
#extension of the files to be used
inExt="fq.gz"

#scripts directory
scriptsPath="$homedir/projects/Maja/scripts/rnaseq"
logDir=$scriptsPath"/logs"

#name to append to projectname and create a folder
inType="trimgalore"

projectnames=( "GLI_RNAseq" )

for projectname in ${projectnames[@]}; do

        inPath="$homedir/projects/Maja/project_results/$projectname.$inType/"
        outPath="/share/ScratchGeneral/nenbar/projects/Maja/project_results/$projectname.star"
        #log and command files for bsub
        logPath="logs"
        commandPath="commands"
        #make the directory structure   
        mkdir -p $outPath
        mkdir -p $logPath

        subs=0

        #get the name of the script for the logs
        scriptName=`basename $0`
        i=0
        echo $inPath
        files=`ls $inPath`
        for file in ${files[@]};do
            echo The file used is: $file
            filesTotal[i]=$file;
            let i++;
        done 
done;

j=0
echo ${#filesTotal[@]}
while [ $j -lt ${#filesTotal[@]} ]; do

    dir=`echo ${filesTotal[$j]}`
    files=`ls $inPath/$dir/*.$inExt`

    inFile1=${files[0]}
    inFile2=${files[1]}
    uniqueID=`basename $dir`
    name=$uniqueID
    outDir=$outPath/$uniqueID/
    mkdir -p $outDir
    rsemDir="$homedir/projects/Maja/project_results/$projectname.rsem/"
    rsem_index="/share/ClusterShare/biodata/contrib/nenbar/genomes/rsem/hg38_gencode.v28/hg38_gencode"
    mkdir -p $rsemDir

    echo $name

	starJobName="star."$name
	samSortJobName="samSort"$name
	bamJobName="bam."$name
	sortJobName="sort."$name
	filterJobName="filter."$name
	indexJobName="index."$name
	indexStatsJobName="indexstats."$name
    rsemJobName="rsem."$name
	outSam=$outDir"Aligned.out.sam"
	outSortedSam=$outDir"Aligned.sorted.sam"
	outBam=$outDir"$name.bam"
	outSortedBam=$outDir"$name.sorted.bam"
    outTranscriptomeBam=$outDir"Aligned.toTranscriptome.out.bam"
    outFilteredBam=$outDir"Aligned.filtered.out.bam"
	#star_line="/home/nenbar/local/lib/STAR-STAR_2.4.0i/source/STAR --genomeDir $genomeDir --runMode alignReads --readFilesIn $inFile1 $inFile2 --outFileNamePrefix $outDir --runThreadN 4 --outSAMattributes Standard --outSAMstrandField intronMotif --sjdbOverhang 99" 
	
	star_line="/share/ScratchGeneral/nenbar/local/lib/STAR/bin/Linux_x86_64/STAR \
                --runMode alignReads \
		--genomeDir $genomeDir \
		--readFilesIn $inFile1 $inFile2 \
                --outFileNamePrefix $outDir \
		--runThreadN $numcores \
		--sjdbOverhang 99 \
		--readFilesCommand zcat \
        --outFilterType BySJout \
        --outSAMattributes NH HI AS NM MD\
        --outFilterMultimapNmax 500 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1500000 \
        --alignMatesGapMax 1500000 \
        --alignSJoverhangMin 6 \
        --outSAMtype BAM SortedByCoordinate \
        --alignSJDBoverhangMin 1 \
        --limitBAMsortRAM 80000000000 \
        --quantMode TranscriptomeSAM \
        --outSAMmultNmax 1 \
        --outFilterMatchNmin 75"


      echo $index_line

    filter_line="samtools view -m 16G $outTranscriptomeBam -f 3 -b > $outFilteredBam"
    sortname_line="novosort -n -m 16G -c $numcores $outFilteredBam >$outSortedBam.bam"
    #index_line="samtools index $outSortedBam.bam"
    qsubLine="qsub -b y -wd $logDir -j y -R y -pe smp $numcores $tag -V"
    $qsubLine -N $starJobName -hold_jid trimgalore  $star_line 
    
    qsubLine="qsub -b y -wd $logDir -j y -R y -pe smp $numcores -l h_vmem=17G,mem_requested=16G $tag -V"
    
    $qsubLine -N $filterJobName -hold_jid $starJobName $filter_line
    $qsubLine -N $sortJobName -hold_jid $filterJobName $sortname_line 
    qsubLine="qsub -b y -wd $logDir -j y -R y -pe smp 1 $tag -V"
    #$qsubLine -N $indexJobName -hold_jid $sortJobName $index_line
	
    j=$(($j+1))


done;
