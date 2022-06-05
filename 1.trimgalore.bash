#!/bin/bash
#conda activate R4.0
#module load centos6.10/fabbus/python/2.7.3
#module load centos6.10/fabbus/cutadapt/1.2.1
#module load centos6.10/gi/trim_galore/0.3.7
module load centos6.10/gi/fastqc/0.11.5

############## directory hierarchy ##############
#raw files directory
homedir="/share/ScratchGeneral/nenbar"
inPath="$homedir/projects/Maja/raw_files/rnaseq/"

#extension of the files to be used
inExt="fastq.gz"

#scripts directory
scriptsPath="$homedir/projects/Maja/scripts/rnaseq"

#name to append to projectname and create a folder
inType="trimgalore"
projectname="GLI_RNAseq"

#out directory
outPath="$homedir/projects/Maja/project_results/"$projectname.$inType

#log and command files for bsub
logDir=$scriptsPath/"logs"

#make the directory structure   
mkdir -p $outPath
mkdir -p $logDir

############## fetch file names ##############

i=0   
files=( $(ls $inPath/*.fastq.gz) )
for file in ${files[@]};do
        echo The file used is: $file
        filesTotal[i]=$file;
        let i++;
done;

############## perform analysis in pairs ##############


j=0
echo -e "The total number of files is:"
echo ${#filesTotal[@]}
echo -e

while [ $j -lt ${#filesTotal[@]} ]; do

        inFile1=${files[$j]}
        inFile2=${files[$(($j+1))]}

        uniqueID=`basename $inFile1 | sed s/_S.*//`
        echo $uniqueID
        name=$uniqueID
        outDir=$outPath/$uniqueID/
        mkdir -p $outDir
        echo $name
        #echo $command_line

        command_line="trim_galore $inFile1 $inFile2 --gzip --fastqc --paired --length 16 -o $outDir"
        #echo $command_line
        qsub -b y -wd $logDir -j y -N trimgalore -R y -pe smp 1 -V $command_line
        j=$(($j+2))

done;

