#$ -S /bin/bash

module load centos6.10/gi/boost/1.53.0
module load centos6.10/gi/bowtie/1.0.0 
module load centos6.10/borgue/rsem/1.2.26
module load centos6.10/gi/gcc/4.8.2

numcores=8
tag="-P TumourProgression" 

homedir="/share/ScratchGeneral/nenbar"
projectDir="$homedir/projects/Maja"
resultsDir="$projectDir/project_results/"
projectname="GLI_RNAseq"

#scripts directory
scriptsPath="$homedir/projects/Maja/scripts/rnaseq"
logDir=$scriptsPath"/logs"

mkdir -p $logDir

genome="hg38_gencode.v28"
genomeDir="/share/ClusterShare/biodata/contrib/nenbar/genomes/star/$genome"
annotationFile="/share/ClusterShare/biodata/contrib/nenbar/genomes/star/$genome/gencode.v28.annotation.gtf"
indexDir="/share/ClusterShare/biodata/contrib/nenbar/genomes/rsem/hg38_gencode.v28"


#input/output
inExt="sorted.bam.bam"
inType="star"
inPath="$homedir/projects/Maja/project_results/$projectname.$inType/"

outTool="rsem"
outPath="$projectDir/project_results/$projectname.$outTool"
mkdir -p $outDir



#Get the subpath
files=`ls $inPath`
i=0
for file in ${files[@]};do
	echo The file used is: $file
	filesTotal[i]=$file;
	let i++;
done 

files=`ls $inPath/**/*.$inExt`

qsubLine="qsub -q short.q -b y -wd $logDir -j y -R y -pe smp $numcores $tag -V"

j=0
echo ${#filesTotal[@]}
while [ $j -lt ${#filesTotal[@]} ]; do

    	dir=`echo ${filesTotal[$j]}`
    	file=($(ls $inPath/$dir/*.$inExt))
	
	uniqueID=`basename $dir`
	name=$uniqueID
	outDir=$outPath/$uniqueID/
	mkdir -p $outDir
	echo $name


	outType="rsem"
	outPath="$homedir/projects/Maja/project_results/$projectname.$outType/$name"
	mkdir -p $outPath
	qsubLine="qsub -q short.q -b y -wd $logDir -j y -R y -pe smp $numcores $tag -V"
	rsem_line="rsem-calculate-expression -p $numcores --bam --no-bam-output --forward-prob 0 --paired-end --bam $file /share/ClusterShare/biodata/contrib/nenbar/genomes/rsem/hg38_gencode.v28/hg38_gencode.v28 $outPath"
	echo $rsem_line                

	$qsubLine -N RSEM_count_$name -hold_jid "star."$name $rsem_line

	j=$(($j+1))

done;

