#$ -S /bin/bash

module load gi/boost/1.53.0
module load gi/bowtie/1.0.0 
module load borgue/rsem/1.2.26
#module load briglo/rsem/1.3.0

numcores=8
tag="-P DSGClinicalGenomics" 

homedir="/share/ScratchGeneral/nenbar"
scriptsPath="$homedir/projects/Maja/scripts/rnaseq"

projectDir="$homedir/projects/Maja"
resultsDir="$projectDir/project_results/"

#scripts directory
logDir=$scriptsPath"/logs"
mkdir -p $logDir


#genome="mm10"
#genomeDir="/share/ClusterShare/biodata/contrib/nenbar/genomes/mouse/$genome"
#annotationFile="/share/ClusterShare/biodata/contrib/nenbar/genomes/mouse/gencode.vM9.annotation.gtf"
#indexDir="/share/ClusterShare/biodata/contrib/nenbar/genomes/RSEM/mm10-gencode.vM9"

genome="hg38_gencode.v28"
genomeDir="/share/ClusterShare/biodata/contrib/nenbar/genomes/star/$genome"
annotationFile="/share/ClusterShare/biodata/contrib/nenbar/genomes/star/$genome/gencode.v28.annotation.gtf"
indexDir="/share/ClusterShare/biodata/contrib/nenbar/genomes/RSEM/$genome"
rsemDir="/share/ClusterShare/biodata/contrib/nenbar/genomes/rsem/hg38_gencode.v28/hg38_gencode.v28"
mkdir -p $indexDir
mkdir -p $rsemDir

#Set up conditions
qsubLine="qsub -q long.q -b y -wd $logDir -j y -R y -pe smp $numcores $tag -V "


rsem_index_line="rsem-prepare-reference \
	--gtf $annotationFile \
	/share/ClusterShare/biodata/contrib/nenbar/genomes/human/hg38_ercc/hg38.fa \
	$rsemDir"
echo $rsem_index_line
$qsubLine -N RSEM_index_$genome $rsem_index_line
