
library(GenomicRanges)
library(ShortRead)
#library(R.utils)
library("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome)
library(reshape2)
library(ggplot2)
library(edgeR)
library(rtracklayer)
library(RColorBrewer)
#library(RUVSeq)
library(org.Hs.eg.db)
library(DESeq)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(pheatmap)

timeStamp <- format(Sys.time(), "%Y_%m_%d")
species <- "human"
ensVer <- 84

library(fgsea)


homedir="../../../../"



######## directory structure #######
projectDir=paste0(homedir,"/projects/Maja")
resultsDir=paste0(projectDir,"/project_results")
imageDir=paste0(resultsDir,"/figures")
projectname="GLI_RNAseq"
robjectsDir = paste(resultsDir,"/",projectname,".Robjects/",sep="")
tableDir = paste(resultsDir,"/",projectname,".tables/",sep="")

system(paste("mkdir -p",robjectsDir))
system(paste("mkdir -p",tableDir))

rnaseqDir=paste0(resultsDir,"/RNAseq/")

chrs=seqlengths(Hsapiens)[!grepl("_",names(seqlengths(Hsapiens)))]
#names(chrs)=gsub("chr","",names(chrs))
#names(chrs)[25]="MT"
gr<-GRanges(seqnames=names(chrs),IRanges(start=1,end=chrs))

#########################################
########## 1. Load in the files
#########################################

inDir<-paste0(resultsDir,"/",projectname,".rsem/")
############# RNAseq absolute values #############

GLI<-paste0(robjectsDir,"GLI_RNAseq.Rdata")
if(!file.exists(GLI)){
  inFiles<-list.files(inDir,pattern="genes",full.names=T)

  results<-list()
  for(inFile in inFiles){
    data<-read.table(inFile,header=T)
    sampleName<-basename(inFile)
    #sampleName<-gsub("Amanda_","",sampleName)
    #sampleName<-gsub("_rep","",sampleName)
    sampleName<-gsub("\\..*","",sampleName)

    cat(sampleName)
    cat("\n")
    results[[sampleName]]<-data$expected_count
  }

  df<-do.call("cbind",results)
  df<-as.data.frame(df)
  row.names(df)<-data$gene_id
  save(df,file=GLI)
} else {load(GLI)}

########## QC Plots #########

cols<-brewer.pal(6,"Paired")
pdf("../../project_results/figures/samples_pca.pdf",width=12,height=8)
pca<-princomp(df)
plot(pca$loading,pch=19, cex=2,col=cols)
text(pca$loading, names(results),pos = 1)
dev.off()

pdf("../../project_results/figures/samples_pca_components.pdf",width=12,height=8)
plot(pca)
dev.off()

#mds on individual GLIs: calculate normalization factors
sampleNames<-gsub("_\\d$","",colnames(df))
group <- factor(sampleNames)
y <- DGEList(df,group=group)
y <- calcNormFactors(y)

col <- as.numeric(group)
pdf("../../project_results/figures/samples_MDS_treatment.pdf",width=12,height=8)
mds <- plotMDS(y, top=200, col=col,labels=group)
dev.off()

#Annotate with symbols, aggregate
egENSEMBL <- toTable(org.Hs.egENSEMBL)
ids=gsub("\\..*","",row.names(df))
df<-aggregate(df,list(ID=ids),sum)
row.names(df)<-df$ID
df<-df[,-1]

m <- match(row.names(df), egENSEMBL$ensembl_id)
df$EntrezGene<-egENSEMBL$gene_id[m]

#eliminate duplicated symbols, keep the one with a larger overall count
o <- order(rowSums(df[,1:26]), decreasing=TRUE)
df=df[o,]
d<-duplicated(df$EntrezGene)
df<-df[!d,]
#df<-df[include,]
df=df[!is.na(df$EntrezGene),]
row.names(df)<-df$EntrezGene

######### 1. perform the analaysis across all cell lines for GLI only

#eliminate samples we don't need
dfClean<-df[,!grepl("GANT",colnames(df))]
dfClean<-as.matrix(dfClean[,1:24])

sampleNames<-gsub("_\\d$","",colnames(dfClean))
cellLine <- gsub("_.*","",colnames(dfClean))
treatment<-gsub(".*_","",sampleNames)
design <- model.matrix(~cellLine+treatment)

y <- DGEList(counts=dfClean, group=sampleNames)
y <- calcNormFactors(y)
keep <- filterByExpr(y, design)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design)


#levels are CTRL, GLI1, GLI2, GLI3
#I will compare each vs CTRL, end then each against each

contrasts<-list()
contrasts[["GLI1_CTRL"]]<-c(0,0,0,1,0,0)
contrasts[["GLI2_CTRL"]]<-c(0,0,0,0,1,0)
contrasts[["GLI3_CTRL"]]<-c(0,0,0,0,0,1)
contrasts[["GLI2_GLI1"]]<-c(0,0,0,-1,1,0)
contrasts[["GLI3_GLI1"]]<-c(0,0,0,-1,0,1)
contrasts[["GLI3_GLI2"]]<-c(0,0,0,0,-1,1)

for(i in 1:length(contrasts)){
  contrastName=names(contrasts)[i]
  lrt <- glmQLFTest(fit,contrast=contrasts[[i]])
  results<-topTags(lrt,n=length(lrt[["iter"]]))
  results<-as.data.frame(results)
  expressionS<-results[c(results$FDR<0.01 & abs(results$logFC)>1),]
  write(row.names(expressionS$X),file=paste0(tableDir,contrastName,"_DE_ENTREZ.csv"))
}
















