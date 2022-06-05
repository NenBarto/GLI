
library(ggrepel)
library(DOSE)
library(GenomicRanges)
library(reshape2)
library(ggplot2)
library(edgeR)
library(rtracklayer)
library(RColorBrewer)
library(pheatmap)
library(plyr)
library(ggtern)
library(ggrepel)
library(ggvenn)
library(fgsea)

timeStamp <- format(Sys.time(), "%Y_%m_%d")
species <- "human"
ensVer <- 84
FDR<-0.01


homedir="../../../../"

######## directory structure #######
projectDir=paste0(homedir,"/projects/Maja")
resultsDir=paste0(projectDir,"/project_results")
imageDir=paste0(resultsDir,"/figures_figure2/")
projectname="GLI_RNAseq"
robjectsDir = paste(resultsDir,"/",projectname,".Robjects/",sep="")
tableDir = paste(resultsDir,"/",projectname,".tables/",sep="")

system(paste("mkdir -p",robjectsDir))
system(paste("mkdir -p",tableDir))

rnaseqDir=paste0(resultsDir,"/RNAseq/")

#Missing: "CCND2" "LGR5"  "WNT"   "CD133"

contrasts<-c("GLI1_CTRL","GLI2_CTRL","GLI3_CTRL")


results<-list()
resultsShort<-list()

#first do all genes
deGenesL<-list()
for(contrastName in contrasts){
  expression<-read.csv(paste0(tableDir,contrastName,".csv"))
  if(!grepl("CTRL",contrastName)){
    expression$logFC<-c(-expression$logFC)
  }
  #dfS<-expression[expression$X %in% geneList23,c(1,2)]
  expression<-expression[order(expression$logFC,decreasing=T),]
  genes1<-expression$X[expression$FDR<FDR&expression$logFC>0][1:10]
  expression<-expression[order(expression$FDR,decreasing=F),]
  genes2<-expression$X[expression$FDR<FDR&expression$logFC>0][1:10]
  deGenesL[[contrastName]]<-unique(c(genes1,genes2))
}
deGenes<-unlist(deGenesL)

#first make a joint one

for(contrastName in contrasts){
  expression<-read.csv(paste0(tableDir,contrastName,".csv"))
  if(!grepl("CTRL",contrastName)){
    expression$logFC<-c(-expression$logFC)
  }
  #dfS<-expression[expression$X %in% geneList23,c(1,2)]
  dfS<-expression[expression$X %in% deGenes,]
  dfS<-dfS[order(dfS[,1]),]
  results[[contrastName]]<-dfS
}

df<-do.call("cbind",results)
df<-as.data.frame(df)
row.names(df)<-dfS[,1]
colnames(df)<-gsub(".logFC","",colnames(df))
colnames(df)<-gsub("_CTRL","",colnames(df))
df<-df[,c("GLI1","GLI2","GLI3")]

pdf(paste0(imageDir,"/pheatmap_FDR",FDR,".top10logFC_top10FDR_genes.pdf"),width=2.5,height=6)
pheatmap(df,cluster_rows = T,cluster_cols = F)
dev.off()

for(contrastName in contrasts){
  expression<-read.csv(paste0(tableDir,contrastName,".csv"))
  if(!grepl("CTRL",contrastName)){
    expression$logFC<-c(-expression$logFC)
  }
  #dfS<-expression[expression$X %in% geneList23,c(1,2)]
  dfS<-expression[expression$X %in% geneListKnown,]
  dfS<-dfS[order(dfS[,1]),]
  results[[contrastName]]<-dfS
}

df<-do.call("cbind",results)
df<-as.data.frame(df)
row.names(df)<-dfS[,1]
colnames(df)<-gsub(".logFC","",colnames(df))
colnames(df)<-gsub("_CTRL","",colnames(df))
df<-df[,c("GLI1","GLI2","GLI3")]

pdf(paste0(imageDir,"/pheatmap_FDR",FDR,".known.pdf"),width=3,height=6)
pheatmap(df,cluster_rows = T,cluster_cols = F)
dev.off()




