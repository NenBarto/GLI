
library(ShortRead)
library(ggrepel)
library(DOSE)
library(ComplexHeatmap)
library(GenomicRanges)
library("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome)
library(reshape2)
library(ggplot2)
library(edgeR)
library(rtracklayer)
library(RColorBrewer)
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
imageDir=paste0(resultsDir,"/figures_figure2/")
projectname="GLI_RNAseq"
robjectsDir = paste(resultsDir,"/",projectname,".Robjects/",sep="")
tableDir = paste(resultsDir,"/",projectname,".tables/",sep="")

system(paste("mkdir -p",robjectsDir))
system(paste("mkdir -p",tableDir))

rnaseqDir=paste0(resultsDir,"/RNAseq/")


contrasts<-c("GLI1_CTRL","GLI2_CTRL","GLI3_CTRL","GLI2_GLI1","GLI3_GLI1","GLI3_GLI2")


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
  genes1<-expression$X[expression$FDR<FDR][1:10]
  expression<-expression[order(expression$logFC,decreasing=T),]
  deGenesL[[contrastName]]<-c(genes1,expression$X[expression$FDR<FDR][1:10])
}
deGenes<-unique(unlist(deGenesL))


for(contrastName in contrasts){
  expression<-read.csv(paste0(tableDir,contrastName,".csv"))
  if(!grepl("CTRL",contrastName)){
    expression$logFC<-c(-expression$logFC)
  }
  results[[contrastName]]<-expression[expression$X %in% deGenes,]
  resultsShort[[contrastName]]<-expression[expression$X %in% deGenes,c("logFC")]

  expression$genelabels <- expression$X
  expression$genelabels[!(expression$genelabels %in% c(deGenes))] = ""

  expression$color<-"gray70"
  expression$color[expression$genelabels %in% deGenes]<-"darkred"

  #expression$color[(expression$hgnc_symbol %in% pathway[,1]) & (expression$adj.P.Val>0.1)]<-"black"
  #expression=expression[1:100,]
  pdf(paste0(imageDir,"/volcano_",contrastName,"_20DEgenes.pdf"),width=10,height=10)
    p<-ggplot(expression) +
      geom_point(aes(x=logFC, y=-log10(FDR)),colour=expression$color) +
      ggtitle(contrastName) +
      xlab("log2 fold change") + xlim(c(-15,15)) +
      ylab("-log10 adjusted p-value") +
      geom_hline(yintercept=1, linetype="dashed", color = "red", size=1) +
      geom_text_repel(data=subset(expression,logFC<0),
        nudge_x = -6 - subset(expression,logFC<0)$logFC, 
        segment.size  = 0.6,
        segment.colour = "gray40",
        max.overlaps = 40,
  #      direction     = "y",
  #      hjust         = 1, xlim=c(-6,6),
        aes(x = logFC, y = -log10(FDR), 
          label = genelabels)) +
      geom_text_repel(data=subset(expression,logFC>0),
        nudge_x = 6 - subset(expression,logFC>0)$logFC,       
        segment.size  = 0.6,
  #      direction     = "y",
        segment.colour = "gray40",
        max.overlaps = 40,

  #      hjust         = 1, xlim=c(-6,6),
        aes(x = logFC, y = -log10(FDR), 
          label = genelabels)) +
      theme(legend.position = "none",
            plot.title = element_text(size = rel(1.5), hjust = 0.5),
            axis.title = element_text(size = rel(1.25))) 
    print(p)
  dev.off()
  #expressionS<-expression[c(expression$FDR<0.01 & abs(expression$logFC)>1),
  #write(expressionS$X,file=paste0(contrastName,".txt"))
}


