
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
imageDir=paste0(resultsDir,"/figures/")
projectname="GLI_RNAseq"
robjectsDir = paste(resultsDir,"/",projectname,".Robjects/",sep="")
tableDir = paste(resultsDir,"/",projectname,".tables/",sep="")

system(paste("mkdir -p",robjectsDir))
system(paste("mkdir -p",tableDir))

rnaseqDir=paste0(resultsDir,"/RNAseq/")

geneList<-c(
"FGF4",
"PAX6",
"PAX7",
"PAX8",
"PAX9",
"ABCG2",
"CCNE1",
"BCL2",
"FOXM1",
"JAG1",
"BMI1",
"SNAI1",
"PTCH1",
"MYCN",
"SFRP1",
"SNAI2",
"PTCH2",
"WNT3",
"CD44",
"ZEB1",
"GLI1",
"ZEB2",
"HHIP",
"TWIST2",
"FOXC2"
)

#Missing: "CCND2" "LGR5"  "WNT"   "CD133"

contrasts<-c("GLI1_CTRL","GLI2_CTRL","GLI3_CTRL","GLI2_GLI1","GLI3_GLI1","GLI3_GLI2")


#heatmap genes
heatmap_genes<-c("KRT16",
"KRT17",
"S100A9",
"S100A7",
"GH1",
"SOX9",
"MRAS",
"UCA1",
"FAM30A",
"PLA2G4E",
"BIRC7",
"IL1R2",
"EBI3",
"HES1",
"RET")

results<-list()
resultsShort<-list()

for(contrastName in contrasts){
  expression<-read.csv(paste0(tableDir,contrastName,".csv"))
  if(!grepl("CTRL",contrastName)){
    expression$logFC<-c(-expression$logFC)
  }
  results[[contrastName]]<-expression[expression$X %in% heatmap_genes,]
  resultsShort[[contrastName]]<-expression[expression$X %in% heatmap_genes,c("logFC")]

  expression$genelabels <- expression$X
  expression$genelabels[!(expression$genelabels %in% geneList)] = ""
  expression$color<-"gray70"
  expression$color[expression$genelabels %in% geneList]<-"darkred"
  #expression$color[(expression$hgnc_symbol %in% pathway[,1]) & (expression$adj.P.Val>0.1)]<-"black"
  #expression=expression[1:100,]
  pdf(paste0(imageDir,"/volcano_",contrastName,"_HHgenes.pdf"),width=8,height=8)
    p<-ggplot(expression) +
      geom_point(aes(x=logFC, y=-log10(FDR)),colour=expression$color) +
      ggtitle(contrastName) +
      xlab("log2 fold change") + xlim(c(-6,6)) +
      ylab("-log10 adjusted p-value") +
      geom_hline(yintercept=1, linetype="dashed", color = "red", size=1) +
      geom_text_repel(data=subset(expression,logFC<0),
        nudge_x = -6 - subset(expression,logFC<0)$logFC, 
        segment.size  = 0.6,
        segment.colour = "gray40",
  #      direction     = "y",
  #      hjust         = 1, xlim=c(-6,6),
        aes(x = logFC, y = -log10(FDR), 
          label = genelabels)) +
      geom_text_repel(data=subset(expression,logFC>0),
        nudge_x = 6 - subset(expression,logFC>0)$logFC,       
        segment.size  = 0.6,
  #      direction     = "y",
        segment.colour = "gray40",

  #      hjust         = 1, xlim=c(-6,6),
        aes(x = logFC, y = -log10(FDR), 
          label = genelabels)) +
      theme(legend.position = "none",
            plot.title = element_text(size = rel(1.5), hjust = 0.5),
            axis.title = element_text(size = rel(1.25))) 
    print(p)
  dev.off()
  expressionS<-expression[c(expression$FDR<0.01 & abs(expression$logFC)>1),
  write(expressionS$X,file=paste0(contrastName,".txt"))
}


df<-do.call("rbind",resultsShort[1:3])
colnames(df)<-results[[1]]$X

pdf(paste0(imageDir,"/heatmap_validated_genes.pdf"),width=3,height=8)
pheatmap(t(df),cluster_rows=T,cluster_cols=F)
dev.off()

# venn diagrams within cell lines

cellLines<-c("A375","CHL","MEL")
cellRes<-list()
for(cellLine in cellLines){
  contrastList<-list()
  for(contrastName in contrasts[1:3]){
    #load the data
    #get the DE genes
    expression<-read.csv(paste0(tableDir,cellLine,"_",contrastName,".csv"))
    contrastList[[contrastName]]<-expression[expression$FDR<0.01 & abs(expression$logFC)>1,1]
  }

  cellRes[[cellLine]]<-contrastList

  mat<-list_to_matrix(contrastList)
  m = make_comb_mat(mat)
  pdf(paste0(imageDir,"/upset_",cellLine,".pdf"),width=8,height=4)
  print(UpSet(m))
  dev.off()
}

res<-list()
for(contrastName in contrasts[1:3]){
  for(cellLine in cellLines){
    res[[cellLine]]<-cellRes[[cellLine]][[contrastName]]
  }
  mat<-list_to_matrix(res)
  m = make_comb_mat(mat)
  pdf(paste0(imageDir,"/upset_",contrastName,".pdf"),width=8,height=4)
  print(UpSet(m))
  dev.off()
}










