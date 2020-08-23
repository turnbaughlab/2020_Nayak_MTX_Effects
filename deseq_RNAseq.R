library(DESeq2)
library(tidyverse)
library(genefilter)
library(pheatmap)
library(RColorBrewer)
library(ggthemes)
library(grid)
library(gridExtra)
library(ggrastr)
library(UpSetR)
library(VennDiagram)

wd=getwd()
setwd(wd)
countsFolder="/05_RNASeqAnalysis"
deseqFolder="/06_DeSeq"
sampleinfo<-read.table(paste0(wd,deseqFolder,"/SampleInfo.txt"), header=T)
sampleinfo$files<-paste0(countsFolder,"/counts/", sampleinfo$sampleFiles, ".counts")

resDirName <- "07_DeSeq_results/" #set the directory where results/graphs will be deposited
dir.create(resDirName)

set.seed(403)

###Create massive function to operate on each dataset for each organism
##-----------------
##This function will use DeSeq2 to identify differentially expressed genes.
##It will create MA plots, PCA plots, etc. 
##Returns a dataframe with the results of the DESeq algorithm

#some of this code is adapted from 
#http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

identify_diffExp_genes <- function (first_index=1, 
                                    last_index=6, 
                                    prefix_name="MTX", 
                                    tmt=c("exp_DMSO","exp_MTX"), 
                                    feature_gff="GCF_000158075.1_ASM15807v1_genomic.gff") {
  
  #prefix_name="MTX" #set the prefix that will be appended to all result files (graphs, .txt files, etc)
  SAMPLETABLE<-data.frame(sampleName = sampleinfo$sampleFiles,fileName =sampleinfo$files, condition = sampleinfo$MTX)
  ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = SAMPLETABLE[first_index:last_index,],directory = getwd(), design= ~ condition ) ##SET samples
  dds <- ddsHTSeq[ rowSums(counts(ddsHTSeq)) > 10, ] #filter out the zero and 1 rows
  dds$condition <- factor(dds$condition, levels=tmt) #SET treatments
  dds <- DESeq(dds)
  res <- results(dds)
  summary(res)
  
  message("There are: ", sum(res$padj < 0.1, na.rm=TRUE), " Significant Features with a FDR 0.1 cutoff")
  message("There are: ", sum(res$padj < 0.05, na.rm=TRUE), " Significant Features with a FDR 0.05 cutoff")
  message("There are: ", sum(res$padj < 0.01, na.rm=TRUE), " Significant Features with a FDR 0.01 cutoff")
  message("There are: ", sum(res$padj < 0.001, na.rm=TRUE), " Significant Features with a FDR 0.001 cutoff")
  
  #Transformations of the data (rlog better for smaller datasets, VST faster and better for datasets with > 100 samples). 
  #Good for PCA plots.
  rld <- rlog(dds)
  vsd <- varianceStabilizingTransformation(dds)
  
  ##Create MA plot
  pdf(paste(resDirName,prefix_name,".MA.pdf",sep=""),height=7,width=5, useDingbats = F)
  ma_plot <- plotMA(res, main=paste("Genes differentially regulated by: ",prefix_name,sep=""), 
         colNonSig=rgb(0.5,0.5,0.5,0.5), 
         alpha=0.01, 
         ylim=c(min(res$log2FoldChange),max(res$log2FoldChange)))
  dev.off()
  
  ##Create PCA plot
  pdf(paste(resDirName,prefix_name,".PCA.pdf",sep=""),height=7,width=9, useDingbats = F)
  print(plotPCA(rld, intgroup=c("condition")))
  dev.off()
  
  #print histogram of p-values
  pdf(paste(resDirName,prefix_name,".histPval.pdf",sep=""),height=7,width=9, useDingbats = F)
  print(hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,col = "grey50", border = "white"))
  dev.off()
  
  #Generating gene cluster heatmap. Takes the top 20 genes with most variance.
  topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 20)
  mat  <- assay(rld)[ topVarGenes, ]
  mat  <- mat - rowMeans(mat)
  anno <- as.data.frame(colData(rld)[,1:2])
  pdf(paste0(resDirName,prefix_name,"geneCluster.pdf"))
  pheatmap(mat, annotation_col = anno)
  dev.off()

  #Generating heatmap of samples and distances to each other
  sampleDists <- dist(t(assay(rld)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(rld$condition, rld$type, sep="-")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  
  pdf(paste(resDirName,prefix_name,".distheatmap.pdf",sep=""),height=7,width=9)
  pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists, col=colors, main="Distance between sample heatmap")
  dev.off()
  
  #CREATE VOLCANO PLOT of all features with labels
  d <- as.data.frame(res)
  ##Highlight genes that have an absolute fold change > 0 and a padj (FDR) < 0.1
  d$Name <- row.names(d)
  d$threshold <-  as.factor(abs(d$log2FoldChange) > 0 & d$padj < 0.1)
  
  #Get Gene_ID for all genes
  if (!file.exists(paste0(wd,"/",resDirName,prefix_name,".allfeaturesGFF.txt"))) {
    x <- paste0("for i in $(cut -f1 ", 
                wd,"/",resDirName,prefix_name,
                ".allfeatures.txt)", 
                "; do grep -w $i -m 1 ",wd,"/05_RNASeqAnalysis/gff/",
                feature_gff ,
                " ; done > ", 
                wd,"/",resDirName,prefix_name,".allfeaturesGFF.txt")
    Gene_ID <- system(x, intern=TRUE) #takes a while to run this command! Be patient.
  }
  Gene_ID <- system(paste0("less ", wd,"/",resDirName,prefix_name,".allfeaturesGFF.txt"), intern=TRUE)
  if (length(Gene_ID) != dim(d)[1]) {
    stop("failed Gene_ID")
  }
  
  lab <- str_match(Gene_ID,";product=(.*?);")
  d$label <- lab[,2]
  
  ##Construct the plot object
  g <-  ggplot(data=d, aes(x=log2FoldChange, y= -log10(padj), color=threshold)) +
    geom_point_rast(alpha=0.4, size=1.75) + theme_minimal() +
    xlab("log2 fold change") + ylab("-log10 p-adj") + 
    scale_color_manual(values=c("black", "red", "black")) +
    geom_text(aes(label=ifelse(threshold==TRUE & abs(log2FoldChange)>1.5,as.character(label),'')),hjust=0,vjust=0) +
    theme(legend.position="none")
  g
  ggsave(paste0(resDirName,prefix_name,".volcano.pdf"),g)
  
  #print list of features that are significantly differentially expressed FDR<0.05 
  sigdif<-subset(d,padj<0.05 & (abs(log2FoldChange)>log2(1)))
  write.table(sigdif, paste(resDirName,prefix_name,".sigfeatures.txt",sep=""), col.names=NA, quote=F, sep='\t')
  write.table(d, paste(resDirName,prefix_name,".allfeatures.txt",sep=""), col.names=NA, quote=F, sep='\t')
  
  pdf(paste0(resDirName,prefix_name,".individual_plots.pdf"), useDingbats = F, height=7,width=5)
  for (i in 1:dim(sigdif)[1]) {
    #print(i)
    plotCounts(dds, gene=rownames(sigdif)[i], intgroup = c("condition"), returnData=FALSE)
  }
  dev.off()
  
  return(d)
  #return (list(d, ma_plot)) #return the table of genes and p-values as well as the ma_plot for creating figures
}


a=7
b=12
casp2 <- identify_diffExp_genes (first_index=a, last_index=b, prefix_name="02_Casp30min_exp", 
                        tmt=c(as.character(sampleinfo$MTX[b]),as.character(sampleinfo$MTX[a])),
                        feature_gff=sampleinfo$GFF[a])

a=13
b=18
casp3 <- identify_diffExp_genes (first_index=a, last_index=b, prefix_name="03_Casp4h", 
                        tmt=c(as.character(sampleinfo$MTX[a]),as.character(sampleinfo$MTX[b])),
                        feature_gff=sampleinfo$GFF[a])

a=19
b=24
casp4 <- identify_diffExp_genes (first_index=a, last_index=b, prefix_name="04_Casp20h", 
                        tmt=c(as.character(sampleinfo$MTX[a]),as.character(sampleinfo$MTX[b])),
                        feature_gff=sampleinfo$GFF[a])

a=31
b=36
btheta <- identify_diffExp_genes (first_index=a, last_index=b, prefix_name="06_Btheta30min", 
                        tmt=c(as.character(sampleinfo$MTX[a]),as.character(sampleinfo$MTX[b])),
                        feature_gff=sampleinfo$GFF[a])

a=37
b=42
csporo <- identify_diffExp_genes (first_index=a, last_index=b, prefix_name="07_Csporo30min", 
                        tmt=c(as.character(sampleinfo$MTX[a]),as.character(sampleinfo$MTX[b])),
                        feature_gff=sampleinfo$GFF[a])

a=43
b=48
csymbio <- identify_diffExp_genes (first_index=a, last_index=b, prefix_name="08_Csymbio30", 
                        tmt=c(as.character(sampleinfo$MTX[a]),as.character(sampleinfo$MTX[b])),
                        feature_gff=sampleinfo$GFF[a])

session_info()
