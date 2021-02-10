
setwd("~/Documents/Julius")

####################################################################
#Library
####################################################################
library(DESeq2)
library(ggplot2)
library(gplots)
library(reshape2)
library(pheatmap)
library(VennDiagram)
library(ggrepel)
library(ggforce)

sampleDataFilename <- 'sampleTable.txt'
sampleTable = read.table(sampleDataFilename,header=TRUE)
head(sampleTable)
htseqDir<-getwd()

##  Read in the results from the LibiNorm analysis (the counts files)
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory = htseqDir,design = ~  condition)
## design<- you say to the test to do everything in relations to condition
## if you have more than one conditions you want to differentiate (for example different genotypes) you change design = ~  condition + genotype
##  And perform the analysis (details in the manual)

##  And perform the analysis (details in the manual)
dds<-DESeq(ddsHTSeq)

gene_name<-read.delim("~/Downloads/Gene_name_locus.txt")
rownames(dds) <- gene_name[,2]
####################################################################
# Do PCA
####################################################################
#principal component analysis

vst = vst(dds)

v <- plotPCA(vst, intgroup=c("condition"))
v<- v+ geom_label_repel(aes(label = name))
v
pcaData <- DESeq2::plotPCA(vst, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
#pdf("PCA_parents.pdf", height = 6, width = 6)
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=condition)) +
  geom_point(size=3) + geom_mark_ellipse(aes(fill=condition))+
  #scale_colour_manual(name="",values = c("a12"="goldenrod2", "gd33"="darkslateblue", "f1"="saddlebrown"))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  geom_label_repel(aes(label = name)) +
  coord_fixed()+theme_classic()
####################################################################
#Plotting Reps
####################################################################
plot_reps =  function(dds,x=1,y=2,cond_choice=1, cond='condition'){
  ##  Estimate the size factors for normalisation
  dds<-estimateSizeFactors(dds)
  
  ## Extract the normalised counts for the condition you want
  rep_values<- counts(dds, normalized=TRUE)[,dds[[cond]]==cond_choice]
  
  # Take logs of these values
  vals <- log2(rep_values[,c(x,y)] + 0.5)
  # And plot
  plot(vals,pch=16, cex=0.4,xlab=paste('rep',x),ylab=paste('rep',y))
  grid(col = "darkgray", lty = "solid",lwd = par("lwd"), equilogs = TRUE)
  title(paste("Comparison of",cond_choice,"replicates"))
}

par(mfrow = c(3,1))

plot_reps(dds, x=1, y=2, cond_choice="E")
plot_reps(dds, x=1, y=3, cond_choice="E")
plot_reps(dds, x=2, y=3, cond_choice="E")

plot_reps(dds, x=1, y=2, cond_choice="M")
plot_reps(dds, x=1, y=3, cond_choice="M")
plot_reps(dds, x=2, y=3, cond_choice="M")

####################################################################
#DEGs
####################################################################
filter_degs <- function(res){
  summary(res)
  res2 = res[!(is.na(res$padj)),]
  res2 = res2[res2$padj < 0.05,]
  return(res2)
}

resultsNames(dds)
E_M_DEGs = results(dds, contrast= c("condition", "E", "M"), alpha = 0.05, pAdjustMethod = "BH")
E_M_DEG = filter_degs(E_M_DEGs)

summary(E_M_DEG)
head(E_M_DEG)

write.table(rownames(E_M_DEG),"E_M_DEG",quote=F,row.names = F,col.names = F)

####################################################################
#Up and Downregulation
####################################################################

E_M_DEG_up <- E_M_DEG[E_M_DEG[,2]>0,]
E_M_DEG_down <- E_M_DEG[E_M_DEG[,2]<0,]

####################################################################
#MA Plots
####################################################################
par(mfrow = c(1,1))

DESeq2::plotMA(E_M_DEGs, ylim=c(-10,15), main='E_M_DEGs')

####################################################################
#Volcano Plots
####################################################################

library(EnhancedVolcano)

EnhancedVolcano(E_M_DEGs,
                lab = rownames(E_M_DEGs),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-5, 8),
                ylim = c(0,60))

####################################################################
#Heatmap
####################################################################

counts = counts(dds, normalized = TRUE)
counts <-  counts[apply(counts, MARGIN = 1, FUN = function(x) sd(x) != 0 ),]#it removes genes that are not express and have no variance
colnames(counts) <- c("E1","E2","E3","M1","M2","M3")

counts <- counts[rownames(counts) %in% rownames(E_M_DEG),]
pheatmap((log2(counts+1)), scale = "row",border_color=NA,show_rownames = F,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),main = 'DEGs expression across samples',cluster_rows = T, cluster_cols = T)

###################################################################
#TF terms
###################################################################

library(goseq)
library(tidyr)
library(dplyr)

TFs<- read.delim("./families_data.txt")

###################################################################
#Gage and list construction
###################################################################
library(gage)

ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory = htseqDir,design = ~  condition)

dds<-DESeq(ddsHTSeq)

#Exclude lowly expressed genes for GSEA
DESeq2_negative_gene_IDs <- is.na(as.data.frame(E_M_DEGs$log2FoldChange))

###################################################################
list <- list()
for(i in 1:52){
  
  TF_class <- as.character(unique(TFs$TF))
  TF_class_name <- TF_class[i]
  
  list[[i]] <- TFs[grep(paste(TF_class_name),TFs$TF),1]
}
names(list)<-TF_class[1:52]

#Run GAGE command for all leaky and induced expressed transgenics

Enriched <- gage(counts(dds)[!DESeq2_negative_gene_IDs,],list,ref=c(4:6),samp=c(1:3),
                 rank.test = T,
                 set.size=c(1,800), compare="unpaired",same.dir = T)

Enriched_greater <- Enriched$greater[1:51,1:5]
Enriched_lesser <- Enriched$less[1:51,1:5]

q.val <- -log10(Enriched_greater[,4])

data<-data.frame(rownames(Enriched_greater), q.val)
colnames(data) <- c("TF","q.val")


library(ggplot2)
ggplot(data[1:51,], aes(x=TF, y=q.val)) +geom_bar(stat="identity") +
  geom_col(aes(fill = q.val)) + 
  scale_fill_gradient2(low = "blue", 
                       high = "red", 
                       mid ="yellow",
                       midpoint = median(data$q.val)) +
  xlab("TF class") +
  ylab("-log10(q.val)") +
  ggtitle("Gene set enrichment analysis of TF classes upregulated") +
  theme_bw(base_size=10) + 
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=16, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=10, face="bold", hjust=1.10),
    axis.text.y=element_text(angle=0, size=10, face="bold", vjust=0.5),
    axis.title=element_text(size=12, face="bold"),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=14),  #Text size
    title=element_text(size=14)) +
  guides(colour=guide_legend(override.aes=list(size=2.5)))+
  geom_hline(yintercept=1.3,linetype="dashed", color = "red") +
  ylim(0,1.5)+
  coord_flip()

q.val <- -log10(Enriched_lesser[,4])

data<-data.frame(rownames(Enriched_lesser), q.val)
colnames(data) <- c("TF","q.val")


library(ggplot2)
ggplot(data[1:51,], aes(x=TF, y=q.val)) +geom_bar(stat="identity") +
  geom_col(aes(fill = q.val)) + 
  scale_fill_gradient2(low = "blue", 
                       high = "red", 
                       mid ="yellow",
                       midpoint = median(data$q.val)) +
  xlab("TF class") +
  ylab("-log10(q.val)") +
  ggtitle("Gene set enrichment analysis of TF classes downregulated") +
  theme_bw(base_size=10) + 
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=16, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=10, face="bold", hjust=1.10),
    axis.text.y=element_text(angle=0, size=10, face="bold", vjust=0.5),
    axis.title=element_text(size=12, face="bold"),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=14),  #Text size
    title=element_text(size=14)) +
  guides(colour=guide_legend(override.aes=list(size=2.5)))+
  geom_hline(yintercept=1.3,linetype="dashed", color = "red") +
  ylim(0,5)+
  coord_flip()

###################################################################

library("biomaRt")
library(topGO)
#collect gene names from biomart
mart <- biomaRt::useMart(biomart = "plants_mart",
                         dataset = "athaliana_eg_gene",
                         host = 'plants.ensembl.org')
# Get ensembl gene ids and GO terms
GTOGO <- biomaRt::getBM(attributes = c( "ensembl_gene_id",
                                        "go_id"), mart = mart)
#examine result
head (GTOGO)
#Remove blank entries
GTOGO <- GTOGO[GTOGO$go_id != '',]
# convert from table format to list format
geneID2GO <- by(GTOGO$go_id,
                GTOGO$ensembl_gene_id,
                function(x) as.character(x))
#examine result
head (geneID2GO)

all.genes <- sort(unique(as.character(GTOGO$ensembl_gene_id)))
int.genes <- rownames(E_M_DEG) # some random genes 
int.genes <- factor(as.integer(all.genes %in% int.genes))
names(int.genes) = all.genes

go.obj <- new("topGOdata", ontology='BP'
              , allGenes = int.genes
              , annot = annFUN.gene2GO
              , gene2GO = geneID2GO
)

resultsFisher <- runTest(go.obj, algorithm = "elim", statistic = "fisher")

allRes <- GenTable(go.obj, classic = resultsFisher,
            orderBy = "Fisher", ranksOf = "classic", topNodes = 17)

plot_go = function(goterms,name){
  goterms$percquery = allRes$Significant*100
  goterms$percback = allRes$Expected*100
  filtered_go = goterms[allRes$classic < 0.05,]
  #filtered_go = filtered_go[filtered_go$term_type == "P",]
  filtered_go_perc = cbind(filtered_go$percquery, filtered_go$percback)
  colnames(filtered_go_perc) = c("query","background")
  row.names(filtered_go_perc) = paste(filtered_go$Term,filtered_go$GO_acc,sep ="-->")
  meled = melt(filtered_go_perc)
  x = ggplot(meled, aes(Var1, value, fill=Var2)) +
    geom_bar(stat="identity", position="dodge")+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    xlab("sig GO Term") +
    ylab("Ratio genes with term in list") +
    ggtitle(name)+
    coord_flip()
  plot(x)
  return(x)
}

plot_go(allRes, "E vs M DEG GO terms")

showSigOfNodes(go.obj, score(resultsFisher), firstSigNodes = 20, useInfo = 'pval')
printGraph(go.obj, resultsFisher, firstSigNodes = 17, fn.prefix = "tGO", useInfo = "def", pdfSW = TRUE)

AgriGo <- read.delim('AgriGOv2_table.txt')
###################################################################
#go
####################################################################
plot_go = function(goterms,name){
  goterms$percquery = goterms$queryitem/goterms$querytotal*100
  goterms$percback = goterms$bgitem/goterms$bgtotal*100
  filtered_go = goterms[goterms$FDR < 0.05,]
  #filtered_go = filtered_go[filtered_go$term_type == "P",]
  filtered_go_perc = cbind(filtered_go$percquery, filtered_go$percback)
  colnames(filtered_go_perc) = c("query","background")
  row.names(filtered_go_perc) = paste(filtered_go$Term,filtered_go$GO_acc,sep ="-->")
  meled = melt(filtered_go_perc)
  x = ggplot(meled, aes(Var1, value, fill=Var2)) +
    geom_bar(stat="identity", position="dodge")+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    xlab("sig GO Term") +
    ylab("Ratio genes with term in list") +
    ggtitle(name)+
    coord_flip()
  plot(x)
  return(x)
}

v<-plot_go(AgriGo, "0")

plot_go(over_allRes[1:50,],"Top50 over-represented microspore GO term analysis 0.05")
plot_go(under_allRes[1:50,],"Top50 under-represented microspore GO term analysis 0.05")

