
setwd("~/Documents/Casparian Strip")
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
library(RColorBrewer)

sampleDataFilename <- 'sampleTable_no_blank_neg.txt'
sampleTable = read.table(sampleDataFilename,header=TRUE)
head(sampleTable)
htseqDir<-getwd()

##  Read in the results from the LibiNorm analysis (the counts files)
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory = htseqDir,design = ~condition)
## design<- you say to the test to do everything in relations to condition
## if you have more than one conditions you want to differentiate (for example different genotypes) you change design = ~  condition + genotype
##  And perform the analysis (details in the manual)

##  And perform the analysis (details in the manual)
dds<-DESeq(ddsHTSeq)

####################################################################
# Do PCA
####################################################################
#principal component analysis
vst = vst(dds,fitType="parametric")

v <- plotPCA(vst, intgroup=c("condition"),ntop=500)
v<- v+ geom_label_repel(aes(label = name))
v
pcaData <- DESeq2::plotPCA(vst, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
#pdf("PCA_parents.pdf", height = 6, width = 6)
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=condition)) +
  geom_point(size=3) + geom_mark_ellipse(aes(fill=condition,label = condition),con. = 0)+
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

plot_reps(dds, x=1, y=2, cond_choice="Col0")
plot_reps(dds, x=2, y=3, cond_choice="Col0")
plot_reps(dds, x=1, y=3, cond_choice="Col0")

par(mfrow = c(3,1))

plot_reps(dds, x=1, y=2, cond_choice="esb1")
plot_reps(dds, x=2, y=3, cond_choice="esb1")
plot_reps(dds, x=1, y=3, cond_choice="esb1")

par(mfrow = c(3,1))

plot_reps(dds, x=1, y=2, cond_choice="myb36")
plot_reps(dds, x=2, y=3, cond_choice="myb36")
plot_reps(dds, x=1, y=3, cond_choice="myb36")

par(mfrow = c(3,1))

plot_reps(dds, x=1, y=2, cond_choice="sgn3")
plot_reps(dds, x=2, y=3, cond_choice="sgn3")
plot_reps(dds, x=1, y=3, cond_choice="sgn3")

#dev.off()

####################################################################
#DEGs
####################################################################

##  This filter_degs function must be changed whenever the alpha value in the results() function changes as it filters out DEGs from based on p-value

filter_degs <- function(res){
  summary(res)
  res2 = res[!(is.na(res$padj)),]
  res2 = res2[res2$padj < 0.05,]
  return(res2)
}

##  alpha=p.value, lfcThreshold=log fold change, pAdjust..=hypothesis correction
##  only really appropriate to change lfcThreshold in my opinion, maybe pAdjust but p-value is just changing the rate of false positives

resultsNames(dds)
Esb1_DEGs = results(dds, contrast= c("condition", "esb1", "Col0"), alpha = 0.05, pAdjustMethod = "BH")
Esb1_DEG = filter_degs(Esb1_DEGs)

Sgn3_DEGs = results(dds, contrast= c("condition", "sgn3", "Col0"), alpha = 0.05, pAdjustMethod = "BH")
Sgn3_DEG = filter_degs(Sgn3_DEGs)

Myb36_DEGs = results(dds, contrast= c("condition", "myb36", "Col0"), alpha = 0.05, pAdjustMethod = "BH")
Myb36_DEG = filter_degs(Myb36_DEGs)

summary(TE_vs_NE_DEG)
summary(TE_vs_TM_DEG)
summary(TM_vs_NE_DEG)

Esb1_DEG_IDs <- rownames(Esb1_DEG)
Myb36_DEG_IDs <- rownames(Myb36_DEG)
Sgn3_DEG_IDs <- rownames(Sgn3_DEG)

Esb1_Myb36_DEG_IDs <- Esb1_DEG_IDs[Esb1_DEG_IDs %in% Myb36_DEG_IDs]
Shared_DEGs_IDs <- Esb1_Myb36_DEG_IDs[Esb1_Myb36_DEG_IDs %in% Sgn3_DEG_IDs]

####################################################################
#DEGs
####################################################################
par(mfrow=c(1,1))

log10.pval <- -log10(Esb1_DEGs$padj)
log2.fc    <- Esb1_DEGs$log2FoldChange
plot(log2.fc,log10.pval,
     xlab="log2 (fold change)",
     ylab="-log10 (p-value)",
     col="black",
     xlim=c(-10,10),
     ylim=c(0,20),
     main='Esb1_DEGs')
abline(h = -log10(0.05),col='red',lwd=1.5)
abline(v=-log2(2),col='blue',lwd=1.5)
abline(v=log2(2),col='blue',lwd=1.5)

log10.pval <- -log10(Myb36_DEGs$padj)
log2.fc    <- Myb36_DEGs$log2FoldChange
plot(log2.fc,log10.pval,
     xlab="log2 (fold change)",
     ylab="-log10 (p-value)",
     col="black",
     xlim=c(-10,10),
     ylim=c(0,20),
     main='Myb36_DEGs')
abline(h = -log10(0.05),col='red',lwd=1.5)
abline(v=-log2(2),col='blue',lwd=1.5)
abline(v=log2(2),col='blue',lwd=1.5)

log10.pval <- -log10(Sgn3_DEGs$padj)
log2.fc    <- Sgn3_DEGs$log2FoldChange
plot(log2.fc,log10.pval,
     xlab="log2 (fold change)",
     ylab="-log10 (p-value)",
     col="black",
     xlim=c(-10,10),
     ylim=c(0,20),
     main='Sgn3_DEGs')
abline(h = -log10(0.05),col='red',lwd=1.5)
abline(v=-log2(2),col='blue',lwd=1.5)
abline(v=log2(2),col='blue',lwd=1.5)

####################################################################
#Heatmap
####################################################################

counts = counts(dds , normalized = TRUE)

##  Removes genes that are not express and have no variance, will give error and heatmaps are meant to compare distinct expression profiles 

counts <-  counts[apply(counts, MARGIN = 1, FUN = function(x) sd(x) != 0 ),]
##  Select genes that we are interested in looking at across samples 
heatmap <- (counts[rownames(Esb1_DEG),])

pheatmap((log2(heatmap+1)), scale = "row",border_color=NA,show_rownames = F, main = 'Esb1_DEGs expression across samples',
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))

heatmap <- (counts[rownames(Sgn3_DEG),])

pheatmap((log2(heatmap+1)), scale = "row",border_color=NA,show_rownames = F, main = 'Sgn3_DEGs expression across samples',
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))

heatmap <- (counts[rownames(Myb36_DEG),])

pheatmap((log2(heatmap+1)), scale = "row",border_color=NA,show_rownames = F, main = 'Myb36_DEGs expression across samples',
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))

heatmap <- (counts[Shared_DEGs_IDs,])

pheatmap((log2(heatmap+1)), scale = "row",border_color=NA,show_rownames = F, main = 'Shared DEGs expression across samples',
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))

####################################################################
#Write counts
####################################################################
write.table(rownames(Esb1_DEG),'Esb1_DEGs.txt',quote=F,row.names = F,col.names = F)
write.table(rownames(Sgn3_DEG),'Sgn3_DEGs.txt',quote=F,row.names = F,col.names = F)
write.table(rownames(Myb36_DEG),'Myb36_DEGs.txt',quote=F,row.names = F,col.names = F)

####################################################################
#Random forest
####################################################################

