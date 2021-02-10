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
  geom_point(size=3) +
  #scale_colour_manual(name="",values = c("a12"="goldenrod2", "gd33"="darkslateblue", "f1"="saddlebrown"))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+theme_bw()

####################################################################
# Do PCA
####################################################################

library(ggplot2)
library(ggrepel)

ggplot.PCA = function(data, PCx = 1, PCy = 2, title = "", square = F){
  # perform a PCA on the data
  pca <- prcomp(t(rlog(data))
                
                # the contribution to the total variance for each component
                percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
                d <- data.frame(PCa=pca$x[,PCx], PCb=pca$x[,PCy])
                
                pc_plot <- ggplot(data=d, aes(x=PCa, y=PCb)) + 
                  geom_hline(yintercept = 0, linetype="dashed") +
                  geom_vline(xintercept = 0, linetype="dashed") +
                  geom_point(size=2) +
                  geom_text_repel(aes(label = colnames(data)), data = d, force = 10) +
                  labs(title = title) +
                  xlab(paste0("PC",PCx,": ",round(percentVar[PCx] * 100),"% variance")) +
                  ylab(paste0("PC",PCy,": ",round(percentVar[PCy] * 100),"% variance")) +
                  theme_linedraw() +
                  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
                
                if(!square) { # if they should be proportional to variance
                  pc_plot = pc_plot +coord_fixed() 
                }
                
                return(pc_plot)
}

counts <- counts(dds)
ggplot.PCA(counts, 1,2)
ggplot.PCA(counts, 3,2)

## Modified plotPCA from DESeq2 package. Shows the Names of the Samples (the first col of SampleTable), and uses ggrepel pkg to plot them conveniently.
# @SA 10.02.2017 
dds<-dds[-5794,]

library(genefilter)
library(ggplot2)
library(ggrepel)

plotPCA.san <- function (object, intgroup = "condition", ntop = 1000, returnData = FALSE) 
{
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC3 = pca$x[,3], PC4 = pca$x[,4], group = colData(object)[[intgroup]], intgroup.df, name = colData(dds)[,1])
  if (returnData) {
    attr(d, "percentVar") <- percentVar[c(3,4)]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC3", y = "PC4", color = "group", label = "name")) + geom_point(size = 3) + xlab(paste0("PC3: ", round(percentVar[3] * 100), "% variance")) + ylab(paste0("PC4: ", round(percentVar[4] * 100), "% variance")) + coord_fixed() + geom_text_repel(size=3) 
  
}