#' Differential gene expression analysis pipeline of oxLDL-treated M2 macrophages vs untreated M2 macrophages


# Clear environment
rm(list = ls())


#### Libraries ####

# ## Data wrangling
# library(tidyverse)         
# 
# ## DGE analysis
library(DESeq2)             
library(factoextra)       # Screeplots (PCA) 
# library(ashr)             # LFC shrinkage  
# 
# ## Visualisation
library(ggplot2)        # main package for visualisation
library(vsn)            # for meanSdPlot()
library(pheatmap)       # Heatmaps
# library(gridExtra)        # to arrange heatmaps in a grid. Alternative to patchwork
# library(ggrepel)          # Labeling      
# library(RColorBrewer)       
# library(viridis)          # colour scheme (incl Inferno)  
# library(VennDiagram)
# library(patchwork)        # display of multiple plots and simple arrangements: p1 + p2
# library(scales)           # for scaling values (used for visualisations, etc.)
# 
# 
# ## Annotation
# #library(tximport)         # For transcript database?
# library(GenomicFeatures)  # Import GTF file data for annotation 
# library(biomaRt)          # Main annotation package
# #library(AnnotationDbi)    # Annotation - mapIds() function  # replace by biomart



##### Utility functions ####







#### Loading data ####
#' Processing of RNAseq data prior to import into DESeq2 using STAR 2.7.3a and featureCounts (from subread 2.0.2). 
#' No trimming was performed. QC performed by novogone

# Metadata/ Sample Table 
colData <- data.frame(
  sample = factor(c('A03_M2', 'A03_MOX', 
                    'A04_M2', 'A04_MOX',
                    'S02_M2', 'S02_MOX')),
  condition = factor(rep(c('M2',"MOX"), 3)),
  donor = factor(c('A03', 'A03', 
                   'A04', 'A04', 
                   'S02', 'S02'))
)
rownames(colData) <- colData$sample 

# Count matrix
cts <- as.matrix(read.csv("../data/countmatrix.txt", sep="\t",row.names="Geneid"))
colnames(cts) <- colData$sample # change the long sample names to short sample names

print(paste("Full data set:", all(rownames(colData) %in% colnames(cts)))) # Checking if the sampleTable and count matrix annotations correspond to each other



#### Running differential expression analysis
# Parsing count matrix to DESeq2 dataset
ddsCts <- DESeqDataSetFromMatrix(countData = cts,
                                 colData = colData,
                                 design = ~donor + condition)
# Running DESeq analysis
dds <- DESeq(ddsCts)


#### Assessing data
# Transformation for data set: Variance stabilized transformation
vsd <- vst(dds, blind=TRUE)
vsd_mat <- assay(vsd)


# Per gene standard deviation
# Plotting the per-gene standard deviation (SD, taken across samples) against the rank of the mean for the the variance-stabilizing transformation(vst)
vsd_msd <- meanSdPlot(vsd_mat, ylab="Standard deviation (vst)", plot=FALSE)

png(paste0("../Results/",Sys.Date(),"_VST-Transformation.png"), width=400, height=350 )
vsd_msd$gg + ggtitle("Per gene standard deviation (VST)")
dev.off()


# Hierarchical clustering of samples
vsd_cor  <- cor(vsd_mat)#Computing the pairwise correlation values for transformed counts.
colnames(vsd_cor)
annotation <- data.frame(Subset = c(rep(c("M2","MOX"), 3) ))# Preparing annotations and labeling.
rownames(annotation) <- colnames(vsd_cor) #all matrices have the same rownames. Note that the rownames are colnames in the correlation matrix
ann_colors = list(Subset = c(M2 ="#ebb573", MOX ="#eb9282"))

pheatmap(vsd_cor, annotation = annotation, annotation_colors = ann_colors, main = "VST transformed",  silent=F)


# Principal component analysis
pca <- prcomp(t(vsd_mat)) # compute PCA
fviz_screeplot(pca, addlabels = T) + ylim(0,100) + theme_classic() + labs(title = "Variances - PCA (VST transformated counts)\n(Full)", x = "Principal components", y = "% of variance") # Screeplot visualizes the variances accross the different PCs



df.pca <- cbind(colData, pca$x)# Create data frame with metadata form the sampleTable, this contains the different PCs. This is used to plot the PCA

# get % of Variance
var.pc1 <- round(summary(pca)$importance[2,1]*100 , 1) # % variance of PC1
var.pc2 <- round(summary(pca)$importance[2,2]*100 , 1) # % variance of PC2

# PC1 and PC2 according to subset
pca_plot1 <- ggplot(df.pca, aes(x=PC1,y=PC2, color = condition))
p1 <- pca_plot1 + labs(color = "Condition") + geom_point(size=5) + ggtitle("PCA by subset") +
  scale_color_manual(values=c("#8cafc8", "#eb9282")) +
  labs(fill = "Condition") +
  xlab(paste0("PC1 (", var.pc1,"% variance)")) +    
  ylab(paste0("PC2 (", var.pc2,"% variance)")) +   
  coord_fixed() +
  xlim(-100,100) +
  ylim(-100,100) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p1

# PC1 and PC2 according to donor
pca_plot2 <- ggplot(df.pca, aes(x=PC1,y=PC2, color = donor))
p2 <- pca_plot2 + labs(color = "Donor") + geom_point(size=5) + ggtitle("PCA by donor") +
  labs(fill = "Donor") +
  xlab(paste0("PC1 (", var.pc1,"% variance)")) +    
  ylab(paste0("PC2 (", var.pc2,"% variance)")) +  
  coord_fixed() +
  xlim(-100,100) +
  ylim(-100,100) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p2
