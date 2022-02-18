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

makePCAplot <- function(mat, sampleTable, PC.x = 1, PC.y = 2, colorby="Condition", main="PCA-plot"){
  #' Perform and plot PCA of count matrix.
  #' Requires transformed count matrix (e.g. VST transformation)
  #' @param mat Array, 2D matrix with GeneIDs as rows and samples as columns 
  #' @param sampleTable DataFrame, Table containing the sample information/meta data of the sequencing. 
  #' @param PC.x Numeric, Principal component on x-axis. [Default: 1]
  #' @param PC.y Numeric, Principal component on y-axis. [Default: 2]
  #' @param colorby Character, [Default: "condition"]
  #' @param main Character, title of plot. [Default: "PCA-plot"]
  #' @examples 
  #' makePCAplot(vsd_matrix, PC.x = 2, PC.y = 3, colorby = "treatment", main = "PCA plot by Treatment")
  
  # Prepare PCA data
  pca <- prcomp(t(mat)) # compute PCA
  df.pca <- cbind(sampleTable, pca$x)
  
  
  # get proportion of Variance
  var.pcx <- round(summary(pca)$importance[2,PC.x]*100 , 1) # % variance of PCx
  var.pcy <- round(summary(pca)$importance[2,PC.y]*100 , 1) # % variance of PCy
  
  # plot PCA results
  PC.x = paste0("PC", PC.x)
  PC.y = paste0("PC", PC.y)
  pca_plot <- ggplot(df.pca, aes_string(x = PC.x, y = PC.y, color= colorby))
  pca_plot + labs(color = colorby) + geom_point(size=5) + ggtitle(main) +
    labs(fill = colorby) +
    xlab(paste0(PC.x, " (", var.pcx,"% variance)")) +    
    ylab(paste0(PC.y ," (", var.pcy,"% variance)")) +  
    coord_fixed() +
    xlim(-100,100) +
    ylim(-100,100) +
    theme_bw() +
    scale_color_manual(values=c("#8cafc8", "#eb9282", "#ebb573")) +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
}





#### Loading data ####
#' Processing of RNAseq data prior to import into DESeq2 using STAR 2.7.3a and featureCounts (from subread 2.0.2). 
#' No trimming was performed. QC performed by novogone

# Metadata/ Sample Table 
colData <- data.frame(
  Sample = factor(c('A03_M2', 'A03_MOX', 
                    'A04_M2', 'A04_MOX',
                    'S02_M2', 'S02_MOX')),
  Condition = factor(rep(c('M2',"MOX"), 3)),
  Donor = factor(c('A03', 'A03', 
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
                                 design = ~Donor + Condition)
# Running DESeq analysis
dds <- DESeq(ddsCts)


#### Assessing data
## Transformation for data set: Variance stabilized transformation
vsd <- vst(dds, blind=TRUE)
vsd_mat <- assay(vsd)


## Per gene standard deviation
# Plotting the per-gene standard deviation (SD, taken across samples) against the rank of the mean for the the variance-stabilizing transformation(vst)
vsd_msd <- meanSdPlot(vsd_mat, ylab="Standard deviation (vst)", plot=FALSE)

png(paste0("../Results/",Sys.Date(),"_VST-Transformation.png"), width=400, height=350 )
vsd_msd$gg + ggtitle("Per gene standard deviation (VST)")
dev.off()


## Hierarchical clustering of samples
vsd_cor  <- cor(vsd_mat)#Computing the pairwise correlation values for transformed counts.
colnames(vsd_cor)
annotation <- data.frame(Subset = c(rep(c("M2","MOX"), 3) ))# Preparing annotations and labeling.
rownames(annotation) <- colnames(vsd_cor) #all matrices have the same rownames. Note that the rownames are colnames in the correlation matrix
ann_colors = list(Subset = c(M2 ="#ebb573", MOX ="#eb9282"))

pheatmap(vsd_cor, annotation = annotation, annotation_colors = ann_colors, main = "VST transformed",  silent=F)


## Principal component analysis
pca <- prcomp(t(vsd_mat)) # compute PCA
fviz_screeplot(pca, addlabels = T) + ylim(0,100) + theme_classic() + labs(title = "Variances - PCA (VST transformated counts)\n(Full)", x = "Principal components", y = "% of variance") # Screeplot visualizes the variances accross the different PCs


# Plots of genes contributing to variance for the individual PCs
# PC1
fviz_contrib(pca, choice = "var", axes = 1, top = 15) + 
  theme_classic() + 
  labs(title = "Contribution to PC1", x = "Genes", y = "Contribution (%)") +
  theme(axis.text.x = element_text(angle = 90))

# PC2
fviz_contrib(pca, choice = "var", axes = 2, top = 15) + 
  theme_classic() + 
  labs(title = "Contribution to PC2", x = "Genes", y = "Contribution (%)") +
  theme(axis.text.x = element_text(angle = 90))

# PC3
fviz_contrib(pca, choice = "var", axes = 3, top = 15) + 
  theme_classic() + 
  labs(title = "Contribution to PC3", x = "Genes", y = "Contribution (%)") +
  theme(axis.text.x = element_text(angle = 90))


# PCA plots
png(paste0("../Results/",Sys.Date(),"PCA_PC1+PC2_condition.png"), width=300, height=275 )
makePCAplot(mat=vsd_mat, sampleTable = colData, colorby="Condition", PC.x = 1, PC.y = 2, main = "PC1 and PC2 by condition")
dev.off()

png(paste0("../Results/",Sys.Date(),"PCA_PC1+PC2_donor.png"), width=300, height=275 )
makePCAplot(mat=vsd_mat, sampleTable = colData, colorby="Donor", PC.x = 1, PC.y = 2, main = "PC1 and PC2 by donor")
dev.off()

png(paste0("../Results/",Sys.Date(),"PCA_PC1+PC3_condition.png"), width=300, height=275 )
makePCAplot(mat=vsd_mat, sampleTable = colData, colorby="Condition", PC.x = 1, PC.y = 3, main = "PC1 and PC3 by condition")
dev.off()

png(paste0("../Results/",Sys.Date(),"PCA_PC1+PC3_donor.png"), width=300, height=275 )
makePCAplot(mat=vsd_mat, sampleTable = colData, colorby="Donor", PC.x = 1, PC.y = 3, main = "PC1 and PC3 by donor")
dev.off()




## Diagnostic plots
# Cooks distances
png(paste0("../Results/",Sys.Date(),"_Cooks-distance.png"), width=300, height=275 )
par(mar=c(5,5,1,1))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=1, main="Cook's distance")
dev.off()

# Dispersion estimates
png(paste0("../Results/",Sys.Date(),"_Dispersion-estimates.png"), width=300, height=275 )
par(mar=c(5,5,2,2))
plotDispEsts(dds) 
title("Dispersion estimates")
dev.off()


