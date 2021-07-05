library(SpatialExperiment)
library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)
library(Seurat)
library(cowplot)
seurat_spatialObj@meta.data$row = seurat_spatialObj@images$slice1@coordinates$row
seurat_spatialObj@meta.data$col = seurat_spatialObj@images$slice1@coordinates$col
sce <- as.SingleCellExperiment(seurat_spatialObj)
sce <- scater::logNormCounts(sce)
sce <- spatialPreprocess(sce, platform="ST", 
                              n.PCs=7, n.HVGs=2000, log.normalize=FALSE)
sce <- spatialCluster(sce, q=10, platform="ST", d=7,
                           init.method="mclust", model="t", gamma=2,
                           nrep=10000, burn.in=100, save.chain=TRUE)
p3 <-SpatialDimPlot(seo, label = T, label.size = 3)+ labs(title="Seurat") +DarkTheme()
p4 <- clusterPlot(sce,label = "spatial.cluster" )+ labs(title="BayesSpace")+DarkTheme()
plot_grid(p3, p4)
