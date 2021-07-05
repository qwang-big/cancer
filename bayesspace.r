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
p3 <-SpatialDimPlot(seurat_spatialObj, label = T, label.size = 3) +DarkTheme()
p4 <- clusterPlot(sce,label = "spatial.cluster",color=NA)+ DarkTheme()
cowplot::plot_grid(p3, p4)

c1=scran::findMarkers(sce, seurat_spatialObj@meta.data$seurat_cluster, pval.type="all")
c2=scran::findMarkers(sce, colData(sce)$spatial.cluster, pval.type="all")
n=1:10
df=data.frame(lapply(n,function(i){
x=rownames(c1[[i]][c1[[i]]$p.value < 0.05,])
unlist(lapply(n,function(j){
y=rownames(c2[[j]][c2[[j]]$p.value < 0.05,])
length(intersect(x,y))/length(union(x,y))
}))
}))
colnames(df)=n
heatmap.2(as.matrix(df),Rowv=F,Colv=F,trace='none')
