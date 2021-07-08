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

exe=function(f, nrep=10000){
load(f)
seurat_spatialObj@meta.data$row = seurat_spatialObj@images$slice1@coordinates$row
seurat_spatialObj@meta.data$col = seurat_spatialObj@images$slice1@coordinates$col
sce <- as.SingleCellExperiment(seurat_spatialObj)
sce <- scater::logNormCounts(sce)
sce <- spatialPreprocess(sce, platform="Visium", n.PCs=7, n.HVGs=2000, log.normalize=FALSE)
q = length(unique(seurat_spatialObj@meta.data$seurat_cluster)) - 1
sce <- spatialCluster(sce, q=q, platform="Visium", d=7, init.method="mclust", model="t", gamma=2, nrep=nrep, burn.in=100, save.chain=TRUE)
p3 <-SpatialDimPlot(seurat_spatialObj, label = T, label.size = 3) +DarkTheme()
p4 <- clusterPlot(sce,label = "spatial.cluster",color=NA)+ DarkTheme()+scale_y_continuous(trans="reverse")
p=cowplot::plot_grid(p3, p4)
ggplot2::ggsave(p,file=paste0("/home/wangqi9/tmp/",f,'.png'),width=20,height=10)
c1=scran::findMarkers(sce, seurat_spatialObj@meta.data$seurat_cluster, pval.type="all")
c2=scran::findMarkers(sce, colData(sce)$spatial.cluster, pval.type="all")
saveRDS(sce,file=paste0("/home/wangqi9/tmp/",f,'.rds'))
d1=plyr::ldply(lapply(seq_len(q),function(j)rownames(c1[[j]][c1[[j]]$p.value < 0.05,])),rbind)
d2=plyr::ldply(lapply(seq_len(q),function(j)rownames(c2[[j]][c2[[j]]$p.value < 0.05,])),rbind)
write.table(d1,sep='\t',file=paste0("/home/wangqi9/tmp/",f,'.csv'),col.names=F)
write.table(d2,sep='\t',file=paste0("/home/wangqi9/tmp/",f,'.csv'),col.names=F)
}
