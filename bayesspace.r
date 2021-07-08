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
sig=c('C1QTNF6','GTSE1','STRIP2','FAM72B','RS1','MKI67','SULT1A1','C11orf16','SHCBP1','LOXL2','AGER','C2orf91',
'TLR10','MS4A1','P2RX2','CA4','FCRL2','GPR78','DONSON','TNFSF11','VIPR1','TNFRSF13C','FDCSP','RSPH9','PI15',
'EMP2','NDC80','KIF4A','FCRLA','ADAM12','CD52','COL11A1','CYB5A','AC025594.2','SYT2','PTPRN','GRAMD2A','BLK',
'SLC6A8','E2F8','SPP1','PIFO','SEC14L4','PAX5','COL5A2','TUBB3','CHAD','CD19','DOK2','CILP2','IBSP','FFAR4',
'SLC1A1','AK4','SERPINH1','AGMAT','PITX1','VPREB3','COL1A1','ANKRD66','TCTE1','MMP14','KNTC1','CASP14','TRPA1')
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
df=do.call(rbind,lapply(seq_len(length(c2)),function(i){
d=data.frame(fc=-abs(c2[[i]][rownames(c2[[i]]) %in% sig,'summary.logFC']),cl=i,method='bayes')
rbind(d,data.frame(fc=-abs(c1[[i]][rownames(c1[[i]]) %in% sig,'summary.logFC']),cl=i,method='seurat'))
}))
p=ggplot(df, aes(fc,color=method)) + stat_ecdf(geom = "step")+facet_wrap(~cl)+labs(title=f,x='fold change rank',y='ecdf')
ggplot2::ggsave(p,file=paste0("/home/wangqi9/tmp/",f,'-cl.png'),width=8,height=6)
}
