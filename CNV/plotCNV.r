library(Seurat)
args <- commandArgs(trailingOnly = TRUE)
n=args[1]
ff=dir()
lw=function(d)length(which(d))
df=data.frame(V1=ff[grep('rds$',ff)])
df$f=unlist(lapply(df[,1],function(d) substr(d,1,13)))
dd=split(df,df$f)
dd=dd[1]
lapply(dd, function(d){tryCatch({
load(paste0(substr(d[1,1],1,4),'.rda'))
Tumor=unlist(lapply(d[,1], function(s){
cnv_obj=readRDS(s)
tree = as.data.frame(cutree(cnv_obj$cluster$hclust,k = n))
colnames(tree) = 'Cluster'
tree$barcode = rownames(tree)
tree = tree[cnv_obj$target,]
i=which.min(lapply(seq(1,max(tree$Cluster)),function(i)lw(tree$Cluster==i)))
Tumor = tree[tree$Cluster==i,]$barcode
}))
tt=table(Tumor)
Tumor=names(tt[tt>floor(nrow(d)/3)])
Tumor=unique(Tumor)
Tumor = Tumor[gsub('_.*','',Tumor) == 'Tumor']
print(length(Tumor))
seurat_spatialObj@meta.data$cnv='low'
seurat_spatialObj@meta.data[gsub('.*_','',Tumor),]$cnv='high'
png(paste0(d[1,2],'.png'))
p=SpatialDimPlot(seurat_spatialObj,group.by = 'cnv',pt.size.factor = 2.5)
print(p)
dev.off()}, error=function(e){})
})

