library(Seurat)
args <- commandArgs(trailingOnly = TRUE)
n=6
ff=dir('..')
lw=function(d)length(which(d))
df=data.frame(V1=ff[grep('rds$',ff)])
df$f=unlist(lapply(df[,1],function(d) substr(d,1,13)))
dd=split(df,df$f)
lapply(dd, function(d){tryCatch({
load(paste0('../',substr(d[1,1],1,4),'.rda'))
Tumor=unlist(lapply(d[,1], function(s){
cnv_obj=readRDS(paste0('../',s))
cl = cutree(cnv_obj$cluster$hclust, k=n)
unique(unlist(lapply(seq(1,max(cl)),function(i) {
tb=cl[cl==i]
if (lw(cnv_obj$reference_obs %in% names(tb)) <100 & length(tb)>1000){
return(names(tb))}else{return(c())}
})))
}))
tt=table(Tumor)
Tumor=names(tt[tt>floor(nrow(d)/3)])
Tumor=unique(Tumor)
Tumor = Tumor[gsub('_.*','',Tumor) == 'Tumor']
print(length(Tumor))
seurat_spatialObj@meta.data$cnv='low'
seurat_spatialObj@meta.data[gsub('.*_','',Tumor),]$cnv='high'
p=SpatialDimPlot(seurat_spatialObj,group.by = 'cnv',pt.size.factor = 2.5)
png(paste0('~/b/pic/cnv/',gsub('rda','vs',d[1,2]),'.png'), width=10*max(abs(p$data$imagerow)), height=10*max(abs(p$data$imagecol)))
print(p)
dev.off()}, error=function(e){})
})
