## Process seurat obj
```r
pbmc=readRDS('striatum_Astrocytes.rds')
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
colnames(pbmc@assays$SCT@data) = pbmc@meta.data[colnames(pbmc@assays$SCT@data),'age_group']
nm=pbmc@meta.data$age_group
unm=unique(nm)
qs=lapply(unm,function(n){
apply(pbmc@assays$SCT@data[,which(nm==n)],1,quantile,seq(0,1,.1))
})
names(qs)=unm
m=do.call(rbind,lapply(qs,function(d) d[8,]))
mm=t(m[,apply(m,2,function(d) length(which(d!=0))>=3)])
write.table(mm[,c(3:1,4)],file='mm.csv',sep=',',quote=F)

pbmc=readRDS("cortex.Microglia")
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap",group.by='age_group')
pbmc=RunTSNE(pbmc, dims = 1:20)
DimPlot(pbmc, reduction = "tsne",group.by='age_group')
pbmc=RunTSNE(pbmc, features=ll)
DimPlot(pbmc, reduction = "tsne",group.by='age_group')
ll=toupper(ll[which(toupper(ll) %in% rownames(pbmc@assays$RNA@data))])

for( d in levels(pbmc@meta.data$celltype)){
  x=subset(x = pbmc, subset =celltype== d)
  saveRDS(x,file=paste0('/home/wangqi9/b/data/monkey_aging/hippocampus.',d))
}

colnames(pbmc@assays$RNA@data) = pbmc@meta.data[colnames(pbmc@assays$RNA@data),'sample']
nm=pbmc@meta.data$sample
unm=unique(nm)
qs=do.call(cbind,lapply(unm,function(n){
apply(pbmc@assays$RNA@data[,which(nm==n)],1,sum)
}))

df2n=function(d){r=d[,2];names(r)=d[,1];r}
info=unique(pbmc@meta.data[c('sample','age_group')])
info=df2n(info)
cl=RColorBrewer::brewer.pal(7, "BrBG")
names(cl)=sort(unique(info))
d=t(apply(qs[ll,],1,function(x) (x-min(x))/(max(x)-min(x))))
colnames(d)=unm
gplots::heatmap.2(na.omit(d),ColSideColors=cl[info[unm]],col=bluered(11),trace='none',srtCol=45)
r=c("Monkey3b_cortex", "Monkey5g_cortex","Monkey5f_cortex")
rmv = pbmc@meta.data[pbmc@meta.data$sample %in% r,]
pbmc <- pbmc[,!colnames(pbmc) %in% rownames(rmv)]

```
