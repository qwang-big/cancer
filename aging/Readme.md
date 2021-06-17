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
```
