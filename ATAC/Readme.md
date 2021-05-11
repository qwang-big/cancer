## ArchR file format
```py
def h5py_dataset_iterator(g, prefix=''):
    for key in g.keys():
        item = g[key]
        path = '{}/{}'.format(prefix, key)
        if isinstance(item, h5py.Dataset): # test for dataset
            yield (path, item)
        elif isinstance(item, h5py.Group): # test for group (go down)
            yield from h5py_dataset_iterator(item, path)


f= h5py.File(s, 'r')
for (path, dset) in h5py_dataset_iterator(f):
	print(path, dset)
```
## parse Arrow file
```
library(Seurat)
library(GenomicRanges)
library(ArchR)
library(Signac)
setwd('/home/wangqi9/mnt/data/monkey/Save-ProjMonkey2/ArrowFiles/')
ff=dir()
ff=ff[grep('arrow$',ff)]
for (ArrowFile in ff){
for (i in c(1:20,'X')) {
chr=paste0('chr',i)
fname1 = function(d,n=6) substr(d,1,nchar(d)-n)
output <- h5read(ArrowFile, paste0("Fragments/",chr,"/Ranges")) %>% 
      {IRanges(start = .[,1], width = .[,2])}
#sampleName <- h5read(ArrowFile, paste0("Metadata/Sample"))

mcols(output)$RG <- Rle(
      values = h5read(ArrowFile, paste0("Fragments/",chr,"/RGValues")), 
      lengths = h5read(ArrowFile, paste0("Fragments/",chr,"/RGLengths"))
    )

h5closeAll()

output <- GRanges(seqnames = i, ranges(output), RG = mcols(output)$RG) 

write.table(output, file=paste(fname1(ArrowFile),chr,"bed",sep='.'), quote=F, sep="\t", row.names=F, col.names=F)
}
}
```
## ATAC peaks to gene features
```
zcat Macaca_fascicularis.Macaca_fascicularis_5.0.102.chr.gtf.gz |g protein_cod|awk '$3=="transcript"'|perl -F'\t' -lne '$F[1]=$F[6] eq "-"?($F[4]-100)."\t".($F[4]+1000):($F[3]-1000)."\t".($F[3]+100);$F[2]=$1 if /gene_name \"(\S+)\"/;print join "\t", @F[0..2]'|g -v transcript > Macaca_fascicularis_prom.bed
#for f in *.arrow;do cat ${f%.*}.chr*.bed|sort -k1,1 -k2,2n > ${f%.*}.bed;done

#! /bin/sh

/home/wangqi9/d/bin/bedtools2/bin/bedtools intersect -a $1 -b /home/wangqi9/d/data/monkey/Save-ProjMonkey2/Macaca_fascicularis_prom.bed -wo |cut -f6,10|sort|uniq -c > $1.ctx

qsub -clear -cwd -l vf=5g,p=1 -binding linear:1 -q st.q -V "Neocortex_Monkey2-19.9_20200819_AY12.chr11.bed" -P P20Z10200N0059 bedt.sh 

for f in *.arrow;do cat ${f%.*}.chr*.bed.ctx > ${f%.*}.mtx;done
cat *.mtx|perl -F'\s+' -lne '$F[2]=$1 if /BC(\d+)/;print join "\t", @F[1..3]' > all.mat
```
## scATAC and scRNA-seq integration
```
library(reshape2)
library(Seurat)
library(Signac)
library(Matrix)
library(ggplot2)
fname1 = function(d,n=4) substr(d,1,nchar(d)-n)
rna=readRDS('../../20200120_Final_1217_Neocortex.rds')

ff=dir('.','*.mtx$')
for(f in ff){
x=read.table(f)
f=fname1(f)
x1=dcast(x,V2~V3, value.var = 'V1', fun.aggregate = sum)
m = Matrix(as.matrix(x1[,-1]),sparse=TRUE)
rownames(m)=x1[,1]
atac <- CreateSeuratObject(counts = t(m))
atac <- NormalizeData(atac)
atac <- FindVariableFeatures(atac)
atac <- ScaleData(atac)
#atac <- FindTopFeatures(atac, min.cutoff = "q0")
atac <- RunPCA(atac)
atac <- RunUMAP(atac, reduction = "pca", dims = 1:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

#rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna)
#rna <- ScaleData(rna)

anchors <- FindTransferAnchors(reference = rna, query = atac)
predictions <- TransferData(anchorset = anchors, refdata = rna$Celltype, k.weight=10)
atac <- AddMetaData(object = atac, metadata = predictions)
save(atac,file=paste0(f,'.rda'))
p1 <- DimPlot(rna, group.by = "Celltype", label = TRUE) + NoLegend() + ggtitle("RNA")
p2 <- DimPlot(atac, group.by = "predicted.id", label = TRUE) + NoLegend() + ggtitle(f)
p=p1+p2
ggsave(p,file=paste0(f,'.pdf'))
}
```
