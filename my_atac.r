library(reshape2)
library(Seurat)
library(Signac)
library(Matrix)
library(ggplot2)
library(ArchR)
fname1 = function(d,n=4) substr(d,1,nchar(d)-n)
rna=readRDS('../../20200120_Final_1217_Neocortex.rds')

load('genescore.rda')
gs <- assays(genescore)$GeneScoreMatrix
rownames(gs) <- rowData(genescore)$name
#gs = gs[,seq(1,ncol(gs),4)]
atac <- CreateSeuratObject(counts = gs)
atac <- NormalizeData(atac)
atac <- FindVariableFeatures(atac)
atac <- ScaleData(atac)
#atac <- FindTopFeatures(atac, min.cutoff = "q0")
#VariableFeatures(atac) <- names(which(Matrix::rowSums(atac) > 100))
atac <- RunPCA(atac)
atac <- RunUMAP(atac, reduction = "pca", dims = 1:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
rna <- ScaleData(rna)

anchors <- FindTransferAnchors(reference = rna, query = atac)
predictions <- TransferData(anchorset = anchors, refdata = rna$Celltype, k.weight=10)
atac <- AddMetaData(object = atac, metadata = predictions)
save(atac,file=paste0('all.rda'))
p1 <- DimPlot(rna, group.by = "Celltype", label = TRUE) + NoLegend() + ggtitle("RNA")
p2 <- DimPlot(atac, group.by = "predicted.id", label = TRUE) + NoLegend() + ggtitle('ATAC')
p=p1+p2
ggsave(p,file=paste0('all.pdf'))

