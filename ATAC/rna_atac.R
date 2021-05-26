library(Seurat)
library(ggplot2)
library(patchwork)
library(ArchR)
library(cowplot)
library(Signac)

atac <- readRDS("/hwfssz5/ST_PRECISION/TOMCAT/xuyoujia/Macaca_eye/scATAC/test/M1_panc/Save-projMonkey4/Save-ArchR-Project.rds")
atac@sampleColData$ArrowFiles=dir("ArrowFiles")
names(atac@sampleColData$ArrowFiles)=substr(atac@sampleColData$ArrowFiles,1,nchar(atac@sampleColData$ArrowFiles)-6)

peaks <- getMatrixFromProject(atac, useMatrix='PeakMatrix')
genescore <- getMatrixFromProject(atac, useMatrix='GeneScoreMatrix')
peak_grange <- getPeakSet(atac)
a <- peak_grange
a$rowData <- paste(a@seqnames,a@ranges, sep = '-')

counts <- assays(peaks)$PeakMatrix
rownames(counts) <- a$rowData
atac <- CreateSeuratObject(counts = counts, assay = "ATAC", project = "ATAC")

gs <- assays(genescore)$GeneScoreMatrix
#rowname <- paste0(rowData(genescore)$seqnames,rowData(genescore)$start,"-",rowData(genescore)$end)
rowname <- rowData(genescore)$name
rownames(gs) <- rowname
atac[["ACTIVITY"]] <- CreateAssayObject(counts = gs)
atac$tech <- "atac"

DefaultAssay(atac) <- "ACTIVITY"
atac <- FindVariableFeatures(atac)
atac <- NormalizeData(atac)
atac <- ScaleData(atac)
DefaultAssay(atac) <- "ATAC"
VariableFeatures(atac) <- names(which(Matrix::rowSums(atac) > 100))
atac <- RunLSI(atac, n = 50, scale.max = NULL)
atac <- RunUMAP(atac, reduction = "lsi", dims = 1:50)

all <-readRDS("/hwfssz5/ST_PRECISION/TOMCAT/xuyoujia/Macaca_eye/iDrop_scRNA_Macaque_IN_EX/rds/M1_panc_mt5_celltype.rds")
rna <- CreateSeuratObject(counts = all@assays$RNA@counts, project = "rna", meta.data = all@meta.data)
rna$tech = "rna"
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(rna)
rna <- ScaleData(rna, features = all.genes)
rna <- RunPCA(rna)
rna <- FindNeighbors(rna, dims = 1:10)
rna <- FindClusters(rna,resolution = 0.5)
rna <- RunUMAP(rna, dims = 1:10)
saveRDS(rna, "/hwfssz5/ST_PRECISION/TOMCAT/xuyoujia/Macaca_eye/scATAC/test/M1_panc/M1_panc_rna_umap.rds")

pdf("/hwfssz5/ST_PRECISION/TOMCAT/xuyoujia/Macaca_eye/scATAC/test/M1_panc/atac_rna.pdf", width = 10, height = 5)
p1 <- DimPlot(atac, reduction = "umap") + NoLegend() + ggtitle("scATAC-seq")
p2 <- DimPlot(rna, group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("scRNA-seq")
plot_grid(p1, p2)

transfer.anchors <- FindTransferAnchors(reference = rna, query = atac, features = VariableFeatures(object = rna),
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = rna$celltype,
                                     weight.reduction = atac[["lsi"]])
atac <- AddMetaData(atac, metadata = celltype.predictions)


genes.use <- VariableFeatures(rna)
refdata <- GetAssayData(rna, assay = "RNA", slot = "data")[genes.use, ]

imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = atac[["lsi"]])

atac[["RNA"]] <- imputation
coembed <- merge(x = rna, y = atac)

coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)
coembed$celltype <- ifelse(!is.na(coembed$celltype), coembed$celltype, coembed$predicted.id)

p5 <- DimPlot(coembed, group.by = "tech")
p6 <- DimPlot(coembed, group.by = "celltype", label = TRUE, repel = TRUE)
plot_grid(p5, p6)
DimPlot(coembed, split.by = "tech", group.by = "celltype", label = TRUE)
dev.off()


