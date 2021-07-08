library(Seurat)
library(ggplot2)
library(patchwork)
library(ArchR)
library(cowplot)
library(Signac)
setwd("/mnt/11/monkey_aging_brain/scATAC/cortex/cortex_TSS4_3000_OUT_NEW/") #将工作目录设置到ATAC文件夹下
atac <- readRDS("/mnt/11/monkey_aging_brain/scATAC/cortex/cortex_TSS4_3000_OUT_NEW/Save-ArchR-Project.rds")
atac@sampleColData$ArrowFiles=dir("ArrowFiles")
atac@sampleColData$ArrowFiles=paste0('ArrowFiles/',atac@sampleColData$ArrowFiles)
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
atac <- RunUMAP(atac, reduction = "lsi", dims = 2:30)


rna <- readRDS("/mnt/11/monkey_aging_brain/test_RNA_ATAC/cortex_test_RNA.rds")
atac$tech <- "rna"

atac$Sample <- atac$orig.ident
atac$group <- unlist(lapply(X = atac$Sample, FUN = function(x) {return(strsplit(as.character(x), split = "-")[[1]][[1]])}))
atac$organ <- "cortex"
sample_list <- unique(data.frame(age_group = as.character(rna$age_group),
  group = as.character(rna$group),
  age = as.character(rna$age)))
atac$age_group <- atac$group
atac$age_group <- sample_list$age_group[match(atac$group, sample_list$group)]
atac$group2 <- paste(atac$age_group, atac$group, sep = "_")
atac$sample <- paste(atac$group, atac$organ, sep = "_")
atac$age<- atac$group
atac$age <- sample_list$age[match(atac$group, sample_list$group)]


saveRDS(atac, "/mnt/11/monkey_aging_brain/test_RNA_ATAC/cortex_atac_signac_test.rds")





pdf("/mnt/11/monkey_aging_brain/test_RNA_ATAC/atac_rna.pdf")
p1 <- DimPlot(atac, reduction = "umap") + NoLegend() + ggtitle("scATAC-seq")
p2 <- DimPlot(rna, group.by = "subtype_SCT", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("scRNA-seq")
plot_grid(p1, p2)

transfer.anchors <- FindTransferAnchors(reference = rna, query = atac, features = VariableFeatures(object = rna),
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = rna$subtype_SCT,
                                     weight.reduction = atac[["lsi"]])

atac <- AddMetaData(atac, metadata = celltype.predictions)


genes.use <- VariableFeatures(rna)
refdata <- GetAssayData(rna, assay = "RNA", slot = "data")[genes.use, ]

imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = atac[["lsi"]])

atac[["RNA"]] <- imputation
coembed <- merge(x = rna, y = atac)
remove(rna,atac)
saveRDS(coembed, "cortex_ATAC_RNA_merge_log_raw.rds")
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)
coembed$subtype_SCT <- ifelse(!is.na(coembed$subtype_SCT), coembed$subtype_SCT, coembed$predicted.id)

saveRDS(coembed, "cortex_ATAC_RNA_merge_log_newest.rds")

path <- "/mnt/11/monkey_aging_brain/test_RNA_ATAC/"
p5 <- DimPlot(coembed, group.by = "tech", label.size = 5)
#ggsave(paste(path, "cortex_merge_plot_", "tech.png"), p5)
p6 <- DimPlot(coembed, group.by = "subtype_SCT", label = TRUE, repel = TRUE,label.size = 8)
#ggsave(paste(path, "cortex_merge_", "subtype_SCT.png"), p6)
#plot_grid(p5, p6)
p7 <- DimPlot(coembed, split.by = "tech", group.by = "subtype_SCT", label = TRUE)
p8 <- DimPlot(coembed, split.by = "tech", group.by = "age_group", label = TRUE)
#ggsave(paste(path, "cortex_merge_", "subtype_SCT.png"), p7)
#ggsave(paste(path, "cortex_merge_", "age_group.png"), p8)
p9 <- DimPlot(coembed, label = TRUE, label.size = 8)
#ggsave(paste(path, "cortex_merge_", "total.png"), p9)
p10 <- DimPlot(coembed, label = TRUE, label.size = 8, split.by = "group2")
#ggsave(paste(path, "cortex_merge_", "group2.png"), p10, width = 7, height = 80, limitsize = F)

dev.off()

#对RNA ATAC 分别画图
atac <- subset(rds, tech == "atac")

rna <- subset(rds, tech == "rna")

p11 <- DimPlot(rna, label = TRUE, label.size = 6, split.by = "group2")
ggsave(paste(path, "RNA_cortex_merge_", "group2.png"), p11, width = 80, height = 7, limitsize = F)
p12 <- DimPlot(atac, label = TRUE, label.size = 6, split.by = "group2")
ggsave(paste(path, "ATAC_cortex_merge_", "group2.png"), p12, width = 80, height = 7, limitsize = F)











