library(Seurat)
library(ggplot2)
library(patchwork)
library(ArchR)
library(cowplot)
library(Signac)
setwd("snATAC-seq/snATAC_cortex_TSS5_Frag1000/")
atac <- readRDS("Save-ArchR-Project.rds")
atac@sampleColData$ArrowFiles=dir("ArrowFiles")
atac@sampleColData$ArrowFiles=paste0('ArrowFiles/',atac@sampleColData$ArrowFiles)
names(atac@sampleColData$ArrowFiles)=substr(atac@sampleColData$ArrowFiles,1,nchar(atac@sampleColData$ArrowFiles)-6)
x=read.csv('~/b/cortex_merge_meta.csv')
ac=atac$cellNames
atac = subsetArchRProject(
  ArchRProj = atac,
  cells = ac[ac %in% x[,1]],
  outputDirectory = "ArchRSubset",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = FALSE
)
saveRDS(atac, "Save-ArchR-Project.rds")
rna <- readRDS("../../snRNA-seq/cortex_MTgene1.0%_UMI500_annotation.rds")
dd = rownames(rna@meta.data)
rna <- subset(rna, cells=dd[dd %in% x[,1]])
dd = rownames(rna@meta.data)
dm = x$subtype_SCT
names(dm) = x[,1]
rna@meta.data$subtype_SCT = dm[dd]
x=split(x,x$subtype_SCT)
ll=lapply(x,function(d){
    SimpleList(
        ATAC = ac[ac %in% d[,1]],
        RNA = dd[dd %in% d[,1]]
    )
  })
names(ll)=names(x)
groupList <- SimpleList(ll)
atac <- addGeneIntegrationMatrix(
    ArchRProj = atac, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = rna,
    addToArrow = FALSE, 
    groupList = groupList,
    groupRNA = "subtype_SCT"
)
