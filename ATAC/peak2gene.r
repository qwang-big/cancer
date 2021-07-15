library(Seurat)
library(ggplot2)
library(patchwork)
library(ArchR)
library(cowplot)
library(Signac)
setwd("snATAC-seq/snATAC_cortex_TSS5_Frag1000/")
atac <- readRDS("Save-ArchR-Project.rds")
ff=dir("ArrowFiles")
atac@sampleColData <- DataFrame(row.names = substr(ff,1,nchar(ff)-6), ArrowFiles = paste0('ArrowFiles/',ff))
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
    reducedDims = "IterativeLSI1",
    seRNA = rna,
    addToArrow = T, force = TRUE,
    groupList = NULL,logFile='/dev/null',
    groupRNA = "subtype_SCT",threads = 1
)
atac <- addPeak2GeneLinks(atac, reducedDims = "IterativeLSI1")
p2g <- getPeak2GeneLinks(
    ArchRProj = atac,
    corCutOff = 0.45,
    resolution = 1e4,
    returnLoops = T
)
mat = plotPeak2GeneHeatmap(atac, groupBy = "age_group", returnMatrices=T)
png(paste0('~/b/pic/Astro_1.png'))
plotPeak2GeneHeatmap(atac, groupBy = "age_group")
dev.off()
markerGenes <- mat[['Peak2GeneLinks']][mat[['ATAC']][['kmeansId']]==17,'gene']
p <- plotBrowserTrack(
    ArchRProj = atac, 
    groupBy = "age_group", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000,
    loops = p2g
)
for(f in markerGenes){
png(paste0('~/b/pic/',f,'.png'),width=800,height=800)
grid::grid.newpage()
grid::grid.draw(p[[f]])
dev.off()
}
pk = mat[['Peak2GeneLinks']][mat[['ATAC']][['kmeansId']]==17,]
gr = read.table(text=gsub(':','-',pk$peak), sep='-', , col.names=c("chr", "start", "end"))
gr$gene = pk$gene
gr = makeGRangesFromDataFrame(gr, keep.extra.columns=TRUE)
table(queryHits(findOverlaps(gr, p2g[[1]])))
