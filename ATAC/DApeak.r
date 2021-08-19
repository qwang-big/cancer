x=read.csv('~/b/hippo_ATAC_cellColdata.csv')
rownames(x)=x[,1]
atac@cellColData$celltype=x[rownames(atac@cellColData),'merge_subtype_SCT']
ff=dir("ArrowFiles")
atac@sampleColData <- DataFrame(row.names = substr(ff,1,nchar(ff)-6), ArrowFiles = paste0('ArrowFiles/',ff))
dd=unique(atac@cellColData$celltype)
atac@cellColData$age_group=as.character(atac@cellColData$age_group)
atac@projectMetadata$outputDirectory='./tmp'
dir.create('./tmp')
library(BSgenome.Mfascicularis.NCBI.5.0)
seqlevels(BSgenome.Mfascicularis.NCBI.5.0) = sub('MFA','chr',seqlevels(BSgenome.Mfascicularis.NCBI.5.0))
atac <- addGroupCoverages(ArchRProj = atac, groupBy = "age_group")
atac <- addReproduciblePeakSet(
  ArchRProj = atac, groupBy = 'age_group', peakMethod = 'Macs2', genomeSize=2.85e9,
  pathToMacs2 = '/hwfssz5/ST_PRECISION/OCG/wangqi/miniconda3/envs/py36/bin/macs2'
)
atac=addPeakMatrix(atac)
for (d in dd[-1]) {
atac <- readRDS("Save-ArchR-Project.rds")
atac@cellColData$celltype=x[rownames(atac@cellColData),'merge_subtype_SCT']
ff=dir("ArrowFiles")
atac@sampleColData <- DataFrame(row.names = substr(ff,1,nchar(ff)-6), ArrowFiles = paste0('ArrowFiles/',ff))
dd=unique(atac@cellColData$celltype)
ac=getCellNames(atac)
a = subsetArchRProject(
  ArchRProj = atac,
  cells = ac[ac %in% x[x$merge_subtype_SCT==d,1]],
  outputDirectory = d,
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = T
)
markersPeaks <- getMarkerFeatures(
    ArchRProj = a, 
    useMatrix = "PeakMatrix", 
    groupBy = "age_group",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  normBy='nFrags', testMethod = "wilcoxon"
)
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1.5")
for(i in 1:4)
write.table(markerList[[i]],file=paste0(d,'/',i,'.bed'),sep='\t',quote=F,col.names=F)
}
for (d in dd) {
setwd(paste0('../',d))
atac <- readRDS("Save-ArchR-Project.rds")
ff=dir("ArrowFiles")
atac@sampleColData <- DataFrame(row.names = substr(ff,1,nchar(ff)-6), ArrowFiles = paste0('ArrowFiles/',ff))
atac@cellColData$age_group=as.factor(atac@cellColData$age_group)
markersPeaks <- getMarkerFeatures(
    ArchRProj = atac, 
    useMatrix = "PeakMatrix", 
    groupBy = "age_group",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 1")
for(i in 1:4)
write.table(markerList[[i]],file=paste0(i,'.bed'),sep='\t',quote=F,col.names=F)
}
