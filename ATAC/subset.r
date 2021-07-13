library(Seurat)
library(ArchR)
setwd("snATAC-seq/snATAC_cortex_TSS5_Frag1000/")
atac <- readRDS("Save-ArchR-Project.rds")
atac@sampleColData$ArrowFiles=dir("ArrowFiles")
atac@sampleColData$ArrowFiles=paste0('ArrowFiles/',atac@sampleColData$ArrowFiles)
names(atac@sampleColData$ArrowFiles)=substr(atac@sampleColData$ArrowFiles,1,nchar(atac@sampleColData$ArrowFiles)-6)
x=read.csv('~/b/cortex_merge_meta.csv')
ff=unique(x$subtype_SCT)
ac=getCellNames(atac)
for(f in ff[2:28]){
atac1 = subsetArchRProject(
  ArchRProj = atac,
  cells = ac[ac %in% x[x$subtype_SCT==f,1]],
  outputDirectory = f,
  dropCells = TRUE,
  logFile = NULL,
  threads = 6,
  force = T
)
saveRDS(atac1, paste0(f,"/Save-ArchR-Project.rds"))
}
