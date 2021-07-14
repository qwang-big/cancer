library(Seurat)
library(ArchR)
setwd("snATAC-seq/snATAC_cortex_TSS5_Frag1000/")
atac <- readRDS("Save-ArchR-Project.rds")
atac@sampleColData$ArrowFiles=dir("ArrowFiles")
atac@sampleColData$ArrowFiles=paste0('ArrowFiles/',atac@sampleColData$ArrowFiles)
names(atac@sampleColData$ArrowFiles)=substr(atac@sampleColData$ArrowFiles,1,nchar(atac@sampleColData$ArrowFiles)-6)
s="/hwfssz1-tmp/ST_PRECISION/USER/wangqi/mk_brain/snATAC-seq/snATAC_cortex_TSS5_Frag1000/ArrowFiles/"
ff=dir(s)
cells <- unlist(lapply(ff,function(f) paste0(substr(f,1,nchar(f)-6),'#',h5read(paste0(s,f), "Metadata/CellNames"))))
ac=getCellNames(atac)
x=read.csv('~/b/cortex_merge_meta.csv')
x=x[x[,1] %in% cells,]
ff=table(x$subtype_SCT)
ff=names(ff)[ff>100]
for(f in ff[2:28]){
atac@projectMetadata$outputDirectory = paste0(f,'/tmp')
dir.create(paste0(f,'/tmp'))
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
library(Seurat)
rna <- readRDS("../../../snRNA-seq/cortex_MTgene1.0%_UMI500_annotation.rds")
x=read.csv('~/b/cortex_merge_meta.csv')
dd = rownames(rna@meta.data)
f='Astro_1'
r <- subset(rna, cells=dd[dd %in% x[x$subtype_SCT==f,1]])
saveRDS(r, paste0(f,'/seurat.rds'))
