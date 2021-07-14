library(Seurat)
library(ArchR)
atac <- readRDS("Save-ArchR-Project.rds")
ff=dir("ArrowFiles")
atac@sampleColData <- DataFrame(row.names = substr(ff,1,nchar(ff)-6), ArrowFiles = paste0('ArrowFiles/',ff))
#cells <- unlist(lapply(atac@sampleColData$ArrowFiles,function(f) paste0(substr(f,1,nchar(f)-6),'#',h5read(f, "Metadata/CellNames"))))
#x=x[x[,1] %in% basename(cells),]
ac=getCellNames(atac)
x=read.csv('~/b/cortex_merge_meta.csv')
x=x[x[,1] %in% ac,]
ff=table(x$subtype_SCT)
ff=names(ff)[ff>100]
f="Astro_1"
#for (f in ff){
atac@projectMetadata$outputDirectory = paste0(f,'/tmp')
dir.create(paste0(f,'/tmp'))
atac = subsetArchRProject(
  ArchRProj = atac,
  cells = ac[ac %in% x[x$subtype_SCT==f,1]],
  outputDirectory = f,
  dropCells = TRUE,
  logFile = NULL,
  threads = 39,
  force = T
)
saveRDS(atac, paste0(f,"/Save-ArchR-Project.rds"))
}
library(Seurat)
rna <- readRDS("../../../snRNA-seq/cortex_MTgene1.0%_UMI500_annotation.rds")
x=read.csv('~/b/cortex_merge_meta.csv')
dd = rownames(rna@meta.data)
f='Astro_1'
r <- subset(rna, cells=dd[dd %in% x[x$subtype_SCT==f,1]])
saveRDS(r, paste0(f,'/seurat.rds'))
