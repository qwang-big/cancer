x=read.csv('~/b/cortex_merge_meta.csv')
ac=atac$cellNames
dd=rownames(rna@meta.data)
atac = subsetArchRProject(
  ArchRProj = atac,
  cells = ac[ac %in% x[,1]],
  outputDirectory = "ArchRSubset",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = FALSE
)
rna <- readRDS("../../snRNA-seq/cortex_MTgene1.0%_UMI500_annotation.rds")
rna <- subset(rna, cells=dd[dd %in% x[,1]])
ll=lapply(split(x,x$subtype_SCT),function(d){
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
    groupRNA = "BioClassification"
)
