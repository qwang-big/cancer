suppressPackageStartupMessages({
library(RColorBrewer)
library(ArchR)
library(dplyr)
library(ggplot2)
library(BSgenome.Mfascicularis.NCBI.5.0)
library(GenomicFeatures)
library(OrganismDbi)

})


args <- dir()[grep("gz$", dir())]
name <- gsub("_scATAC.fragments.tsv.gz", "",args )

inputFiles <- args
names(inputFiles) <- name
names(inputFiles)

genomeAnnotation <- readRDS("/mnt/3/ywlai_genome/genome/Macaca/Macaca_genomeAnnotation")
geneAnnotation <- readRDS("/mnt/3/ywlai_genome/genome/Macaca/Macaca_geneAnnotation")


ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = name,
  filterTSS = 5, #Dont set this too high because you can always increase later
  filterFrags = 3000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation
)


doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1
)

  proj <- ArchRProject(
    ArrowFiles = ArrowFiles, 
    outputDirectory = "./striatum_scATAC_20210217",
    copyArrows = T,
    geneAnnotation = geneAnnotation,
    genomeAnnotation = genomeAnnotation
  )


  proj <- filterDoublets(proj, filterRatio = 2)

#proj$group <- gsub("_pigv2.fragments_new", "", proj$Sample)

name <- "thulums_scATAC_20210308"

p1 <- plotGroups(
    ArchRProj = proj, groupBy = "Sample",
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "ridges"
   )
ggsave(paste(name, ".qc1.png", sep =""), p1)

p2 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
ggsave(paste(name, ".qc2.png", sep =""), p2)

p3 <- plotFragmentSizes(ArchRProj = proj)
ggsave(paste(name, ".qc3.png", sep =""), p3)
proj$log10_uniqueFrags <- log10(proj$nMonoFrags)
p4 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10_uniqueFrags",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
ggsave(paste(name, ".qc4.png", sep =""), p4)


proj <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI1", 
    iterations = 5, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30
)

proj <- addUMAP(
    ArchRProj = proj, 
    reducedDims = "IterativeLSI1", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine", force = T
)
  
proj <- addClusters(
    input = proj,
    reducedDims = "IterativeLSI1",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8
)

proj <- addImputeWeights(
  ArchRProj = proj,
  reducedDims = "IterativeLSI1",
  dimsToUse = NULL,
  scaleDims = NULL,
  corCutOff = 0.75,
  td = 3,
  ka = 4,
  sampleCells = 5000,
  nRep = 2,
  k = 15,
  epsilon = 1,
  useHdf5 = TRUE,
  randomSuffix = FALSE,
  threads = getArchRThreads(),
  seed = 1,
  verbose = TRUE,
  logFile = createLogFile("addImputeWeights")
)

seRNA <- readRDS("/mnt/11/monkey_aging_brain/striatum/annotation/striatum_annotation.rds")
seRNA@active.assay <- "RNA"
proj@reducedDims$IterativeLSI <- proj@reducedDims$IterativeLSI1

proj <- addGeneIntegrationMatrix(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    reduction = "cca",
    addToArrow = FALSE,
    sampleCellsATAC = 10000,
    groupRNA = "celltype",
    nameCell = "Cell_celltype",
    nameGroup = "Group_celltype",
    nameScore = "Score_celltype"
)

proj$group <- unlist(lapply(X = proj$Sample, FUN = function(x) {return(strsplit(x, split = "-")[[1]][[1]])}))
proj$organ <- "cortex"
sample_list <- unique(data.frame(age_group = as.character(seRNA$age_group),
  group = as.character(seRNA$group),
  age = as.character(seRNA$age)))
proj$age_group <- proj$group
proj$age_group <- sample_list$age_group[match(proj$group, sample_list$group)]



umap_plot1 <- plotEmbedding(proj, colorBy = "cellColData", name = "group",  embedding = "UMAP")
ggsave(paste(path,"group.umap1_harmony.png", sep = ""), umap_plot1)
umap_plot2 <- plotEmbedding(proj, colorBy = "cellColData", name = "age_group",  embedding = "UMAP")
ggsave(paste(path,"age_group.umap2_harmony.png", sep = ""), umap_plot2)
umap_plot3 <- plotEmbedding(proj, colorBy = "cellColData", name = "Clusters",  embedding = "UMAP")
ggsave(paste(path,"Clusters.umap3_harmony.png", sep = ""), umap_plot3)
umap_plot4 <- plotEmbedding(proj, colorBy = "cellColData", name = "Group_celltype",  embedding = "UMAP")
ggsave(paste(name,".umap4_harmony.png", sep = ""), umap_plot4)


saveArchRProject(ArchRProj = proj, outputDirectory = name, load = T, overwrite = T)

gene <- c("SLC17A7", "NFEH", "SLC1A2",  # excitory
		 "GRM4", "RBFOX3", "FAT2", # granule
		 "GAD1", # inhibitory
		 "FLT1", "DUSP1", # endo
		 "MOG",  ## Oligodendrocyte
		 "NEU4", "CSPG4", "PDGFRA", # Opc, oligodendrocyte percusor
	 	 "PDGFRB", "COBLL1",   # pericyte
	 	 "APBB1IP", "P2RY12", # Microglia
	      "SLC4A4", "GPC5", "GRIA1", "SLC1A3", # Astrocyte
	      "ALDH1A1",  # Astrocyte_cerebellar
	      "RYR1", "GLCE", "SORCS3" # Purkinje neurons
	      )  

gene <- c("PROX1", "PID1", "PCP4", "SULF2", "SLC1A3", "CSF1R",
 "SLC38A11", "NR2F2", "LHX6", "MOG", "KCNIP4", "RYR3", "RXFP1", "PPFIA2", "ABCA12")

for (f in gene){
    try(p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = f, 
    size = 1,
    embedding = "UMAP",
    quantCut = c(0.01, 0.95),
    imputeWeights = NULL
))
    try(ggsave(paste(name,f,".png", sep = ""), p))
}

p <- plotBrowserTrack(
    ArchRProj = proj, 
    groupBy = "Group_celltype", 
    geneSymbol = gene, 
    upstream = 50000,
    downstream = 50000
)


grid::grid.newpage()
for (i in seq_along(1:length(gene))){
pdf(paste(name, gene[i],"browserplot.pdf", sep = "_"), width=12, height=6)
  grid::grid.draw(p[[i]])
dev.off()
}

## peaks
pathToMacs2 <- findMacs2()
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "Group_celltype", force = T)
proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "Group_celltype", 
    geneAnnotation= geneAnnotation,
    pathToMacs2 = pathToMacs2,
    genomeAnnotation = genomeAnnotation,
    genomeSize = 2.5e9,
    force = TRUE
)
proj <- addPeakMatrix(proj)



proj <- addBgdPeaks(proj,force = T)

#proj <- addMotifAnnotations(ArchRProj = proj, force = T, species = getGenome(proj))
proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "JASPAR2020", name = "JASPAR2020_Motif", force = T)
#proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "homer", name = "homer_Motif", force = T)
#proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "encode", name = "encode_Motif", force = T)


#for (g in c("JASPAR2020_Motif", "homer_Motif", "encode_Motif")){
for (g in c("JASPAR2020_Motif")){
matrixname <- paste(g, "dev",sep = "_")
proj <- addDeviationsMatrix(
  ArchRProj = proj,
  peakAnnotation = g,
  matches = NULL,
  matrixName = matrixname,
  out = c("z", "deviations"),
  binarize = FALSE,
  threads = getArchRThreads(),
  verbose = TRUE,
  parallelParam = NULL,
  logFile = createLogFile("addDeviationsMatrix"),
  force = T
)
}

#plot Heatmap


proj$annotation <- proj$predictedGroup_Un
proj$annotation <- unlist(lapply(X = proj$annotation, FUN = function(x) {return(strsplit(x, split = "_")[[1]][[1]])}))
proj$group <- proj$Sample
proj$group <- unlist(lapply(X = proj$group, FUN = function(x) {return(strsplit(x, split = "_")[[1]][[1]])}))

sample_group <- unlist(tapply(rownames(proj), proj$annotation,function(i)sample(i,100)))
g="homer_Motif"


  matrixname <- paste(g, "dev",sep = "_")
  cisbp_matrix <- getMatrixFromProject(
    ArchRProj = subsetCells(ArchRProj = proj, cellNames = sample_group),
    useMatrix = matrixname,
    useSeqnames = NULL,
    verbose = TRUE,
    binarize = FALSE,
    threads = getArchRThreads(),
    logFile = createLogFile("getMatrixFromProject")
  )

  cisbp_meta <- cisbp_matrix@colData
  cisbp_matrix_de <- cisbp_matrix@assays@data$deviations
  sd <- rowSds(as.matrix(cisbp_matrix_de))  
  cisbp_matrix_de <- cisbp_matrix_de[order(sd, decreasing= T),]
  cisbp_matrix_de <- cisbp_matrix_de[1:100,]


  ## plot cell type

  cisbp_meta$new_ident <- paste(cisbp_meta$group, cisbp_meta$annotation, sep = "_")
  cisbp_meta <- cisbp_meta[order(cisbp_meta$new_ident),]     
  order_rn <- rownames(cisbp_meta)
  order_rn_match <- match(order_rn , colnames(cisbp_matrix_de))
  cisbp_matrix_de <- cisbp_matrix_de[,order_rn_match]


  col_fun = viridis(length(unique(cisbp_meta$new_ident)))
  names(col_fun) = unique(cisbp_meta$new_ident)

  col_fun2 = brewer.pal(n=3, "Set2")
  col_fun2 = col_fun2[1:2]
  names(col_fun2) = unique(cisbp_meta$group)

  ha1 = HeatmapAnnotation(
      cell_type = cisbp_meta$new_ident,
      cell_group = cisbp_meta$group,
      col = list(cell_type = col_fun,
        cell_group = col_fun2
      ),
      na_col = "black"
  )

  hm_motif1 <-  Heatmap(as.matrix(cisbp_matrix_de), 
       cluster_columns = F, cluster_rows = T,
      show_column_names = FALSE, top_annotation = ha1)
  hm_motif2 <-  Heatmap(as.matrix(cisbp_matrix_de), 
       cluster_columns = T, cluster_rows = T,
      show_column_names = FALSE, top_annotation = ha1)


  png(paste(Sys.time(), "hmplot", g, "motif.png", sep = "_"), width = 2400, height = 2300)
  hm_motif1
  dev.off()
  png(paste(Sys.time(), "hmplot", g, "motif_cluster.png", sep = "_"), width = 2400, height = 2300)
  hm_motif2
  dev.off()







