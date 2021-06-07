### Get the parameters
parser = argparse::ArgumentParser(description="Script to QC scATAC data by ArchR")
parser$add_argument('-I','--inputFiles', help='input fragment file')
parser$add_argument('-D','--id', help='sample ID')
parser$add_argument('-T','--filtertss', help='tss threhold')
parser$add_argument('-F','--filterFrags', help='fragment number threhold')
parser$add_argument('-O','--out', help='out directory')
parser$add_argument('-B','--bsgenome', help='BSgenome')
parser$add_argument('-G','--txDb', help='TxDb')
parser$add_argument('-R','--org', help='org')
args = parser$parse_args()

###

library("ArchR")
library("GenomicRanges")


addArchRThreads(threads = 1)
######Creating Arrow Files

library(args$bsgenome,character.only = T)
library(args$txDb,character.only = T)
library(args$org,character.only = T)

r1=get(args$txDb)
r2=get(args$org)

genomeAnnotation <- createGenomeAnnotation(genome = args$bsgenome,filter = TRUE)
geneAnnotation <- createGeneAnnotation(TxDb = r1, OrgDb = r2)

#################################################################################################################
#####------------------------------------------------ArchR--------------------------------------------------#####
#################################################################################################################

inputFiles=c(args$inputFiles)
names(inputFiles)=c(args$id)

setwd(args$out)

#########-----------------------------Creating Arrow Files---------------------------------#########
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = as.numeric(args$filtertss), #Dont set this too high because you can always increase later
  filterFrags = as.numeric(args$filterFrags),
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  force=T,
  nChunk=1,
  subThreading=FALSE
)
#########-----------------------------Inferring scATAC-seq Doublets with ArchR-------------#########

addArchRThreads(threads = 1)

doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1
)
### Get the parameters
parser = argparse::ArgumentParser(description="Script to Cluster scATAC data by ArchR")
parser$add_argument('-I','--inputpath', help='input arrow directory')
parser$add_argument('-D','--id', help='tissue ID')
parser$add_argument('-O','--out', help='out directory')
parser$add_argument('-F','--filterRatio', help='filterRatio')
parser$add_argument('-T','--threads', help='threads')
parser$add_argument('-B','--bsgenome', help='BSgenome')
parser$add_argument('-G','--txDb', help='TxDb')
parser$add_argument('-R','--org', help='org')
args = parser$parse_args()

###
library("ArchR")
library("GenomicRanges")

addArchRThreads(threads = as.numeric(args$threads))
######Creating Arrow Files

library(args$bsgenome,character.only = T)
library(args$txDb,character.only = T)
library(args$org,character.only = T)

r1=get(args$txDb)
r2=get(args$org)

genomeAnnotation <- createGenomeAnnotation(genome = args$bsgenome,filter = TRUE)
geneAnnotation <- createGeneAnnotation(TxDb = r1, OrgDb = r2)

#################################################################################################################
#####------------------------------------------------ArchR--------------------------------------------------#####
#################################################################################################################

#########-----------------------------Creating Arrow Files---------------------------------#########
setwd(args$inputpath)
ArrowFiles <- list.files(pattern = args$id)
ArrowFiles_1=paste(args$inputpath,ArrowFiles,sep="/")

setwd(args$out)

########------------------------------Creating an ArchRProject----------------------------########
proj1 <- ArchRProject(
  ArrowFiles = ArrowFiles_1, 
  outputDirectory = "Results",
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  copyArrows = FALSE
)
saveArchRProject(ArchRProj = proj1, outputDirectory = "Save-Proj1", load = FALSE)
########-----------------------------Filtering Doublets from an ArchRProject-------------########
proj2 <- filterDoublets(proj1,filterRatio = as.numeric(args$filterRatio))
########-----------------------------Iterative Latent Semantic Indexing (LSI)------------########
proj2 <- addIterativeLSI(
    ArchRProj = proj2,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30
)
########-----------------------------Batch Effect Correction wtih Harmony----------------########
proj2 <- addHarmony(
    ArchRProj = proj2,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample"
)
#######-----------------------------Clustering using Seurat’s FindClusters() function-----#######
proj2 <- addClusters(
    input = proj2,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8
)
######------------------------------Single-cell Embeddings--------------------------------#######
proj2 <- addUMAP(
    ArchRProj = proj2, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)

proj2 <- addUMAP(
    ArchRProj = proj2, 
    reducedDims = "Harmony", 
    name = "UMAPHarmony", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)

p1 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p4 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")

plotPDF(p1,p2,p3,p4, name = paste(args$id,"Plot-UMAP2Harmony-Sample-Clusters.pdf",sep="_"), ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)
######------------------------------Marker Genes Imputation with MAGIC--------------------------------#######
proj2 <- addImputeWeights(proj2)
saveArchRProject(ArchRProj = proj2, outputDirectory = "Save-Proj2", load = FALSE)
######------------------------------Identifying Marker Genes--------------------------------#######
markersGS <- getMarkerFeatures(
  ArchRProj = proj2, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 0")

ml=as.data.frame(markerList)
write.table(ml,paste(args$id,"marker_list.xls",sep="_"),sep = "\t",quote = FALSE,row.names = FALSE)
######------------------------------Marker Genes Imputation with MAGIC--------------------------------#######
proj2 <- addImputeWeights(proj2)
saveArchRProject(ArchRProj = proj2, outputDirectory = "Save-Proj2", load = FALSE)




library(ArchR)

proj <- loadArchRProject(path = "/hwfssz5/ST_PRECISION/OCG/wuliang2_SingleCell/CellOmics/Project/Rat_Cortex_ATAC/1.Cluster/Cortex_F2/Save-Proj2_F2")

cellColData=as.data.frame(proj@cellColData)
coor=proj@embeddings$UMAPHarmony$df
cellColData=cbind(cellColData,coor)
cellColData$Celltype="N"
cellColData$Celltype[cellColData$Clusters=="C1"]="Astrocyte"
cellColData$Celltype[cellColData$Clusters=="C2"]="Astrocyte"
cellColData$Celltype[cellColData$Clusters=="C3"]="Astrocyte"
cellColData$Celltype[cellColData$Clusters=="C4"]="Oligodendrocyte_precursor"
cellColData$Celltype[cellColData$Clusters=="C5"]="Oligodendrocyte"
cellColData$Celltype[cellColData$Clusters=="C6"]="Excitatory_neuron_1"
cellColData$Celltype[cellColData$Clusters=="C7"]="Excitatory_neuron_2"
cellColData$Celltype[cellColData$Clusters=="C8"]="Excitatory_neuron_3"
cellColData$Celltype[cellColData$Clusters=="C9"]="Excitatory_neuron_3"
cellColData$Celltype[cellColData$Clusters=="C10"]="Excitatory_neuron_3"
cellColData$Celltype[cellColData$Clusters=="C11"]="Excitatory_neuron_3"
cellColData$Celltype[cellColData$Clusters=="C12"]="Excitatory_neuron_4"
cellColData$Celltype[cellColData$Clusters=="C13"]="Excitatory_neuron_4"
cellColData$Celltype[cellColData$Clusters=="C14"]="Excitatory_neuron_4"
cellColData$Celltype[cellColData$Clusters=="C15"]="PVALB"
cellColData$Celltype[cellColData$Clusters=="C16"]="SST"
cellColData$Celltype[cellColData$Clusters=="C17"]="Excitatory_neuron_5"
cellColData$Celltype[cellColData$Clusters=="C18"]="Inhibitory_neuron"
cellColData$Celltype[cellColData$Clusters=="C19"]="VIP"
cellColData$Celltype[cellColData$Clusters=="C20"]="Excitatory_neuron_5"
cellColData$Celltype[cellColData$Clusters=="C21"]="Microglia"
cellColData$Celltype[cellColData$Clusters=="C22"]="Microglia"
cellColData$Celltype[cellColData$Clusters=="C23"]="Endothelial_cell"
cellColData$Celltype[cellColData$Clusters=="C24"]="Meningeal_cell"
cellColData$Celltype[cellColData$Clusters=="C25"]="Pericytes"
proj$Celltype=cellColData$Celltype

###

pathToMacs2 <- findMacs2()
proj1 <- addGroupCoverages(ArchRProj = proj, groupBy = "Celltype")
projHeme4 <- addReproduciblePeakSet(
    ArchRProj = proj1, 
    groupBy = "Celltype", 
    pathToMacs2 = pathToMacs2,
    threads = 5,
    genomeSize = 2782028915
)

projHeme5 <- addPeakMatrix(projHeme4)

saveArchRProject(ArchRProj = projHeme5, outputDirectory = "/hwfssz5/ST_PRECISION/OCG/wuliang2_SingleCell/CellOmics/Project/Rat_Cortex_ATAC/2.Callpeak/Save-Proj-MACS2", load = FALSE)

getAvailableMatrices(projHeme5)

markersPeaks <- getMarkerFeatures(
    ArchRProj = projHeme5, 
    useMatrix = "PeakMatrix", 
    groupBy = "Celltype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

peakmarkerList <- getMarkers(markersPeaks, cutOff = "FDR <= 1")
pml= as.data.frame(peakmarkerList)
write.table(pml,paste0("/hwfssz5/ST_PRECISION/OCG/wuliang2_SingleCell/CellOmics/Project/Rat_Cortex_ATAC/2.Callpeak/","Peak_marker_list_1007.xls"),sep="\t",quote=F)



library(ArchR)

proj <- loadArchRProject(path = "/hwfssz5/ST_PRECISION/OCG/wuliang2_SingleCell/CellOmics/Project/Rat_Cortex_ATAC/3.Motif/Save-Proj-MACS2-Motif")

motifPositions <- getPositions(proj)

motifs <- c("NEUROD2", "ASCL1", "SOX9", "LHX4", "SPI1")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
projHeme5 <- addGroupCoverages(ArchRProj = proj, groupBy = "Celltype",force = TRUE)

seFoot <- getFootprints(
  ArchRProj = projHeme5, 
  positions = motifPositions[markerMotifs], 
  groupBy = "Celltype",
  flank = 100
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = projHeme5, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5,
  flank = 100,
  height=12,
  width=8
)

###----------------------library packages---------------------###
library(ArchR)
library(ggplot2)
library(gridExtra)
library(Rtsne)
library(clusterProfiler)
###----------------------load data----------------------------###
proj <- loadArchRProject(path = "Save-Proj2_F2")

cellColData=as.data.frame(proj@cellColData)
coor=proj@embeddings$UMAPHarmony$df
cellColData=cbind(cellColData,coor)
cellColData$Celltype="N"
cellColData$Celltype[cellColData$Clusters=="C1"]="Astrocyte"
cellColData$Celltype[cellColData$Clusters=="C2"]="Astrocyte"
cellColData$Celltype[cellColData$Clusters=="C3"]="Astrocyte"
cellColData$Celltype[cellColData$Clusters=="C4"]="Oligodendrocyte_precursor"
cellColData$Celltype[cellColData$Clusters=="C5"]="Oligodendrocyte"
cellColData$Celltype[cellColData$Clusters=="C6"]="Excitatory_neuron_1"
cellColData$Celltype[cellColData$Clusters=="C7"]="Excitatory_neuron_2"
cellColData$Celltype[cellColData$Clusters=="C8"]="Excitatory_neuron_3"
cellColData$Celltype[cellColData$Clusters=="C9"]="Excitatory_neuron_3"
cellColData$Celltype[cellColData$Clusters=="C10"]="Excitatory_neuron_3"
cellColData$Celltype[cellColData$Clusters=="C11"]="Excitatory_neuron_3"
cellColData$Celltype[cellColData$Clusters=="C12"]="Excitatory_neuron_4"
cellColData$Celltype[cellColData$Clusters=="C13"]="Excitatory_neuron_4"
cellColData$Celltype[cellColData$Clusters=="C14"]="Excitatory_neuron_4"
cellColData$Celltype[cellColData$Clusters=="C15"]="PVALB"
cellColData$Celltype[cellColData$Clusters=="C16"]="SST"
cellColData$Celltype[cellColData$Clusters=="C17"]="Excitatory_neuron_5"
cellColData$Celltype[cellColData$Clusters=="C18"]="Inhibitory_neuron"
cellColData$Celltype[cellColData$Clusters=="C19"]="VIP"
cellColData$Celltype[cellColData$Clusters=="C20"]="Excitatory_neuron_5"
cellColData$Celltype[cellColData$Clusters=="C21"]="Microglia"
cellColData$Celltype[cellColData$Clusters=="C22"]="Microglia"
cellColData$Celltype[cellColData$Clusters=="C23"]="Endothelial_cell"
cellColData$Celltype[cellColData$Clusters=="C24"]="Meningeal_cell"
cellColData$Celltype[cellColData$Clusters=="C25"]="Pericytes"
proj$Celltype=cellColData$Celltype

GeneScoreMatrix <- getMatrixFromProject(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix"
)
gene_score <- assays(GeneScoreMatrix)$GeneScoreMatrix
rownames(gene_score) <- rowData(GeneScoreMatrix)$name

###------------------------------Figure2--------------------------------------###
cellColData$Sample_1=unlist(lapply(strsplit(as.character(cellColData$Sample),"_CL"),"[",1))
cellColData$Region=unlist(lapply(strsplit(as.character(cellColData$Sample),"_2020"),"[",1))
cellColData$Sample_1=gsub("20200611_","",cellColData$Sample_1)
cellColData$Sample_1=gsub("AN11","1",cellColData$Sample_1)
cellColData$Sample_1=gsub("AN12","2",cellColData$Sample_1)
cellColData$Sample_1=gsub("AN16","1",cellColData$Sample_1)
cellColData$Sample_1=gsub("AN17","2",cellColData$Sample_1)
cellColData$Sample_1=gsub("AN18","3",cellColData$Sample_1)
cellColData$Sample_1=gsub("AN4","1",cellColData$Sample_1)
cellColData$Sample_1=gsub("AN5","2",cellColData$Sample_1)
cellColData$Sample_1=gsub("AN6","3",cellColData$Sample_1)
cellColData$Sample_1=gsub("AN22","1",cellColData$Sample_1)
cellColData$Sample_1=gsub("AN23","2",cellColData$Sample_1)
cellColData$Sample_1=gsub("AN24","3",cellColData$Sample_1)

###Figure2A
ggplot(cellColData, aes(x=Sample_1, y=log10(nFrags),fill=Sample_1)) +
  geom_violin() +
  scale_fill_brewer(palette = "Set3")+
  geom_boxplot(fill="grey", width=.2)+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5, hjust = 0.5))
###Figure2B
ggplot(cellColData, aes(x=Sample_1, y=TSSEnrichment,fill=Sample_1)) +
  geom_violin() +
  scale_fill_brewer(palette = "Set3")+
  geom_boxplot(fill="grey", width=.2)+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5, hjust = 0.5))

###------------------------------Figure3--------------------------------------###
###Figure3A
stallion = c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B", "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","19"="#E6C2DC",
             "10"="#90D5E4", "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416","9"="#D8A767","20"="#3D3D3D")

stallion =c ("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B", "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","9"="#D8A767",
             "10"="#90D5E4", "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416",
             "19"="#E6C2DC","20"="#3D3D3D")
names(stallion)=""
ggplot()+
  geom_point(data=cellColData,aes(x=UMAP_1, y =UMAP_2,color=Celltype),alpha=1, size=0.0001)+
  theme_gray() +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size=25),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))+
  scale_color_manual(values = stallion)+
  guides(colour = guide_legend(override.aes = list(size=5)))

###Figure3B
marker=c("SLC1A2","PDGFRA","MOBP","SLC17A7","GAD1","ADGRE1","P2RY12",
         "TMEM119","TGFB1","CLDN5","FLT1","SIX1","MGP","KLF4","NDRG1")
p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "Celltype", 
  geneSymbol = marker, 
  upstream = 50000,
  downstream = 50000
)
plotPDF(plotList = p, 
        name = "Figure3B.pdf", 
        ArchRProj = proj, 
        addDOC = FALSE, width = 5, height = 5)
###Figure3C
colnames(cellColData)[15]="UMAP_1"
colnames(cellColData)[16]="UMAP_2"
plot=list()
#Auditory_Cortex
sub_coor=subset(cellColData,cellColData$Region=="Auditory_Cortex")
plot[[1]]=ggplot()+
  geom_point(data=cellColData,aes(x=UMAP_1, y =UMAP_2),alpha=1, size=0.0001,colour="gainsboro")+
  theme_gray() +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size=25),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  geom_point(data=sub_coor,aes(x=UMAP_1, y =UMAP_2),alpha=0.3, size=0.0001,colour="red")+
  ggtitle("Auditory_Cortex")+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))
#Motor_Cortex
sub_coor=subset(cellColData,cellColData$Region=="Motor_Cortex")
plot[[2]]=ggplot()+
  geom_point(data=cellColData,aes(x=UMAP_1, y =UMAP_2),alpha=1, size=0.0001,colour="gainsboro")+
  theme_gray() +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size=25),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  geom_point(data=sub_coor,aes(x=UMAP_1, y =UMAP_2),alpha=0.3, size=0.0001,colour="red")+
  ggtitle("Motor_Cortex")+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))
#Primary_Visual_Cortex
sub_coor=subset(cellColData,cellColData$Region=="Primary_Visual_Cortex")
plot[[3]]=ggplot()+
  geom_point(data=cellColData,aes(x=UMAP_1, y =UMAP_2),alpha=1, size=0.0001,colour="gainsboro")+
  theme_gray() +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size=25),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  geom_point(data=sub_coor,aes(x=UMAP_1, y =UMAP_2),alpha=0.3, size=0.0001,colour="red")+
  ggtitle("Primary_Visual_Cortex")+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))
#Somatosensory_cortex
sub_coor=subset(cellColData,cellColData$Region=="Somatosensory_cortex")
plot[[4]]=ggplot()+
  geom_point(data=cellColData,aes(x=UMAP_1, y =UMAP_2),alpha=1, size=0.0001,colour="gainsboro")+
  theme_gray() +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size=25),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  geom_point(data=sub_coor,aes(x=UMAP_1, y =UMAP_2),alpha=0.3, size=0.0001,colour="red")+
  ggtitle("Somatosensory_cortex")+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))

do.call(grid.arrange,c(plot,ncol=2))
###Figure3D
ggplot(cellColData,aes(x=Sample_1,fill=Celltype))+
  geom_bar(stat="count",position="fill")+
  scale_fill_manual(values = stallion)+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))+
  theme(axis.text.x = element_text(angle=45,vjust=0.5))

###Figure3E
cell=sort(unique(proj$Celltype))
for (i in seq(length(cell))) {
  idxPass <- which(proj$Celltype == cell[i])
  cellsPass <- proj$cellNames[idxPass]
  proj_sub=proj[cellsPass, ]
  markersGS <- getMarkerFeatures(
    ArchRProj = proj_sub, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Region",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  
  markerList <- getMarkers(markersGS, cutOff = "FDR <= 1")
  
  write.table(markerList$Auditory_Cortex,paste0(cell[i],"__Auditory_Cortex.txt"),sep = "\t",quote=FALSE)
  write.table(markerList$Motor_Cortex,paste0(cell[i],"__Motor_Cortex.txt"),sep = "\t",quote=FALSE)
  write.table(markerList$Primary_Visual_Cortex,paste0(cell[i],"__Primary_Visual_Cortex"),sep = "\t",quote=FALSE)
  write.table(markerList$Somatosensory_cortex,paste0(cell[i],"__Somatosensory_cortex.txt"),sep = "\t",quote=FALSE)
}

###Statastic & Plot
list=list.files(pattern = "__")
sta=as.data.frame(c(1))
colnames(sta)="Group"
sta$DEG_number=0
deg_list=list()

for (i in seq(64)) {
  input=read.table(list[i],header = T)
  name=gsub(".txt","",list[i])
  input$Group=name
  input_sub=subset(input,input$Log2FC>0.5 & input$FDR<0.05)
  sta[i,1]=name
  sta[i,2]=nrow(input_sub)
  deg_list[[i]]=input_sub
}

sta$Celltype=unlist(lapply(strsplit(as.character(sta$Group),"__"),"[",1))
sta$Region=unlist(lapply(strsplit(as.character(sta$Group),"__"),"[",2))

ggplot(sta,aes(x=Celltype,y=DEG_number,fill=Region))+
  geom_bar(stat="identity")+
  scale_fill_brewer(palette = "Set1")+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))+
  theme(axis.text.x = element_text(angle=45,vjust=0.5))
###------------------------------Figure4--------------------------------------###
enrichMotifs=readRDS("/Users/xiaoyuwei/Desktop/1.Research/8.1Rat_Cortex_Analysis/enrichMotifs_Celltype.RDS")
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 10, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")

