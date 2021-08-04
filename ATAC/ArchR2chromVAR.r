markerTest <- getMarkerFeatures(
  ArchRProj = atac, 
  useMatrix = "PeakMatrix",
  groupBy = "age_group",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "4_Geriatric",
  bgdGroups = "3_Old"
)
pma <- markerPlot(seMarker = markerTest, name = "4_Geriatric", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
pma
atac@peakAnnotation$Motif$Positions="Annotations/Motif-Positions-In-Peaks.rds"
atac@peakAnnotation$Motif$Matches="Annotations/Motif-Matches-In-Peaks.rds"
seqlevels(projHeme5@peakSet)<- sub('MFA','chr',seqlevels(projHeme5@peakSet))
motifsUp <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = atac,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
df=cbind(d1,d2[rownames(d1),],d3[rownames(d1),])
dd=df[apply(df[,c(3,6,9)],1,function(d) any(d<=10)),]
rownames(dd)=str_match(dd[,1],"^[A-Za-z0-9]+")
dd=as.matrix(dd[,c(2,5,8)])
colnames(dd)=c('Young','Middle','Old')
heatmap.2(dd,Colv=F,density.info='none',col=colorpanel(10,'black','red'), srtCol=45)

x=readRDS("Annotations/Motif-Matches-In-Peaks.rds")
seqlevels(x)<- sub('chr','MFA',seqlevels(x))
motifs <- TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, list(species=NULL, collection="CORE"))
genome = BSgenome.Mfascicularis.NCBI.5.0::BSgenome.Mfascicularis.NCBI.5.0
dd=unlist(lapply(motifs,TFBSTools::name))
motif_ix <- matchMotifs(motifs, x, genome)
