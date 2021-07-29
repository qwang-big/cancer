library(Seurat)
library(ggplot2)
library(ArchR)
setwd("snATAC-seq/snATAC_cortex_TSS5_Frag1000/")
atac <- readRDS("Save-ArchR-Project.rds")
ff=dir("ArrowFiles")
atac@sampleColData <- DataFrame(row.names = substr(ff,1,nchar(ff)-6), ArrowFiles = paste0('ArrowFiles/',ff))
seqlevels(atac@peakSet)<- sub('chr','MFA',seqlevels(atac@peakSet))
atac=addMotifAnnotations(ArchRProj = atac, motifSet = "JASPAR2020", name = "Motif")
saveRDS(atac, "Save-ArchR-Project.rds")
#change seqnames back
seqlevels(atac@peakSet)<- sub('MFA','chr',seqlevels(atac@peakSet))
x=readRDS("Annotations/Motif-Positions-In-Peaks.rds")
seqlevels(x)<- sub('MFA','chr',seqlevels(x))
saveRDS(x, file="Annotations/Motif-Positions-In-Peaks.rds")
x=readRDS("Annotations/Motif-Matches-In-Peaks.rds")
seqlevels(x)<- sub('MFA','chr',seqlevels(x))
saveRDS(x, file="Annotations/Motif-Matches-In-Peaks.rds")
#compute DeviationsMatrix
atac <- addBgdPeaks(atac)
atac <- addDeviationsMatrix(atac, peakAnnotation = "Motif", force = TRUE)
atac@cellColData$age_group=as.character(atac@cellColData$age_group)
motifs <- c("SNAI1", "ZEB1", "TCF4", "ASCL1", "MYOD1", "Nr2f6")
x=readRDS('Annotations/Motif-Positions-In-Peaks.rds')
x=x[names(x)[unlist(lapply(motifs,function(d) grep(d,names(x))))]]
## export to network
p2g <- metadata(atac@peakSet)$Peak2GeneLinks
p=metadata(p2g)$peakSet[p2g$idxATAC]
g=metadata(p2g)$geneSet[p2g$idxRNA]
df=do.call(rbind,lapply(seq_along(x),function(i) data.frame(names(x)[i],g[queryHits(findOverlaps(p,x[[i]]))]$name)))
write.table(df,file='~/b/df.sif',sep='\t-\t',row.names=F,col.names=F,quote=F)
## correlate to GA 
gm <- getMatrixFromProject(atac, useMatrix='GeneScoreMatrix')
gsm <- assays(gm)$GeneScoreMatrix
rownames(gsm) <- rowData(gm)$name
d=table(g)
r=rowMeans(gsm[names(d),])
df=data.frame(x=r,y=as.integer(d))
df=df[df$x>0,]
m=lm(y~x,df)
p=ggplot(df,aes(x,y))+geom_point(size=3,shape=20)+scale_x_continuous(limits=c(0,10))+labs(x='Gene activity',y='Number of enhancers')+geom_abline(slope = coef(m)[[1]], intercept = coef(m)[[2]],color='red')+ annotate(geom="text", label=paste0('r=',round(cor(df$x,df$y),2)), x=8, y=500, size=8, color='red')

