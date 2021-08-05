library('psupertime')
library(ArchR)
atac=readRDS('Save-ArchR-Project.rds')
m=getMatrixFromProject(atac, useMatrix='GeneScoreMatrix')
sce=as(m, "SingleCellExperiment")
rownames(sce)=rowData(sce)$name
assayNames(sce)='logcounts'
psuper = psupertime(sce, sce$age_group, sel_genes="tf_human")
i=order(psuper$proj_dt$psuper)
m=as.matrix(assay(sce)[which(rownames(sce) %in% tf_human),])
df=data.frame(lapply(split(i, rep(seq_len(ceiling(length(i)/50)), each = 50)[seq_len(length(i))]),function(j)rowMeans(m[,j])))

psuper_plot1 <- plot_labels_over_psupertime(psuper, label_name='group')
psuper_plot2 <- plot_identified_genes_over_psupertime(psuper, label_name='group' , n_to_plot = 50,plot_ratio = 1)
ggsave("~/b/pic/label.png", psuper_plot1)
ggsave("~/b/pic/genes.png", psuper_plot2)

s=readRDS('seurat.rds')
sce=as.SingleCellExperiment(s)
normalize <- function(x) (x-mean(x))/sd(x)
dd=as.matrix(t(apply(df,1,normalize)))
dd=dd[!apply(dd,1,function(d) all(is.na(d))),]
png('~/b/pic/htm.png')
heatmap.2(dd,Colv=F,trace='none',col=greenred(100),breaks=seq(-4,4,length.out=101))
dev.off()

atac@cellColData$pseudotime=0
atac@cellColData$pseudotime[i]=rep(seq_len(ceiling(length(i)/50)), each = 50)[seq_len(length(i))]
atac@cellColData$pseudotime=as.factor(atac@cellColData$pseudotime)
seGroupMotif <- getGroupSE(atac, useMatrix = "MotifMatrix", groupBy = "pseudotime")
