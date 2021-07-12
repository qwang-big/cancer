
load('LungCancer_T_FO13.robj')
pol=list('M1'=m1,'M2'=m2)
s=AddModuleScore(seurat_spatialObj,features=pol)
colnames(s@meta.data)[143]="M1"
colnames(s@meta.data)[144]="M2"
s@meta.data$M1_M2 = s@meta.data$M1-s@meta.data$M2
SpatialFeaturePlot(s,features='M1_M2')
