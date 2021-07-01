library(Seurat)
#.libPaths('/hwfssz1/ST_CANCER/SCOL/SHARE/APP/R/Rlib_zhangcj/R-3.4.1/')
library(infercnvPlus)
library(ComplexHeatmap)
args <- commandArgs(trailingOnly = TRUE)
# Run example with built-in data
load('genomic_pos.rda')
load(args[1])
print(dim(seurat_spatialObj))
Tumor_seurat = CreateSeuratObject(seurat_spatialObj@assays$Spatial@counts,min.cells = 0,min.features = 0)
print(dim(Tumor_seurat))
load(args[2])
print(dim(seurat_spatialObj))
Normal_seurat = CreateSeuratObject(as.matrix(seurat_spatialObj@assays$Spatial@counts)[,sample(colnames(seurat_spatialObj),3000)],min.cells = 0,min.features = 0)
print(dim(Normal_seurat))
all = merge(Tumor_seurat,Normal_seurat,add.cell.ids = c('Tumor','Normal'))
print(dim(all))
all$sample=substr(rownames(all@meta.data),start = 1,stop = 4)
ref_obs=rownames(all@meta.data[all@meta.data$sample=='Norm',])
print(str(ref_obs))
rm(seurat_spatialObj)
rm(Tumor_seurat)
rm(Normal_seurat)
gc()
expr=as.matrix(all@assays$RNA@counts)
nr <- nrow(expr)
n <- ceiling(nr/20)
i=rep(1:20, each=ceiling(nr/n), length.out=nr)
rm(all)
gc()
counts_to_tpm <- function(x) t(t(x)*1e6/colSums(x))
expr <- counts_to_tpm(expr[i==args[3],])
expr_tr <- log2(expr + 1)
expr_tr = round(expr_tr*10)
rm(expr)
gc()
#load('FJ13_CNV.robj')
# Data tranforming: genes(rows) X cells(columns)
## For 10X counts data 
#expr_tr <- umi_to_log2tpm(expr)
## For Smart-seq2 TPM values
#expr_tr <- log2(expr + 1)
# Calucate cnv score
cnv_obj <- inferCNV(data = expr_tr,
                    gene_pos = genomic_pos,
                    cutoff = 1, #0.1, # use 1 for smart-seq, 0.1 for 10x-genomics
                    reference_obs = ref_obs,
                    window_size = 101,
                    out_path = "output_dir", # dir is auto-created for storing outputs
                    noise_filter = NULL,
                    vis_bounds = "-0.5,0.5")
rm(expr_tr)
gc()
#save(cnv_obj,file = 'FJ13_CNV_result_new_0.5.Robj')
# Cluster cells and visualize
cnv_obj <- visualCNV(data = cnv_obj,
                     cutree_k = 2,
                     out_file = paste0(args[1],".plot_cnv_0.5.png"))
# Extract cells from the specific subtrees
cnv_obj <- extractCells(data = cnv_obj,
                        subtrees = 2,
                        lab_to_rm = "ref")
# Get cell barcode
#cells <- cnv_obj$target
saveRDS(cnv_obj,file = paste0(args[1],".",args[2],".",args[3],".rds"))
