library(data.table)
library(Matrix)
library(rjson)
library(Seurat)

ff=dir()
lapply(ff[grep('.gem$',ff)],function(f){
  bs <- 50
  dat <- fread(file = f)
  dat$x <- trunc((dat$x - min(dat$x)) / bs + 1)
  dat$y <- trunc((dat$y - min(dat$y)) / bs + 1)
  dat <- dat[, sum(MIDCounts), by = .(geneID, x, y)]
  dat$bin_ID <- max(dat$x) * (dat$y - 1) + dat$x
  bin.coor <- dat[, sum(V1), by = .(x, y)]
  geneID <- seq(length(unique(dat$geneID)))
  hash.G <- data.frame(row.names = unique(dat$geneID), values = geneID)
  gen <- hash.G[dat$geneID, 'values']
  bin_ID <- unique(dat$bin_ID)
  hash.B <- data.frame(row.names = as.character(bin_ID), values = bin_ID)
  bin <- hash.B[as.character(dat$bin_ID), 'values']
  cnt <- dat$V1
  rm(dat)
  gc()
  tissue_lowres_image <- matrix(1, max(bin.coor$x), max(bin.coor$y))
  tissue_positions_list <- data.frame(row.names = paste('BIN', rownames(hash.B), sep = '.'),
                                      tissue = 1, 
                                      row = bin.coor$y, col = bin.coor$x,
                                      imagerow = -bin.coor$y, imagecol = bin.coor$x)
  scalefactors_json <- toJSON(list(fiducial_diameter_fullres = 1,
                                   tissue_hires_scalef = 1,
                                   tissue_lowres_scalef = 1))
  mat <- sparseMatrix(i = gen, j = bin, x = cnt)
  rownames(mat) <- rownames(hash.G)
  colnames(mat) <- paste('BIN', seq(max(hash.B[, 'values'])), sep = '.')
  seurat_spatialObj <- CreateSeuratObject(mat, project = 'Spatial', assay = 'Spatial')
  generate_spatialObj <- function (image, scale.factors, tissue.positions, filter.matrix = TRUE) {
    if (filter.matrix) {
      tissue.positions <- tissue.positions[which(tissue.positions$tissue == 1), , drop = FALSE]
    }
    
    unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
    
    spot.radius <- unnormalized.radius / max(dim(image))
    
    return(new(Class = 'VisiumV1', 
               image = image, 
               scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef, 
                                            fiducial = scale.factors$fiducial_diameter_fullres, 
                                            hires = scale.factors$tissue_hires_scalef, 
                                            lowres = scale.factors$tissue_lowres_scalef), 
               coordinates = tissue.positions, 
               spot.radius = spot.radius))
  }
  spatialObj <- generate_spatialObj(image = tissue_lowres_image, 
                                    scale.factors = fromJSON(scalefactors_json), 
                                    tissue.positions = tissue_positions_list)
  spatialObj <- spatialObj[Cells(seurat_spatialObj)]
  DefaultAssay(spatialObj) <- 'Spatial'
  seurat_spatialObj[['slice1']] <- spatialObj
  seurat_spatialObj <- subset(seurat_spatialObj, subset = nFeature_Spatial > 0)
  seurat_spatialObj <- SCTransform(seurat_spatialObj, verbose = FALSE,assay = 'Spatial')
  seurat_spatialObj <- RunPCA(seurat_spatialObj, verbose = FALSE, assay = 'SCT')
  seurat_spatialObj <- RunUMAP(seurat_spatialObj, dims = 1:15, verbose = FALSE)
  seurat_spatialObj <- FindNeighbors(seurat_spatialObj, dims = 1:15, verbose = FALSE)
  seurat_spatialObj <- FindClusters(seurat_spatialObj, resolution = 0.3, verbose = FALSE)
  marker_NPC=FindAllMarkers(seurat_spatialObj,only.pos = T)
  save(seurat_spatialObj,marker_NPC,file=paste0('~/c/CNV/',f))
})
