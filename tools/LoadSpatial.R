########
#
#   LoadSpatial.R
#
#   Author: YU Hao     (yuhao@genomics.cn)
#   Author: BAI yinqi  (baiyinqi@genomics.cn)
#
#   Date:   2020-12-27  (v0.1)
#   Date:   2021-02-22  (v0.2) To generate a white background & to readjust the direction of Y-axis
#   Date:   2021-02-23  (v0.3) To resume a black background and to readjust the background XY-axes
#   Date:   2021-02-27  (v0.4) 1) To fix up unfounding of some row, due to their special row name,
#.                                such as BIN.100000, BIN.200000 ...
#                                 (Thank ZOU Xuanxuan and HOU Liangzhen for reporting this indiscoverable bug)
#                              2) To add a new optional option on args[7] for assigning Bin size
#   Date:   2021-02-28  (v0.5) To readjust the calculation method of Y-axis for eliminating influence of boundary
#   Date:   2021-03-04  (v0.6) To resume Y-axis direction
#   Date:   2021-03-09  (v0.7) To purge the Bins carrying 0 umi for making following 'SCTransform' successfully
#   Date:   2021-03-10  (v0.8) To add a function for generating AnnData (h5ad) file for the following analyzing
#                              by ScanPy or SquidPy
#                              (Use '--h5ad' at the last of argument list to active this option)
#
########




library(data.table)
library(Matrix)
library(rjson)
library(Seurat)


##
args <- commandArgs()

bs <- 50

if (length(args) > 6) {if (args[7] != '--h5ad') {bs <- as.numeric(args[7])}}


##
io.dir <- args[6]

dat.file <- file.path(io.dir, 'merge_GetExp_gene.txt')


##
dat <- fread(file = dat.file)

dat$x <- trunc((dat$x - min(dat$x)) / bs + 1)
dat$y <- trunc((dat$y - min(dat$y)) / bs + 1)


##
if ('MIDCounts' %in% colnames(dat)) {
    dat <- dat[, sum(MIDCounts), by = .(geneID, x, y)]
} else {
    dat <- dat[, sum(UMICount), by = .(geneID, x, y)]
}

dat$bin_ID <- max(dat$x) * (dat$y - 1) + dat$x

bin.coor <- dat[, sum(V1), by = .(x, y)]


##
geneID <- seq(length(unique(dat$geneID)))
hash.G <- data.frame(row.names = unique(dat$geneID), values = geneID)
gen <- hash.G[dat$geneID, 'values']


##
bin_ID <- unique(dat$bin_ID)
hash.B <- data.frame(row.names = sprintf('%d', bin_ID), values = bin_ID)
bin <- hash.B[sprintf('%d', dat$bin_ID), 'values']


##
cnt <- dat$V1


##
rm(dat)
gc()


##
tissue_lowres_image <- matrix(0, max(bin.coor$y), max(bin.coor$x))

tissue_positions_list <- data.frame(row.names = paste('BIN', rownames(hash.B), sep = '.'),
                                    tissue = 1, 
                                    row = bin.coor$y, col = bin.coor$x,
                                    imagerow = bin.coor$y, imagecol = bin.coor$x)

scalefactors_json <- toJSON(list(fiducial_diameter_fullres = 1,
                                 tissue_hires_scalef = 1,
                                 tissue_lowres_scalef = 1))


##
mat <- sparseMatrix(i = gen, j = bin, x = cnt)

rownames(mat) <- rownames(hash.G)
colnames(mat) <- paste('BIN', sprintf('%d', seq(max(hash.B[, 'values']))), sep = '.')

seurat_spatialObj <- CreateSeuratObject(mat, project = 'Spatial', assay = 'Spatial')


##
generate_spatialObj <- function(image, scale.factors, tissue.positions, filter.matrix = TRUE) 
{
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


##
spatialObj <- spatialObj[Cells(seurat_spatialObj)]
DefaultAssay(spatialObj) <- 'Spatial'

seurat_spatialObj[['slice1']] <- spatialObj


##
seurat_spatialObj <- subset(seurat_spatialObj, subset = nCount_Spatial > 0)

##
if (args[length(args)] == '--h5ad') {
    write_h5ad.fix <- function(anndata, filename, compression = NULL, compression_opts = NULL, as_dense = list()) 
    {
        anndata$write_h5ad(filename = filename, 
                           compression = compression, 
                           compression_opts = compression_opts, 
                           as_dense = as_dense)
    }
    
    
    genes <- as.data.frame(rownames(seurat_spatialObj), row.names = rownames(seurat_spatialObj))
    names(genes) <- 'Gene'
    
    cells <- as.data.frame(colnames(seurat_spatialObj), row.names = colnames(seurat_spatialObj))
    names(cells) <- 'CellID'
    
    row <- seurat_spatialObj@images$slice1@coordinates$row
    col <- seurat_spatialObj@images$slice1@coordinates$col
    coordinates <- list(matrix(c(row, col), ncol=2)); names(coordinates) <- 'spatial'
    
    ad <- anndata::AnnData(X = seurat_spatialObj@assays$Spatial@counts, obs = genes, var = cells, varm = coordinates)
    ad <- ad$T
    
    DUMP <- write_h5ad.fix(ad, 'squidpy_spatialObj.h5ad')
} else {
    saveRDS(seurat_spatialObj, file = file.path(io.dir, 'seurat_spatialObj.rds'))
}
