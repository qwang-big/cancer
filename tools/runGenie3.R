.libPaths(c("/zfssz2/ST_MCHRI/COHORT/dingxiangning/conda/envs/GRNs/lib/R/library",
            "/zfssz2/ST_MCHRI/COHORT/dingxiangning/conda/envs/test_env/lib/R/library",
            "/zfssz2/ST_MCHRI/COHORT/USER/zhujiacheng/Tools/miniconda3/lib/R/library"))
library(Seurat)
library(tidyverse)
library(patchwork)
library(GENIE3)
library(RcisTarget)
library(AUCell)
library(SCENIC)

scRNA <- readRDS("/hwfssz5/ST_PRECISION/TOMCAT/xuyoujia/Macaca_eye/iDrop_scRNA_Macaque_IN_EX/rds/M1_kid/M1_kid.rds")

##准备细胞meta信息
cellInfo <- data.frame(scRNA@meta.data)
colnames(cellInfo)[which(colnames(cellInfo)=="orig.ident")] <- "sample"
colnames(cellInfo)[which(colnames(cellInfo)=="seurat_clusters")] <- "cluster"
colnames(cellInfo)[which(colnames(cellInfo)=="celltype")] <- "celltype"
cellInfo <- cellInfo[,c("sample","cluster","celltype")]
#saveRDS(cellInfo, file="int/cellInfo.Rds")

##准备表达矩阵
#为了节省计算资源，随机抽取1000个细胞的数据子集
subcell <- sample(colnames(scRNA),1000)
scRNAsub <- scRNA[,subcell]
#saveRDS(scRNAsub, "scRNAsub.rds")
exprMat <- as.matrix(scRNAsub@assays$RNA@counts)

##设置分析环境
mydbDIR <- "cisTarget"
mydbs <- c("hg19-500bp-upstream-7species.mc9nr.feather",
           "hg19-tss-centered-10kb-7species.mc9nr.feather")
names(mydbs) <- c("500bp", "10kb")
scenicOptions <- initializeScenic(org="hgnc", 
                                  nCores=8,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "HNSCC")
#scenicOptions <- readRDS("int/scenicOptions.rds")

##==转录调控网络推断==##
##基因过滤
#过滤标准是基因表达量之和>细胞数*3%，且在1%的细胞中表达
genesKept <- geneFiltering(exprMat, scenicOptions, 
                           minCountsPerGene = 3 * 0.01 * ncol(exprMat), 
                           minSamples = ncol(exprMat) * 0.01)
exprMat_filtered <- exprMat[genesKept, ]

##计算相关性矩阵
runCorrelation(exprMat_filtered, scenicOptions)

##TF-Targets相关性回归分析
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions, nParts = 20)
#这一步消耗的计算资源非常大，个人电脑需要几个小时的运行时间

##推断共表达模块
runSCENIC_1_coexNetwork2modules(scenicOptions)

##推断转录调控网络（regulon）
runSCENIC_2_createRegulons(scenicOptions)
#以上代码可增加参数coexMethod=c("w001", "w005", "top50", "top5perTarget", "top10perTarget", "top50perTarget"))
#默认6种方法的共表达网络都计算，可以少选几种方法以减少计算量

scRNA <- readRDS("/hwfssz5/ST_PRECISION/TOMCAT/xuyoujia/Macaca_eye/iDrop_scRNA_Macaque_IN_EX/rds/M1_kid/M1_kid.rds")
##==regulon活性评分与可视化==##
##regulons计算AUC值并进行下游分析
exprMat_all <- as.matrix(scRNA@assays$RNA@counts)
exprMat_all <- log2(exprMat_all+1)
runSCENIC_3_scoreCells(scenicOptions, exprMat=exprMat_all)

#使用shiny互动调整阈值
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_all)
savedSelections <- shiny::runApp(aucellApp)

#保存调整后的阈值
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

#二进制转换及衍生分析 
runSCENIC_4_aucell_binarize(scenicOptions, exprMat=exprMat_all)

