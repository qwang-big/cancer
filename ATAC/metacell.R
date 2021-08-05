library("metacell")
library("Seurat")
if(!dir.exists("testdb")) dir.create("testdb/")
scdb_init("testdb/", force_reinit=T)
seu = readRDS("xxx.seurat.rds")
sce = as.SingleCellExperiment(seu, assay = "Spatial")
mat = scm_import_sce_to_mat(sce)
mat = scdb_add_mat("test", mat)
mat = scdb_mat("test")
nms = c(rownames(mat@mat), rownames(mat@ignore_gmat))
bad_genes = grep("^(MT|mt)-", nms, ignore.case = T)
mcell_mat_ignore_genes(new_mat_id="test", mat_id="test", bad_genes, reverse=F)
mcell_add_gene_stat(gstat_id="test", mat_id="test", force=T)
mcell_gset_filter_varmean(gset_id="test_feats", gstat_id="test", T_vm=0.08, force_new=T)
mcell_gset_filter_cov(gset_id = "test_feats", gstat_id="test", T_tot=100, T_top3=2)
mcell_add_cgraph_from_mat_bknn(mat_id="test",
                               gset_id = "test_feats",
                               graph_id="test_graph",
                               K=100,
                               dsamp=T)
mcell_coclust_from_graph_resamp(
  coc_id="test_coc500",
  graph_id="test_graph",
  min_mc_size=20,
  p_resamp=0.75, n_resamp=500)
mcell_mc_from_coclust_balanced(
  coc_id="test_coc500",
  mat_id= "test",
  mc_id= "test_mc",
  K=30, min_mc_size=30, alpha=2)
mcell_mc_split_filt(new_mc_id="test_mc_f",
                    mc_id="test_mc",
                    mat_id="test",
                    T_lfc=3, plot_mats=F)
mcell_gset_from_mc_markers(gset_id="test_markers", mc_id="test_mc_f")
mc_f<- scdb_mc("test_mc_f")
mc_f@colors <- colorRampPalette(c("darkgray", "burlywood1", "chocolate4","orange", "red", "purple", "blue","darkgoldenrod3", "cyan"))(ncol(mc_f@mc_fp))
scdb_add_mc("test_mc_f",mc_f)
mc_f <- scdb_mc("test_mc_f")
mcell_gset_from_mc_markers(gset_id="test_markers", mc_id="test_mc_f")
mcell_mc_plot_marks(mc_id="test_mc_f", gset_id="test_markers", mat_id="test",plot_cells = F)
