x=readRDS("Annotations/Motif-Matches-In-Peaks.rds")
seqlevels(x)<- sub('chr','MFA',seqlevels(x))
motifs <- TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, list(species=NULL, collection="CORE"))
genome = BSgenome.Mfascicularis.NCBI.5.0::BSgenome.Mfascicularis.NCBI.5.0
motif_ix <- matchMotifs(motifs, x, genome)
