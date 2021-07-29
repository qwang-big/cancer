af.availableSeqnames = utils::getFromNamespace(".availableSeqnames", "ArchR")
af.getRowSums = utils::getFromNamespace(".getRowSums", "ArchR")
ArrowFiles <- getArrowFiles(atac)
useMatrix <- "PeakMatrix"
availableChr <- af.availableSeqnames(ArrowFiles, useMatrix)
rS <- af.getRowSums(ArrowFiles = ArrowFiles, seqnames = availableChr, useMatrix = useMatrix, filter0 = FALSE)
rS$start <- start(atac@peakSet)
rS$end <- end(atac@peakSet)
rS$GC <- atac@peakSet$GC
se <- SummarizedExperiment::SummarizedExperiment(assays = SimpleList(counts = as.matrix(data.frame(rS$rowSums, 1))), rowData = DataFrame(bias = rS$GC, start = rS$start, end = rS$end))
bgdPeaks <- chromVAR::getBackgroundPeaks(object = se, bias = rowData(se)$bias, iterations = 50, w = 0.1, bs = 50)
motifs <- getJasparMotifs()
motif_ix <- matchMotifs(motifs, counts_filtered, genome = BSgenome.Mfascicularis.NCBI.5.0)
