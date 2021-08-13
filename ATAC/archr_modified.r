addGroupCoverages1=function (ArchRProj = NULL, groupBy = "Clusters", useLabels = TRUE, 
    minCells = 40, maxCells = 500, maxFragments = 25 * 10^6, 
    minReplicates = 2, maxReplicates = 5, sampleRatio = 0.8, 
    kmerLength = 6, threads = getArchRThreads(), returnGroups = FALSE, 
    parallelParam = NULL, force = FALSE, verbose = TRUE, logFile = createLogFile("addGroupCoverages")) 
{
    ArchR:::.validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
    ArchR:::.validInput(input = groupBy, name = "groupBy", valid = c("character"))
    ArchR:::.validInput(input = useLabels, name = "useLabels", valid = c("boolean"))
    ArchR:::.validInput(input = minCells, name = "minCells", valid = c("integer"))
    ArchR:::.validInput(input = maxCells, name = "maxCells", valid = c("integer"))
    ArchR:::.validInput(input = maxFragments, name = "maxFragments", 
        valid = c("integer"))
    ArchR:::.validInput(input = minReplicates, name = "minReplicates", 
        valid = c("integer"))
    ArchR:::.validInput(input = maxReplicates, name = "maxReplicates", 
        valid = c("integer"))
    ArchR:::.validInput(input = sampleRatio, name = "sampleRatio", valid = c("numeric"))
    ArchR:::.validInput(input = kmerLength, name = "kmerLength", valid = c("integer"))
    ArchR:::.validInput(input = threads, name = "threads", valid = c("integer"))
    ArchR:::.validInput(input = returnGroups, name = "returnGroups", 
        valid = c("boolean"))
    ArchR:::.validInput(input = parallelParam, name = "parallelParam", 
        valid = c("parallelparam", "null"))
    ArchR:::.validInput(input = force, name = "force", valid = c("boolean"))
    ArchR:::.validInput(input = verbose, name = "verbose", valid = c("boolean"))
    ArchR:::.validInput(input = logFile, name = "logFile", valid = c("character"))
    if (minReplicates < 2) {
        stop("minReplicates must be at least 2!")
    }
    tstart <- Sys.time()
    ArchR:::.startLogging(logFile = logFile)
    ArchR:::.logThis(mget(names(formals()), sys.frame(sys.nframe())), 
        "addGroupCoverages Input-Parameters", logFile = logFile)
    Params <- SimpleList(groupBy = groupBy, minCells = minCells, 
        maxCells = maxCells, minReplicates = minReplicates, sampleRatio = sampleRatio, 
        kmerLength = kmerLength)
    if (is.null(ArchRProj@projectMetadata$GroupCoverages)) {
        ArchRProj@projectMetadata$GroupCoverages <- SimpleList()
    }
    if (!returnGroups) {
        if (!is.null(ArchRProj@projectMetadata$GroupCoverages[[groupBy]])) {
            if (!force) {
                stop("Group Coverages Already Computed, Set force = TRUE to continue!")
            }
        }
    }
    else {
        if (!is.null(ArchRProj@projectMetadata$GroupCoverages[[groupBy]])) {
            if (!force) {
                message("Group Coverages Already Computed Returning Groups, Set force = TRUE to Recompute!")
                return(ArchRProj@projectMetadata$GroupCoverages[[groupBy]])
            }
        }
    }
    cellNames <- rownames(getCellColData(ArchRProj))
    groups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
    if (any(is.na(groups))) {
        cellNames <- cellNames[!is.na(groups)]
        groups <- groups[!is.na(groups)]
    }
    uniqueGroups <- gtools::mixedsort(unique(groups))
    tableGroups <- table(groups)[uniqueGroups]
    cellGroups <- lapply(seq_along(uniqueGroups), function(x) {
        subColDat <- getCellColData(ArchRProj)[which(groups == 
            uniqueGroups[x]), ]
        cellNamesx <- rownames(subColDat)
        if (useLabels) {
            sampleLabelsx <- paste0(subColDat$Sample)
        }
        else {
            sampleLabelsx <- NULL
        }
        outListx <- ArchR:::.identifyGroupsForPseudoBulk(cells = cellNamesx, 
            sampleLabels = sampleLabelsx, useLabels = useLabels, 
            minCells = minCells, maxCells = maxCells, minReplicates = minReplicates, 
            maxReplicates = maxReplicates, sampleRatio = sampleRatio, 
            prefix = sprintf("%s (%s of %s) :", uniqueGroups[x], 
                x, length(uniqueGroups)), logFile = logFile)
        if (is.null(outListx)) {
            return(NULL)
        }
        if (is.null(names(outListx))) {
            names(outListx) <- paste0("Rep", seq_along(outListx))
        }
        else if (any(names(outListx) == "")) {
            names(outListx)[which(names(outListx) == "")] <- paste0("Rep", 
                which(names(outListx) == ""))
        }
        outListx
    }) %>% SimpleList
    names(cellGroups) <- uniqueGroups
    Params$cellGroups <- cellGroups
    it <- 0
    for (i in seq_along(cellGroups)) {
        for (j in seq_along(cellGroups[[i]])) {
            if (sum(getCellColData(ArchRProj, "nFrags")[cellGroups[[i]][[j]], 
                ]) > maxFragments) {
                it <- it + 1
                nFrags <- getCellColData(ArchRProj, "nFrags")[cellGroups[[i]][[j]], 
                  ]
                cells <- cellGroups[[i]][[j]][order(nFrags)]
                nFrags <- nFrags[order(nFrags)]
                cellGroups[[i]][[j]] <- cells[which(cumsum(nFrags) < 
                  maxFragments)]
            }
        }
    }
    if (it > 0) {
        ArchR:::.logDiffTime(sprintf("Further Sampled %s Groups above the Max Fragments!", 
            it), tstart)
    }
    if (returnGroups) {
        return(cellGroups)
    }
    dir.create(file.path(getOutputDirectory(ArchRProj), "GroupCoverages"), 
        showWarnings = FALSE)
    dir.create(file.path(getOutputDirectory(ArchRProj), "GroupCoverages", 
        groupBy), showWarnings = FALSE)
    unlistGroups <- lapply(seq_along(cellGroups), function(x) {
        if (is.null(cellGroups[[x]])) {
            NULL
        }
        else {
            names(cellGroups[[x]]) <- paste0(names(cellGroups)[x], 
                "._.", names(cellGroups[[x]]))
            cellGroups[[x]]
        }
    }) %>% SimpleList %>% unlist()
    args <- list()
    args$X <- seq_along(unlistGroups)
    args$FUN <- ArchR:::.createCoverages
    args$cellGroups <- unlistGroups
    args$genome <- getGenome(ArchRProj)
    args$kmerLength <- kmerLength
    args$ArrowFiles <- getArrowFiles(ArchRProj)
    args$availableChr <- ArchR:::.availableSeqnames(getArrowFiles(ArchRProj))
    args$chromLengths <- getChromLengths(ArchRProj)
    args$cellsInArrow <- cellsInArrow <- split(rownames(getCellColData(ArchRProj)), 
        stringr::str_split(rownames(getCellColData(ArchRProj)), 
            pattern = "\\#", simplify = TRUE)[, 1])
    args$covDir <- file.path(getOutputDirectory(ArchRProj), "GroupCoverages", 
        groupBy)
    args$parallelParam <- parallelParam
    args$threads <- threads
    args$verbose <- verbose
    args$tstart <- tstart
    args$logFile <- logFile
    args$registryDir <- file.path(getOutputDirectory(ArchRProj), 
        "GroupCoverages", "batchRegistry")
    h5disableFileLocking()
    ArchR:::.logDiffTime(sprintf("Creating Coverage Files!"), tstart, 
        addHeader = FALSE)
    batchOut <- ArchR:::.batchlapply(args)
    coverageFiles <- lapply(seq_along(batchOut), function(x) batchOut[[x]]$covFile) %>% 
        unlist
    nCells <- lapply(seq_along(batchOut), function(x) batchOut[[x]]$nCells) %>% 
        unlist
    nFragments <- lapply(seq_along(batchOut), function(x) batchOut[[x]]$nFragments) %>% 
        unlist
    coverageMetadata <- DataFrame(Group = stringr::str_split(names(unlistGroups), 
        pattern = "\\._.", simplify = TRUE)[, 1], Name = names(unlistGroups), 
        File = coverageFiles, nCells = nCells, nInsertions = nFragments * 
            2)
    ArchR:::.logDiffTime(sprintf("Adding Kmer Bias to Coverage Files!"), 
        tstart, addHeader = FALSE)
    o <- addKmerBiasToCoverage(coverageMetadata = coverageMetadata, 
        genome = getGenome(ArchRProj), kmerLength = kmerLength, 
        threads = threads, verbose = FALSE, logFile = logFile)
    ArchRProj@projectMetadata$GroupCoverages[[groupBy]] <- SimpleList(Params = Params, 
        coverageMetadata = coverageMetadata)
    h5enableFileLocking()
    ArchR:::.logDiffTime(sprintf("Finished Creation of Coverage Files!"), 
        tstart, addHeader = FALSE)
    .endLogging(logFile = logFile)
    ArchRProj
}
addKmerBiasToCoverage=function (coverageMetadata = NULL, genome, kmerLength = NULL, 
    threads = NULL, verbose = TRUE, tstart = NULL, logFile = NULL) 
{
    ArchR:::.logThis(append(args, mget(names(formals()), sys.frame(sys.nframe()))), 
        "kmerBias-Parameters", logFile = logFile)
    ArchR:::.requirePackage("Biostrings", source = "bioc")
    #.requirePackage(genome)
    ArchR:::.requirePackage("Biostrings", source = "bioc")
    BSgenome <- eval(parse(text = genome))
    BSgenome <- validBSgenome(BSgenome)
    if (is.null(tstart)) {
        tstart <- Sys.time()
    }
    coverageFiles <- coverageMetadata$File
    names(coverageFiles) <- coverageMetadata$Name
    availableChr <- ArchR:::.availableSeqnames(coverageFiles, "Coverage")
    biasList <- .safelapply(seq_along(availableChr), function(x) {
        ArchR:::.logMessage(sprintf("Kmer Bias %s (%s of %s)", availableChr[x], 
            x, length(availableChr)), logFile = logFile)
        message(availableChr[x], " ", appendLF = FALSE)
        chrBS <- BSgenome[[availableChr[x]]]
        exp <- Biostrings::oligonucleotideFrequency(chrBS, width = kmerLength)
        obsList <- lapply(seq_along(coverageFiles), function(y) {
            ArchR:::.logMessage(sprintf("Coverage File %s (%s of %s)", 
                availableChr[x], y, length(coverageFiles)), logFile = logFile)
            tryCatch({
                obsx <- .getCoverageInsertionSites(coverageFiles[y], 
                  availableChr[x]) %>% {
                  BSgenome::Views(chrBS, IRanges(start = . - 
                    floor(kmerLength/2), width = kmerLength))
                } %>% {
                  Biostrings::oligonucleotideFrequency(., width = kmerLength, 
                    simplify.as = "collapsed")
                }
                tryCatch({
                  gc()
                }, error = function(e) {
                })
                obsx
            }, error = function(e) {
                errorList <- list(y = y, coverageFile = coverageFiles[y], 
                  chr = availableChr[x], iS = tryCatch({
                    ArchR:::.getCoverageInsertionSites(coverageFiles[y], 
                      availableChr[x])
                  }, error = function(e) {
                    "Error .getCoverageInsertionSites"
                  }))
                ArchR:::.logError(e, fn = ".addKmerBiasToCoverage", info = "", 
                  errorList = errorList, logFile = logFile)
            })
        }) %>% SimpleList
        names(obsList) <- names(coverageFiles)
        SimpleList(expected = exp, observed = obsList)
    }, threads = threads) %>% SimpleList
    names(biasList) <- availableChr
    ArchR:::.logMessage("Completed Kmer Bias Calculation", logFile = logFile)
    for (i in seq_along(biasList)) {
        if (i == 1) {
            expAll <- biasList[[i]]$expected
            obsAll <- biasList[[i]]$observed
        }
        else {
            expAll <- expAll + biasList[[i]]$expected
            for (j in seq_along(obsAll)) {
                obsAll[[j]] <- obsAll[[j]] + biasList[[i]]$observed[[names(obsAll)[j]]]
            }
        }
    }
    for (i in seq_along(coverageFiles)) {
        ArchR:::.logMessage(sprintf("Adding Kmer Bias (%s of %s)", i, 
            length(coverageFiles)), logFile = logFile)
        obsAlli <- obsAll[[names(coverageFiles)[i]]]
        if (!identical(names(expAll), names(obsAlli))) {
            ArchR:::.logMessage("Kmer Names in Exp and Obs not Identical!", 
                logFile = logFile)
            stop("Kmer Names in Exp and Obs not Identical!")
        }
        o <- h5createGroup(coverageFiles[i], "KmerBias")
        o <- h5createGroup(coverageFiles[i], "KmerBias/Info")
        o <- h5write(obj = genome, file = coverageFiles[i], name = "KmerBias/Info/Genome")
        o <- h5write(obj = kmerLength, file = coverageFiles[i], 
            name = "KmerBias/Info/KmerLength")
        o <- h5write(obj = paste0(names(obsAlli)), file = coverageFiles[i], 
            name = "KmerBias/Kmer")
        o <- h5write(obj = obsAlli, file = coverageFiles[i], 
            name = "KmerBias/ObservedKmers")
        o <- h5write(obj = expAll, file = coverageFiles[i], name = "KmerBias/ExpectedKmers")
    }
    return(0)
}

