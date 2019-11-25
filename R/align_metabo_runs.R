#' AlignObj for analytes between a pair for DIA_metabolomics runs.
#'
#' @return A list of AlignObj. Each AlignObj contains alignment path, similarity matrix and related parameters.
#' @export
getMetaboAlignObjs <- function(analytes, runs, dataPath = ".", alignType = "hybrid",
                            query = NULL, oswMerged = TRUE, nameCutPattern = "(.*)(/)(.*)",
                            maxFdrQuery = 0.05, maxFdrLoess = 0.01, spanvalue = 0.1,
                            normalization = "mean", simMeasure = "dotProductMasked",
                            SgolayFiltOrd = 2, SgolayFiltLen = 3,
                            goFactor = 0.125, geFactor = 400,
                            cosAngleThresh = 0.3, OverlapAlignment = TRUE,
                            dotProdThresh = 0.96, gapQuantile = 0.98,
                            hardConstrain = FALSE, samples4gradient = 100,
                            samplingTime = NULL,  RSEdistFactor = 3.5){
  # Check if filter length is odd for Savitzky-Golay filter.
  if( (SgolayFiltLen %% 2) != 1){
    print("SgolayFiltLen can only be odd number.")
    return(NULL)
  }

  if(alignType == "hybrid" && is.null(samplingTime)){
    print("samplingTime cannot be NULL for hybrid alignment.")
    return(NULL)
  }

  # Get filenames from mergedosw file.
  # Check if names are consistent between osw and mzML files. Fetch run names.
  filenames <- getRunNames(dataPath, oswMerged, nameCutPattern)
  print("Following runs will be aligned:")
  rownames(filenames) <- paste0("run", 0:(nrow(filenames)-1), "")

  # Get Precursors from the query and respective chromatogram indices.
  oswFiles <- getMetaboswFiles(filenames, dataPath, peptides = NULL,  query = NULL,
                             oswMerged = oswMerged, maxFdrQuery = maxFdrQuery)
  AnalytesFound <- c()
  for(x in oswFiles){
    AnalytesFound <- x %>% .$transition_group_id %>% dplyr::union(AnalytesFound)
  }

  AnalytesFound <- intersect(analytes, AnalytesFound)
  # Report peptides that are not found
  AnalytesNotFound <- setdiff(analytes, AnalytesFound)
  if(length(AnalytesNotFound)>0){
    message(paste(AnalytesNotFound, "not found."))
    return(NULL)
  }

  ####################### Get XICs ##########################################
  runs <- filenames$runs
  names(runs) <- rownames(filenames)
  # Get Chromatogram for each peptide in each run.
  XICs <- list()
  print("Fetching Extracted-ion chromatograms from runs")
  for(i in 1:length(runs)){
    run <- names(runs)[i]
    mzmlName <- file.path(dataPath, "mzml", paste0(runs[run], ".chrom.mzML"))
    mz <- tryCatch(mzR::openMSfile(mzmlName, backend = "pwiz"),
                   error = function(cond) {
                     c$message <- paste0(c$message,
                                         "If error includes invalid cvParam accession 1002746, use FileConverter from OpenMS to decompress chromatograms")
                     stop(cond)})
    XICs_run <- lapply(1:length(AnalytesFound), function(j){
      chromIndices <- oswFiles[[i]] %>%
        dplyr::filter(transition_group_id == AnalytesFound[j]) %>% .$chromatogramIndex
      if(length(chromIndices) > 0){
        chromIndices <- as.integer(strsplit(chromIndices, split = ",")[[1]])
        XIC_group <- extractXIC_group(mz, chromIndices, SgolayFiltOrd, SgolayFiltLen)
      } else{
        XIC_group <- list()
      }
      return(XIC_group)
    })
    names(XICs_run) <- AnalytesFound
    XICs[[i]] <- XICs_run
    rm(mz)
    print(paste("Fetched Extracted-ion chromatograms from run", runs[run]))
  }
  names(XICs) <- runs

  ####################### Perfrom alignment ##########################################
  AlignObjs <- list()
  loessFits <- list()
  print("Perfroming alignment")
  for(pepIdx in 1:length(AnalytesFound)){
    AlignObjs[[pepIdx]] <- NULL
    analyte <- AnalytesFound[pepIdx]
    # Select reference run based on m-score
    minMscore <- 1; minrunIdx <- NA
    for (runIdx in 1:length(oswFiles)){
      m_score <- oswFiles[[runIdx]] %>%
        dplyr::filter(transition_group_id == analyte & peak_group_rank == 1) %>% .$m_score
      if(length(m_score) == 1){
        if(m_score < minMscore){
          minMscore <- m_score
          minrunIdx <- runIdx
        }
      }
    }
    # Select reference run based on m-score
    # Get XIC_group from reference run
    ref <- names(runs)[minrunIdx]
    exps <- setdiff(names(runs), ref)
    XICs.ref <- XICs[[runs[ref]]][[analyte]]

    # Align experiment run to reference run
    for(eXp in exps){
      # Get XIC_group from experiment run
      XICs.eXp <- XICs[[runs[eXp]]][[analyte]]
      if(length(XICs.eXp) > 0){
        # Get the loess fit for hybrid alignment
        pair <- paste(ref, eXp, sep = "_")
        if(any(pair %in% names(loessFits))){
          Loess.fit <- loessFits[[pair]]
        } else{
          df.ref <-  oswFiles[[ref]] %>% dplyr::filter(m_score <= maxFdrLoess & peak_group_rank == 1) %>%
            dplyr::select(transition_group_id, RT)
          df.eXp <-  oswFiles[[eXp]] %>% dplyr::filter(m_score <= maxFdrLoess & peak_group_rank == 1) %>%
            dplyr::select(transition_group_id, RT)
          RUNS_RT <- dplyr::inner_join(df.ref, df.eXp, by = "transition_group_id", suffix = c(".ref", ".eXp"))
          Loess.fit <- loess(RT.eXp ~ RT.ref, data = RUNS_RT,
                             span = spanvalue,
                             control=loess.control(surface="direct"))
          loessFits[[pair]] <- Loess.fit
        }
        # Set up constraints for penalizing similarity matrix
        rse <- Loess.fit$s
        noBeef <- ceiling(RSEdistFactor*rse/samplingTime)
        tVec.ref <- XICs.ref[[1]][["time"]] # Extracting time component
        tVec.eXp <- XICs.eXp[[1]][["time"]] # Extracting time component
        B1p <- predict(Loess.fit, tVec.ref[1])
        B2p <- predict(Loess.fit, tVec.ref[length(tVec.ref)])
        # Perform dynamic programming for chromatogram alignment
        intensityList.ref <- lapply(XICs.ref, `[[`, 2) # Extracting intensity values
        intensityList.eXp <- lapply(XICs.eXp, `[[`, 2) # Extracting intensity values
        Alignobj <- alignChromatogramsCpp(intensityList.ref, intensityList.eXp,
                                          alignType = alignType, tVec.ref, tVec.eXp,
                                          normalization = normalization, simType = simMeasure,
                                          B1p = B1p, B2p = B2p, noBeef = noBeef,
                                          goFactor = goFactor, geFactor = geFactor,
                                          cosAngleThresh = cosAngleThresh, OverlapAlignment = OverlapAlignment,
                                          dotProdThresh = dotProdThresh, gapQuantile = gapQuantile,
                                          hardConstrain = hardConstrain, samples4gradient = samples4gradient)
        AlignObjs[[pepIdx]] <- list()
        AlignObjs[[pepIdx]][[1]] <- Alignobj
        AlignObjs[[pepIdx]][[runs[ref]]] <- XICs.ref
        AlignObjs[[pepIdx]][[runs[eXp]]] <- XICs.eXp
        AlignObjs[[pepIdx]][[4]] <- oswFiles[[ref]] %>%
          dplyr::filter(transition_group_id == analyte) %>%
          dplyr::select(leftWidth, RT, rightWidth) %>%
          as.vector()
      }
      else {AlignObjs[[pepIdx]] <- NULL}
    }
  }
  names(AlignObjs) <- AnalytesFound
  print("Alignment done. Returning AlignObjs")
  return(AlignObjs)
}
