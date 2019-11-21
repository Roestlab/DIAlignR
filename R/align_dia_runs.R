#' Align XICs of precursors across multiple Targeted-MS runs and outputs quantitative data matrix.
#'
#' @return Saves intensity table in the current directory.
#' @importFrom dplyr %>%
#' @export
alignDIAruns <- function(dataPath, alignType = "hybrid", oswMerged = TRUE, nameCutPattern = "(.*)(/)(.*)",
                              maxFdrQuery = 0.05, maxFdrLoess = 0.01, spanvalue = 0.1, runType = "DIA_Proteomics",
                              normalization = "mean", simMeasure = "dotProductMasked",
                              SgolayFiltOrd = 4, SgolayFiltLen = 9,
                              goFactor = 0.125, geFactor = 40,
                              cosAngleThresh = 0.3, OverlapAlignment = TRUE,
                              dotProdThresh = 0.96, gapQuantile = 0.5,
                              hardConstrain = FALSE, samples4gradient = 100,
                              expRSE = 8.0, samplingTime = 3.4,  RSEdistFactor = 3.5){
  # Check if filter length is odd for Savitzky-Golay filter.
  if( (SgolayFiltLen %% 2) != 1){
    return(stop("SgolayFiltLen can only be odd number"))
  }

  # Get filenames from .merged.osw file.
  # Check if names are consistent between osw and mzML files. Fetch run names.
  filenames <- getRunNames(dataPath, oswMerged, nameCutPattern)
  message("Following runs will be aligned:")
  message(filenames[, "runs"])

  # TODO: Make this part in a separate function. Use Environment for pass-by-referene.
  # Get Precursors from the query and respectve chromatogram indices.
  oswFiles <- list()
  peptides <- c()
  for(i in 1:nrow(filenames)){
    run <- rownames(filenames)[i]
    # Get a query to search against the osw files.
    if(oswMerged == TRUE){
      oswName <- list.files(path = file.path(dataPath, "osw"), pattern="*merged.osw")
      oswName <- file.path(dataPath, "osw", oswName[1])
    } else{
      oswName <- paste0(file.path(dataPath, "osw", filenames$runs[i]), ".osw")
    }
    con <- DBI::dbConnect(RSQLite::SQLite(), dbname = oswName)
    query <- getQuery(maxFdrQuery, oswMerged, peptides = NULL, filename = filenames$filename[i], runType = "DIA_proteomics")
    x <- tryCatch(expr = DBI::dbGetQuery(con, statement = query), finally = DBI::dbDisconnect(con))

    # TODO: change how to go from osw directory to mzML directory.
    mzmlName <- file.path(dataPath, "mzml", paste0(filenames$runs[i], ".chrom.mzML"))
    mz <- tryCatch(expr = mzR::openMSfile(mzmlName, backend = "pwiz"),
                   error = function(cond) {
                     c$message <- paste0(c$message,
                      "If error includes invalid cvParam accession 1002746, use FileConverter from OpenMS to decompress chromatograms")
                     stop(cond)})
    chromHead <- mzR::chromatogramHeader(mz)
    rm(mz)
    chromHead <- chromHead %>%
      dplyr::select(chromatogramId, chromatogramIndex) %>%
      dplyr::mutate(chromatogramId = as.integer(chromatogramId))
    # TODO: Make sure that transition_id has same order across runs. IMO should be specified in query.
    x <- x %>% dplyr::left_join(chromHead, by = c("transition_id" = "chromatogramId"))
    oswFiles[[i]] <- x %>%
      dplyr::group_by(transition_group_id, peak_group_rank) %>%
      dplyr::mutate(transition_ids = paste0(transition_id, collapse = ","),
             chromatogramIndex = paste0(chromatogramIndex, collapse = ",")) %>%
      dplyr::ungroup() %>% dplyr::select(-transition_id) %>% dplyr::distinct()
    peptides <- x %>% dplyr::filter(m_score < 0.01) %>% .$transition_group_id %>% dplyr::union(peptides)
    print(paste0("Fetched chromatogram indices from ", filenames$filename[i]))
  }
  names(oswFiles) <- rownames(filenames)

  runs <- filenames$runs
  names(runs) <- rownames(filenames)
  # Collect all the pointers for each mzML file.
  print("Collecting metadata from mzML files.")
  mzPntrs <- list()
  for(mzMLindex in 1:length(runs)){
    run <- names(runs)[mzMLindex]
    filename <- file.path(dataPath, "mzml", paste0(runs[run], ".chrom.mzML"))
    mzPntrs[[mzMLindex]] <- mzR::openMSfile(filename, backend = "pwiz")
  }
  names(mzPntrs) <- names(runs)
  print("Metadata is collected from mzML files.")

  # Initilize output tables.
  rtTbl <- matrix(NA, nrow = length(peptides), ncol = length(runs))
  intesityTbl <- matrix(NA, nrow = length(peptides), ncol = length(runs))
  lwTbl <- matrix(NA, nrow = length(peptides), ncol = length(runs))
  rwTbl <- matrix(NA, nrow = length(peptides), ncol = length(runs))
  rownames(rtTbl) <- peptides; colnames(rtTbl) <- names(runs)
  rownames(intesityTbl) <- peptides; colnames(intesityTbl) <- names(runs)
  rownames(lwTbl) <- peptides; colnames(lwTbl) <- names(runs)
  rownames(rwTbl) <- peptides; colnames(rwTbl) <- names(runs)

  # Container to save loess fits
  loessFits <- list()

  print("Performing reference-based alignment.")
  start_time <- Sys.time()
  for(pepIdx in 1:length(peptides)){
    peptide <- peptides[pepIdx]
    # Select reference run based on m-score
    minMscore <- 1; minrunIdx <- NA
    for (runIdx in 1:length(oswFiles)){
      m_score <- oswFiles[[runIdx]] %>%
        dplyr::filter(transition_group_id == peptide & peak_group_rank == 1) %>% .$m_score
      if(length(m_score) == 1){
        if(m_score < minMscore){
          minMscore <- m_score
          minrunIdx <- runIdx
        }
      }
    }
    # Get the feature from reference run
    vec <- oswFiles[[minrunIdx]] %>%
      dplyr::filter(transition_group_id == peptide & peak_group_rank == 1) %>%
      dplyr::select(leftWidth, RT, rightWidth, Intensity) %>%
      as.matrix()
    lwTbl[pepIdx, minrunIdx] <- vec[1, "leftWidth"]
    rtTbl[pepIdx, minrunIdx] <- vec[1, "RT"]
    rwTbl[pepIdx, minrunIdx] <- vec[1, "rightWidth"]
    intesityTbl[pepIdx, minrunIdx] <- vec[1, "Intensity"]

    # Get XIC_group from reference run
    ref <- names(runs)[minrunIdx]
    exps <- setdiff(names(runs), ref)
    chromIndices <- oswFiles[[ref]] %>%
      dplyr::filter(transition_group_id == peptide) %>% .$chromatogramIndex
    chromIndices <- as.integer(strsplit(chromIndices, split = ",")[[1]])
    XICs.ref <- extractXIC_group(mzPntrs[[ref]], chromIndices, SgolayFiltOrd, SgolayFiltLen)
    # Align all runs to reference run
    for(eXp in exps){
      # Get XIC_group from experiment run
      chromIndices <- oswFiles[[eXp]] %>%
        dplyr::filter(transition_group_id == peptide) %>% .$chromatogramIndex
      if(length(chromIndices) > 0){
        chromIndices <- as.integer(strsplit(chromIndices, split = ",")[[1]])
        XICs.eXp <- extractXIC_group(mzPntrs[[eXp]], chromIndices)
        # Get the loess fit for hybrid alignment
        pair <- paste(ref,eXp, sep = "_")
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
        rse <- min(Loess.fit$s, expRSE)
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
        AlignedIndices <- cbind(Alignobj@indexA_aligned,
                                Alignobj@indexB_aligned,
                                Alignobj@score)
        colnames(AlignedIndices) <- c("indexAligned.ref", "indexAligned.eXp", "score")
        AlignedIndices[, 1:2][AlignedIndices[, 1:2] == 0] <- NA
        tAligned.ref <- mapIdxToTime(tVec.ref, AlignedIndices[,"indexAligned.ref"])
        tAligned.eXp <- mapIdxToTime(tVec.eXp, AlignedIndices[,"indexAligned.eXp"])
        # Map retention time from reference to eXp.
        rtTbl[pepIdx, eXp] <- tAligned.eXp[which.min(abs(tAligned.ref - rtTbl[pepIdx, ref]))]
        df <- oswFiles[[eXp]] %>% dplyr::filter(transition_group_id == peptide) %>%
          dplyr::select(leftWidth, rightWidth, RT, Intensity, peak_group_rank, m_score)
        adaptiveRT <- RSEdistFactor*rse
        df <- df %>% dplyr::filter(abs(RT - rtTbl[pepIdx, eXp]) <= adaptiveRT)
        df <- df %>% dplyr::filter(m_score < maxFdrQuery & peak_group_rank == min(peak_group_rank)) %>% as.matrix()
        if(nrow(df)==1){
          # A feature is found. Use this feature for quantification.
          lwTbl[pepIdx, eXp] <- df[1, "leftWidth"]
          rtTbl[pepIdx, eXp] <- df[1, "RT"]
          rwTbl[pepIdx, eXp] <- df[1, "rightWidth"]
          intesityTbl[pepIdx, eXp] <- df[1, "Intensity"]
        } else {
          # Feature is not found.}
        }
      }

    }
  }
  end_time <- Sys.time()
  # Report the execution time for hybrid alignment step.
  print(end_time - start_time)
  colnames(rtTbl) <- runs[colnames(rtTbl)]
  write.table(rtTbl,file = "rtTbl.csv", col.names = NA, sep = ",")
  colnames(intesityTbl) <- runs[colnames(intesityTbl)]
  write.table(intesityTbl,file = "intesityTbl.csv", col.names = NA, sep = ",")

  rm(mzPntrs)
  rm(oswFiles, loessFits)
  rm(rtTbl, intesityTbl, lwTbl, rwTbl)
  print("Data matrix is available in the current directory")
  return(1)
}


#' Align XICs of precursors across multiple Targeted-MS runs and outputs quantitative data matrix.
#'
#' @return Saves intensity table in the current directory.
#' @importFrom dplyr %>%
#' @export
alignTargetedruns <- function(dataPath, alignType = "hybrid", oswMerged = TRUE, nameCutPattern = "(.*)(/)(.*)",
                         maxFdrQuery = 0.05, maxFdrLoess = 0.01, spanvalue = 0.1, runType = "DIA_Proteomics",
                         normalization = "mean", simMeasure = "dotProductMasked",
                         SgolayFiltOrd = 4, SgolayFiltLen = 9,
                         goFactor = 0.125, geFactor = 40,
                         cosAngleThresh = 0.3, OverlapAlignment = TRUE,
                         dotProdThresh = 0.96, gapQuantile = 0.5,
                         hardConstrain = FALSE, samples4gradient = 100,
                         expRSE = 8.0, samplingTime = 3.4,  RSEdistFactor = 3.5){
  # Check if filter length is odd for Savitzky-Golay filter.
  if( (SgolayFiltLen %% 2) != 1){
    return(stop("SgolayFiltLen can only be odd number"))
  }

  # Get filenames from .merged.osw file and check if names are consistent between osw and mzML files.
  filenames <- getRunNames(dataPath, oswMerged, nameCutPattern)
  message("Following runs will be aligned:")
  message(filenames[, "runs"])

  ######### Get Precursors from the query and respectve chromatogram indices. ######
  oswFiles <- getOswFiles(dataPath, filenames, maxFdrQuery = 0.05, analyteFDR = 0.01, oswMerged = TRUE,
              peptides = NULL, runType = "DIA_proteomics")

  refAnalytes <- getAnalytesName(oswFiles, analyteFDR, commonAnalytes = FALSE)
  ######### Collect pointers for each mzML file. #######
  runs <- filenames$runs
  names(runs) <- rownames(filenames)
  # Collect all the pointers for each mzML file.
  print("Collecting metadata from mzML files.")
  mzPntrs <- list()
  for(mzMLindex in 1:length(runs)){
    run <- names(runs)[mzMLindex]
    filename <- file.path(dataPath, "mzml", paste0(runs[run], ".chrom.mzML"))
    mzPntrs[[mzMLindex]] <- mzR::openMSfile(filename, backend = "pwiz")
  }
  names(mzPntrs) <- names(runs)
  print("Metadata is collected from mzML files.")

  ######### Initilize output tables. #######
  rtTbl <- matrix(NA, nrow = length(peptides), ncol = length(runs))
  intesityTbl <- matrix(NA, nrow = length(peptides), ncol = length(runs))
  lwTbl <- matrix(NA, nrow = length(peptides), ncol = length(runs))
  rwTbl <- matrix(NA, nrow = length(peptides), ncol = length(runs))
  rownames(rtTbl) <- peptides; colnames(rtTbl) <- names(runs)
  rownames(intesityTbl) <- peptides; colnames(intesityTbl) <- names(runs)
  rownames(lwTbl) <- peptides; colnames(lwTbl) <- names(runs)
  rownames(rwTbl) <- peptides; colnames(rwTbl) <- names(runs)

  ######### Container to save loess fits.  #######
  loessFits <- list()

  print("Performing reference-based alignment.")
  start_time <- Sys.time()
  for(pepIdx in 1:length(peptides)){
    peptide <- peptides[pepIdx]
    # Select reference run based on m-score
    minMscore <- 1; minrunIdx <- NA
    for (runIdx in 1:length(oswFiles)){
      m_score <- oswFiles[[runIdx]] %>%
        dplyr::filter(transition_group_id == peptide & peak_group_rank == 1) %>% .$m_score
      if(length(m_score) == 1){
        if(m_score < minMscore){
          minMscore <- m_score
          minrunIdx <- runIdx
        }
      }
    }
    # Get the feature from reference run
    vec <- oswFiles[[minrunIdx]] %>%
      dplyr::filter(transition_group_id == peptide & peak_group_rank == 1) %>%
      dplyr::select(leftWidth, RT, rightWidth, Intensity) %>%
      as.matrix()
    lwTbl[pepIdx, minrunIdx] <- vec[1, "leftWidth"]
    rtTbl[pepIdx, minrunIdx] <- vec[1, "RT"]
    rwTbl[pepIdx, minrunIdx] <- vec[1, "rightWidth"]
    intesityTbl[pepIdx, minrunIdx] <- vec[1, "Intensity"]

    # Get XIC_group from reference run
    ref <- names(runs)[minrunIdx]
    exps <- setdiff(names(runs), ref)
    chromIndices <- oswFiles[[ref]] %>%
      dplyr::filter(transition_group_id == peptide) %>% .$chromatogramIndex
    chromIndices <- as.integer(strsplit(chromIndices, split = ",")[[1]])
    XICs.ref <- extractXIC_group(mzPntrs[[ref]], chromIndices, SgolayFiltOrd, SgolayFiltLen)
    # Align all runs to reference run
    for(eXp in exps){
      # Get XIC_group from experiment run
      chromIndices <- oswFiles[[eXp]] %>%
        dplyr::filter(transition_group_id == peptide) %>% .$chromatogramIndex
      if(length(chromIndices) > 0){
        chromIndices <- as.integer(strsplit(chromIndices, split = ",")[[1]])
        XICs.eXp <- extractXIC_group(mzPntrs[[eXp]], chromIndices)
        # Get the loess fit for hybrid alignment
        pair <- paste(ref,eXp, sep = "_")
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
        rse <- min(Loess.fit$s, expRSE)
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
        AlignedIndices <- cbind(Alignobj@indexA_aligned,
                                Alignobj@indexB_aligned,
                                Alignobj@score)
        colnames(AlignedIndices) <- c("indexAligned.ref", "indexAligned.eXp", "score")
        AlignedIndices[, 1:2][AlignedIndices[, 1:2] == 0] <- NA
        tAligned.ref <- mapIdxToTime(tVec.ref, AlignedIndices[,"indexAligned.ref"])
        tAligned.eXp <- mapIdxToTime(tVec.eXp, AlignedIndices[,"indexAligned.eXp"])
        # Map retention time from reference to eXp.
        rtTbl[pepIdx, eXp] <- tAligned.eXp[which.min(abs(tAligned.ref - rtTbl[pepIdx, ref]))]
        df <- oswFiles[[eXp]] %>% dplyr::filter(transition_group_id == peptide) %>%
          dplyr::select(leftWidth, rightWidth, RT, Intensity, peak_group_rank, m_score)
        adaptiveRT <- RSEdistFactor*rse
        df <- df %>% dplyr::filter(abs(RT - rtTbl[pepIdx, eXp]) <= adaptiveRT)
        df <- df %>% dplyr::filter(m_score < maxFdrQuery & peak_group_rank == min(peak_group_rank)) %>% as.matrix()
        if(nrow(df)==1){
          # A feature is found. Use this feature for quantification.
          lwTbl[pepIdx, eXp] <- df[1, "leftWidth"]
          rtTbl[pepIdx, eXp] <- df[1, "RT"]
          rwTbl[pepIdx, eXp] <- df[1, "rightWidth"]
          intesityTbl[pepIdx, eXp] <- df[1, "Intensity"]
        } else {
          # Feature is not found.}
        }
      }

    }
  }
  end_time <- Sys.time()
  # Report the execution time for hybrid alignment step.
  print(end_time - start_time)
  colnames(rtTbl) <- runs[colnames(rtTbl)]
  write.table(rtTbl,file = "rtTbl.csv", col.names = NA, sep = ",")
  colnames(intesityTbl) <- runs[colnames(intesityTbl)]
  write.table(intesityTbl,file = "intesityTbl.csv", col.names = NA, sep = ",")

  ######### Cleanup.  #######
  rm(mzPntrs)
  rm(oswFiles, loessFits)
  rm(rtTbl, intesityTbl, lwTbl, rwTbl)
  print("Data matrix is available in the current directory")
  return(1)
}


#' AlignObj for peptides between a pair
#'
#' @return A list of AlignObj. Each AlignObj contains alignment path, similarity matrix and related parameters.
#' @export
getAlignObjs <- function(peptides, runs, dataPath = ".", alignType = "hybrid",
                         query = NULL, oswMerged = TRUE, nameCutPattern = "(.*)(/)(.*)",
                         maxFdrQuery = 0.05, maxFdrLoess = 0.01, spanvalue = 0.1,
                         normalization = "mean", simMeasure = "dotProductMasked",
                         SgolayFiltOrd = 4, SgolayFiltLen = 9,
                         goFactor = 0.125, geFactor = 40,
                         cosAngleThresh = 0.3, OverlapAlignment = TRUE,
                         dotProdThresh = 0.96, gapQuantile = 0.5,
                         hardConstrain = FALSE, samples4gradient = 100,
                         expRSE = 8.0, samplingTime = 3.4,  RSEdistFactor = 3.5){
  if(length(runs)!= 2){
    print("For pairwise alignment, two runs are required.")
    return(NULL)
  }

  if( (SgolayFiltLen %% 2) != 1){
    print("SgolayFiltLen can only be odd number")
    return(NULL)
  }
  # Check if names are consistent between osw and mzML files. Fetch run names.
  filenames <- getRunNames(dataPath, oswMerged, nameCutPattern)
  filenames <- filenames[filenames$runs %in% runs,]
  rownames(filenames) <- paste0("run", 0:(length(runs)-1), "")

  # Get Chromatogram indices for each peptide in each run.
  oswFiles = getOswFiles(filenames, dataPath, NULL, query, oswMerged, maxFdrQuery, nameCutPattern)
  PeptidesFound <- c()
  for(x in oswFiles){
    PeptidesFound <- x %>% .$transition_group_id %>% dplyr::union(PeptidesFound)
  }

  PeptidesFound <- intersect(peptides, PeptidesFound)
  # Report peptides that are not found
  PeptidesNotFound <- setdiff(peptides, PeptidesFound)
  if(length(PeptidesNotFound)>0){
    messsage(paste(PeptidesNotFound, "not found."))
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
    XICs_run <- lapply(1:length(PeptidesFound), function(j){
      chromIndices <- oswFiles[[i]] %>%
        dplyr::filter(transition_group_id == PeptidesFound[j]) %>% .$chromatogramIndex
      if(length(chromIndices) > 0){
        chromIndices <- as.integer(strsplit(chromIndices, split = ",")[[1]])
        XIC_group <- extractXIC_group(mz, chromIndices, SgolayFiltOrd, SgolayFiltLen)
      } else{
        XIC_group <- list()
      }
      return(XIC_group)
    })
    names(XICs_run) <- PeptidesFound
    XICs[[i]] <- XICs_run
    rm(mz)
    print(paste("Fetched Extracted-ion chromatograms from run", runs[run]))
  }
  names(XICs) <- runs

  ####################### Perfrom alignment ##########################################
  AlignObjs <- list()
  loessFits <- list()
  print("Perfroming alignment")
  for(pepIdx in 1:length(peptides)){
    peptide <- peptides[pepIdx]
    # Select reference run based on m-score
    minMscore <- 1; minrunIdx <- NA
    for (runIdx in 1:length(oswFiles)){
      m_score <- oswFiles[[runIdx]] %>%
        dplyr::filter(transition_group_id == peptide & peak_group_rank == 1) %>% .$m_score
      if(length(m_score) == 1){
        if(m_score < minMscore){
          minMscore <- m_score
          minrunIdx <- runIdx
        }
      }
    }
    # Get XIC_group from reference run
    ref <- names(runs)[minrunIdx]
    exps <- setdiff(names(runs), ref)
    XICs.ref <- XICs[[runs[ref]]][[peptide]]

    # Align experiment run to reference run
    for(eXp in exps){
      # Get XIC_group from experiment run
      XICs.eXp <- XICs[[runs[eXp]]][[peptide]]
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
        rse <- min(Loess.fit$s, expRSE)
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
        AlignObjs[[pepIdx]][[4]] <- oswFiles[[minrunIdx]] %>%
          dplyr::filter(transition_group_id == peptide & peak_group_rank == 1) %>%
          dplyr::select(leftWidth, RT, rightWidth) %>%
          as.vector()
      }
      else {AlignObjs[[pepIdx]] <- NULL}
    }
  }
  names(AlignObjs) <- peptides
  print("Alignment done. Returning AlignObjs")
  return(AlignObjs)
}
