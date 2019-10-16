#' Align XICs of precursors across multiple Targeted-MS runs and outputs quantitative data matrix.
#'
#' @return Saves intensity table in the current directory.
#' @importFrom dplyr %>%
#' @export
alignTragetedRuns <- function(dataPath, alignType = "hybrid", oswMerged = TRUE, nameCutPattern = "(.*)(/)(.*)",
                              maxFdrQuery = 0.05, maxFdrLoess = 0.01, spanvalue = 0.1,
                              normalization = "mean", simMeasure = "dotProductMasked",
                              SgolayFiltOrd = 4, SgolayFiltLen = 9,
                              goFactor = 0.125, geFactor = 40,
                              cosAngleThresh = 0.3, OverlapAlignment = TRUE,
                              dotProdThresh = 0.96, gapQuantile = 0.5,
                              hardConstrain = FALSE, samples4gradient = 100,
                              expRSE = 8.0, samplingTime = 3.4,  RSEdistFactor = 3.5){
  # Check if filter length is odd for Savitzky-Golay filter.
  if( (SgolayFiltLen %% 2) != 1){
    print("SgolayFiltLen can only be odd number")
    return(NULL)
  }

  # Get filenames from mergedosw file.
  # Check if names are consistent between osw and mzML files. Fetch run names.
  filenames <- getRunNames(dataPath, oswMerged, nameCutPattern)
  print("Following runs will be aligned:")
  rownames(filenames) <- paste0("run", 0:(nrow(filenames)-1), "")
  print(filenames[, "runs"])

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
    query <- getQuery(maxFdrQuery, oswMerged, filename = filenames$filename[i])
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

#' This is a query that will be used to fetch information from osw files.
#'
#' @return SQL query to be searched.
getQuery <- function(maxFdrQuery, oswMerged = TRUE, peptides = NULL, filename = NULL){
  if(is.null(peptides)){
    selectPeptide <- ""
  } else{
    selectPeptide <- paste0(" AND transition_group_id IN ('", paste(peptides,collapse="','"),"')")
  }

  if(oswMerged){
    matchFilename <- paste0(" AND RUN.FILENAME ='", filename,"'")
  } else{
    matchFilename <- ""
  }

  query <- paste0("SELECT PEPTIDE.MODIFIED_SEQUENCE || '_' || PRECURSOR.CHARGE AS transition_group_id,
  RUN.FILENAME AS filename,
  FEATURE.EXP_RT AS RT,
  FEATURE.DELTA_RT AS delta_rt,
  PRECURSOR.LIBRARY_RT AS assay_RT,
  FEATURE_MS2.AREA_INTENSITY AS Intensity,
  FEATURE.LEFT_WIDTH AS leftWidth,
  FEATURE.RIGHT_WIDTH AS rightWidth,
  SCORE_MS2.RANK AS peak_group_rank,
  SCORE_MS2.QVALUE AS m_score,
  TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID AS transition_id
  FROM PRECURSOR
  INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR.ID = PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID AND PRECURSOR.DECOY=0
  INNER JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID
  INNER JOIN FEATURE ON FEATURE.PRECURSOR_ID = PRECURSOR.ID
  INNER JOIN RUN ON RUN.ID = FEATURE.RUN_ID
  INNER JOIN TRANSITION_PRECURSOR_MAPPING ON TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID = PRECURSOR.ID
  LEFT JOIN FEATURE_MS2 ON FEATURE_MS2.FEATURE_ID = FEATURE.ID
  LEFT JOIN SCORE_MS2 ON SCORE_MS2.FEATURE_ID = FEATURE.ID
  WHERE SCORE_MS2.QVALUE < ", maxFdrQuery, selectPeptide, matchFilename, "
  ORDER BY transition_group_id,
  peak_group_rank;")
  return(query)
}

#' Extract XICs of all transitions requested in chromIndices.
#'
#' @return A list of data-frames. Each data frame has elution time and intensity of fragment-ion XIC.
extractXIC_group <- function(mz, chromIndices, SgolayFiltOrd = 4, SgolayFiltLen = 9){
  XIC_group <- lapply(1:length(chromIndices), function(i) {
    rawChrom <- mzR::chromatograms(mz, chromIndices[i])
    # Savitzky-Golay filter to smooth chromatograms, filter order p = 3, filter length n = 13
    rawChrom[,2] <- signal::sgolayfilt(rawChrom[,2], p = SgolayFiltOrd, n = SgolayFiltLen)
    return(rawChrom)
  })
  return(XIC_group)
}


#' Get names of all runs
#'
#' @return A vector with names of all runs.
getRunNames <- function(dataPath, oswMerged = TRUE, nameCutPattern = "(.*)(/)(.*)"){
  # Get names of osw files.
  if(oswMerged == FALSE){
    print("Looking for .osw files.")
    temp <- list.files(path = file.path(dataPath, "osw"), pattern="*.osw")
    filenames <- sapply(1:length(temp), function(i){
      con <- DBI::dbConnect(RSQLite::SQLite(), dbname = file.path(dataPath, "osw", temp[i]))
      query <- "SELECT DISTINCT RUN.FILENAME AS filename FROM RUN"
      tryCatch(expr = DBI::dbGetQuery(con, statement = query), finally = DBI::dbDisconnect(con))
    })
    filenames <- as.data.frame(unique(unlist(filenames)))
    colnames(filenames) <-  c("filename")
    filenames$filename <- as.character(filenames$filename)
    print(paste0(nrow(filenames), " .osw files are found."))
  } else{
    print("Looking for .merged.osw file.")
    temp <- list.files(path = file.path(dataPath, "osw"), pattern="*merged.osw")
    con <- DBI::dbConnect(RSQLite::SQLite(), dbname = file.path(dataPath, "osw", temp[1]))
    query <- "SELECT DISTINCT RUN.FILENAME AS filename FROM RUN"
    filenames <- tryCatch(expr = DBI::dbGetQuery(con, statement = query), finally = DBI::dbDisconnect(con))
    print(paste0(nrow(filenames), " are in", temp[1], "file"))
  }
  # Get names of mzml files.
  filenames$runs <- sapply(filenames[,"filename"], function(x) strsplit(x, split = ".mzML")[[1]][1])
  filenames$runs <- sapply(filenames$runs, function(x) gsub(nameCutPattern, replacement = "\\3", x))
  temp <- list.files(path = file.path(dataPath, "mzml"), pattern="*.chrom.mzML")
  print(paste0(length(temp), " .chrom.mzML files are found."))
  runs2 <- sapply(temp, function(x) strsplit(x, split = ".chrom.mzML")[[1]][1])
  # Check if osw files have corresponding mzML file.
  runs <- intersect(filenames$runs, runs2)
  filenames <- filenames[filenames$runs %in% runs,]
  rownames(filenames) <- NULL
  if(length(runs) == 0){
    print("Number of osw files and mzml files aren't matching.")
    print("Check if you have correct file names.")
    print("Files in mzml directory should end with .chrom.mzML")
    print("Files in osw directory should end with .osw or .merged.osw")
    return(NULL)
  } else{
    return(filenames)
  }
}

#' Get list of peptides and their chromatogram indices.
#'
#' @importFrom dplyr %>%
#' @return A list of data-frames.
getOswFiles <- function(peptides, filenames, dataPath = ".", query = NULL,
                        oswMerged = TRUE, maxFdrQuery = 1.0, nameCutPattern = "(.*)(/)(.*)"){
  # Get Chromatogram indices for each peptide in each run.
  print("Getting chromatogram indices for each peptide in each run")
  oswFiles <- list()
  for(i in 1:nrow(filenames)){
    run <- rownames(filenames)[i]
    if(oswMerged == TRUE){
      oswName <- list.files(path = file.path(dataPath, "osw"), pattern="*merged.osw")
      oswName <- file.path(dataPath, "osw", oswName[1])
    } else{
      oswName <- paste0(file.path(dataPath, "osw", filenames$runs[i]), ".osw")
    }
    con <- DBI::dbConnect(RSQLite::SQLite(), dbname = oswName)
    # Get a query to search against the osw files.
    if(is.null(query)){
      query <- getQuery(maxFdrQuery, oswMerged, peptides, filename = filenames$filename[i])
    }
    x <- tryCatch(expr = DBI::dbGetQuery(con, statement = query), finally = DBI::dbDisconnect(con))

    # Read chromatogram indices from mzML file
    mzmlName <- file.path(dataPath, "mzml", paste0(filenames$runs[i], ".chrom.mzML"))
    mz <- tryCatch(mzR::openMSfile(mzmlName, backend = "pwiz"),
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
    oswFiles[[i]] <- x %>% dplyr::filter(peak_group_rank == 1) %>%
      dplyr::group_by(transition_group_id, peak_group_rank) %>%
      dplyr::mutate(transition_ids = paste0(transition_id, collapse = ","),
                    chromatogramIndex = paste0(chromatogramIndex, collapse = ",")) %>%
      dplyr::ungroup() %>% dplyr::select(-transition_id, -peak_group_rank, -assay_RT, -delta_rt) %>% dplyr::distinct()
    print(paste("Fetched chromatogram indices from run", filenames$runs[i]))
  }
  names(oswFiles) <- rownames(filenames)
  print("Fetched chromatogram indices for each peptide in each run")
  return(oswFiles)
}

#' Get XICs for a list of peptides
#'
#' @return A list of list. Each list contains XICs for that run.
#' @importFrom dplyr %>%
#' @export
getXICs <- function(peptides, runs, dataPath = ".", maxFdrQuery = 1.0,
                    SgolayFiltOrd = 4, SgolayFiltLen = 9,
                    query = NULL, oswMerged = TRUE, nameCutPattern = "(.*)(/)(.*)"){
  if( (SgolayFiltLen %% 2) != 1){
    print("SgolayFiltLen can only be odd number")
    return(NULL)
  }
  # Check if names are consistent between osw and mzML files. Fetch run names.
  filenames <- getRunNames(dataPath, oswMerged, nameCutPattern)
  filenames <- filenames[filenames$runs %in% runs,]
  rownames(filenames) <- paste0("run", 0:(length(runs)-1), "")

  # Get Chromatogram indices for each peptide in each run.
  oswFiles = getOswFiles(peptides, filenames, dataPath, query, oswMerged, maxFdrQuery, nameCutPattern)
  PeptidesFound <- c()
  for(x in oswFiles){
    PeptidesFound <- x %>% .$transition_group_id %>% dplyr::union(PeptidesFound)
  }

  # Report peptides that are not found
  PeptidesNotFound <- setdiff(peptides, PeptidesFound)
  if(length(PeptidesNotFound)>0){
    messsage(paste(PeptidesNotFound, "not found."))
  }

  runs <- filenames$runs
  names(runs) <- rownames(filenames)
  # Get Chromatogram for each peptide in each run.
  XICs <- list()
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
  }
  names(XICs) <- runs
  return(XICs)
}
