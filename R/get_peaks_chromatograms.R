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

#' Extract XICs of all transitions requested in chromIndices.
#'
#' @return A list of list of data-frames. Each data frame has elution time and intensity of fragment-ion XIC.
#' @export
getXICs4AlignObj <- function(dataPath, runs, oswFiles, analytes,
                             SgolayFiltOrd, SgolayFiltLen){
  mzPntrs <- getMZMLpointers(dataPath, runs)
  XICs <- vector("list", length(runs))
  names(XICs) <- names(runs)
  for(i in seq_along(runs)){
    runname = names(runs)[i]
    XICs[[i]] <- lapply(1:length(analytes), function(j){
      chromIndices <- selectChromIndices(oswFiles, runname = runname, analyte = analytes[j])
      if(is.null(chromIndices)){
        warning("Chromatogram indices for ", analyte, " are missing in ", runs[[runname]])
        message("Skipping ", analyte)
        XIC_group <- NULL
      } else {
        XIC_group <- extractXIC_group(mzPntrs[[runname]], chromIndices, SgolayFiltOrd, SgolayFiltLen)
      }
      XIC_group
    })
    names(XICs[[i]]) <- analytes
  }
  rm(mzPntrs)
  XICs
}

#' Extract XICs of all transitions requested in chromIndices.
#'
#' @return A list of list of data-frames. Each data frame has elution time and intensity of fragment-ion XIC.
getAlignObj <- function(XICs.ref, XICs.eXp, Loess.fit, adaptiveRT, samplingTime,
                        normalization, simType, goFactor, geFactor,
                        cosAngleThresh, OverlapAlignment,
                        dotProdThresh, gapQuantile, hardConstrain,
                        samples4gradient, objType = "light"){
  # Set up constraints for penalizing similarity matrix
  noBeef <- ceiling(adaptiveRT/samplingTime)
  tVec.ref <- XICs.ref[[1]][["time"]] # Extracting time component
  tVec.eXp <- XICs.eXp[[1]][["time"]] # Extracting time component
  B1p <- predict(Loess.fit, tVec.ref[1])
  B2p <- predict(Loess.fit, tVec.ref[length(tVec.ref)])
  # Perform dynamic programming for chromatogram alignment
  intensityList.ref <- lapply(XICs.ref, `[[`, 2) # Extracting intensity values
  intensityList.eXp <- lapply(XICs.eXp, `[[`, 2) # Extracting intensity values
  AlignObj <- alignChromatogramsCpp(intensityList.ref, intensityList.eXp,
                                    alignType = "hybrid", tVec.ref, tVec.eXp,
                                    normalization = normalization, simType = simType,
                                    B1p = B1p, B2p = B2p, noBeef = noBeef,
                                    goFactor = goFactor, geFactor = geFactor,
                                    cosAngleThresh = cosAngleThresh, OverlapAlignment = OverlapAlignment,
                                    dotProdThresh = dotProdThresh, gapQuantile = gapQuantile,
                                    hardConstrain = hardConstrain, samples4gradient = samples4gradient,
                                    objType = objType)
  AlignObj
}

#' Extract XICs of all transitions requested in chromIndices.
#'
#' @return A list of list of data-frames. Each data frame has elution time and intensity of fragment-ion XIC.
#' @export
getMappedRT <- function(refRT, XICs.ref, XICs.eXp, Loess.fit, alignType, adaptiveRT, samplingTime,
                        normalization, simMeasure, goFactor, geFactor, cosAngleThresh,
                        OverlapAlignment, dotProdThresh, gapQuantile, hardConstrain,
                        samples4gradient, objType = "light"){
  AlignObj <- getAlignObj(XICs.ref, XICs.eXp, Loess.fit, adaptiveRT, samplingTime,
                          normalization, simType = simMeasure, goFactor, geFactor,
                          cosAngleThresh, OverlapAlignment,
                          dotProdThresh, gapQuantile, hardConstrain, samples4gradient, objType)
  tVec.ref <- XICs.ref[[1]][["time"]] # Extracting time component
  tVec.eXp <- XICs.eXp[[1]][["time"]] # Extracting time component
  eXpRT <- mappedRTfromAlignObj(refRT, tVec.ref, tVec.eXp, AlignObj)
  eXpRT
}


#' Get XICs for a list of peptides
#'
#' @return A list of list. Each list contains XICs for that run.
#' @importFrom dplyr %>%
#' @export
getXICs <- function(analytes, runs, dataPath = ".", maxFdrQuery = 1.0,
                    SgolayFiltOrd = 4, SgolayFiltLen = 9, runType = "DIA_proteomics",
                    oswMerged = TRUE, nameCutPattern = "(.*)(/)(.*)"){
  if( (SgolayFiltLen %% 2) != 1){
    print("SgolayFiltLen can only be odd number")
    return(NULL)
  }
  # Get filenames from .merged.osw file and check if names are consistent between osw and mzML files.
  filenames <- getRunNames(dataPath, oswMerged, nameCutPattern)
  filenames <- filenames[filenames$runs %in% runs,]

  # Get Chromatogram indices for each peptide in each run.
  oswFiles = getOswFiles(dataPath, filenames, maxFdrQuery = maxFdrQuery, analyteFDR = 1.00,
                         oswMerged = oswMerged, analytes = analytes, runType = runType)
  refAnalytes <- getAnalytesName(oswFiles, commonAnalytes = FALSE)
  analytesFound <- intersect(analytes, refAnalytes)
  analytesNotFound <- setdiff(analytes, analytesFound)
  if(length(analytesNotFound)>0){
    message(paste(analytesNotFound, "not found."))
  }

  ####################### Get XICs ##########################################
  runs <- filenames$runs
  names(runs) <- rownames(filenames)
  # Get Chromatogram for each peptide in each run.
  message("Fetching Extracted-ion chromatograms from runs")
  XICs <- getXICs4AlignObj(dataPath, runs, oswFiles, analytesFound,
                           SgolayFiltOrd, SgolayFiltLen)
  XICs
}


#' Get list of peptides and their chromatogram indices.
#'
#' @importFrom dplyr %>%
#' @return A list of data-frames.
#' @export
getMRMoswFiles <- function(filenames, dataPath = ".", peptides = NULL, query = NULL,
                        oswMerged = FALSE, maxFdrQuery = 1.0){
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
      query <- getQuery(maxFdrQuery, oswMerged, peptides,
                        filename = filenames$filename[i], runType = "MRM_Proteomics")
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
      dplyr::select(chromatogramId, chromatogramIndex, productIsolationWindowTargetMZ) %>%
      dplyr::mutate(chromatogramId = as.integer(chromatogramId),
                    productIsolationWindowTargetMZ = as.integer(floor(productIsolationWindowTargetMZ)))
    # chromatogramId has some characters which may cause some warnings.
    # TODO: Make sure that transition_id has same order across runs. IMO should be specified in query.
    x <- x %>% dplyr::left_join(chromHead, by = c("transition_id" = "chromatogramId"))
    oswFiles[[i]] <- x %>% dplyr::group_by(transition_group_id, RT) %>%
      dplyr::mutate(transition_ids = paste0(transition_id, collapse = ","),
                    chromatogramIndex = paste0(chromatogramIndex, collapse = ","),
                    productMZ = paste0(productIsolationWindowTargetMZ, collapse = ",")) %>%
      dplyr::ungroup() %>%
      dplyr::select(-transition_id, -assay_RT, -delta_rt, -productIsolationWindowTargetMZ) %>% dplyr::distinct()
    print(paste("Fetched chromatogram indices from run", filenames$runs[i]))
  }
  names(oswFiles) <- rownames(filenames)
  print(paste("Fetched chromatogram indices from run", filenames$runs[i]))
  return(oswFiles)
}


#' Get XICs for a list of peptides from MRM runs.
#'
#' @return A list of list. Each list contains XICs for that run.
#' @importFrom dplyr %>%
getMRMXICs <- function(peptides, runs, dataPath = ".", maxFdrQuery = 1.0,
                    SgolayFiltOrd = 2, SgolayFiltLen = 3,
                    query = NULL, oswMerged = FALSE, nameCutPattern = "(.*)(/)(.*)"){
  if( (SgolayFiltLen %% 2) != 1){
    print("SgolayFiltLen can only be odd number")
    return(NULL)
  }
  # Check if names are consistent between osw and mzML files. Fetch run names.
  filenames <- getRunNames(dataPath, oswMerged, nameCutPattern)
  filenames <- filenames[filenames$runs %in% runs,]
  rownames(filenames) <- paste0("run", 0:(length(runs)-1), "")

  # Get Chromatogram indices for each peptide in each run.
  oswFiles = getMRMoswFiles(filenames, dataPath, peptides, query, oswMerged, maxFdrQuery)
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
  return(XICs)
}

