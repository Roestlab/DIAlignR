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





#' Get list of peptides and their chromatogram indices.
#'
#' @importFrom dplyr %>%
#' @return A list of data-frames.
#' @export
get1OswFiles <- function(filenames, dataPath = ".", peptides = NULL,  query = NULL,
                        oswMerged = TRUE, maxFdrQuery = 1.0, nameCutPattern = "(.*)(/)(.*)",
                        runType = "DIA_Proteomics"){
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
                        filename = filenames$filename[i], runType = runType)
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
      dplyr::ungroup() %>% dplyr::select(-transition_id, -assay_RT, -delta_rt) %>% dplyr::distinct()
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
  oswFiles = getOswFiles(filenames, dataPath, peptides, query, oswMerged, maxFdrQuery, nameCutPattern)
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


#' Get list of peptides and their chromatogram indices.
#'
#' @importFrom dplyr %>%
#' @return A list of data-frames.
#' @export
getMRMoswFiles <- function(filenames, dataPath = ".", peptides = NULL,  query = NULL,
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
#' @export
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


#' Get list of analytes and their chromatogram indices.
#'
#' @importFrom dplyr %>%
#' @return A list of data-frames.
#' @export
getMetaboswFiles <- function(filenames, dataPath = ".", peptides = NULL,  query = NULL,
                           oswMerged = FALSE, maxFdrQuery = 1.0){
  # Get Chromatogram indices for each analyte in each run.
  print("Getting chromatogram indices for each analyte in each run")
  oswFiles <- list()
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
    query <- getQuery(maxFdrQuery, oswMerged, peptides = NULL,
                      filename = filenames$filename[i], runType = "DIA_Metabolomics")
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
    print(paste0("Fetched chromatogram indices from ", filenames$filename[i]))
  }
  names(oswFiles) <- rownames(filenames)
  print(paste("Fetched chromatogram indices from run", filenames$runs[i]))
  return(oswFiles)
}


#' Get names of analytes found in all runs.
#'
#' @importFrom dplyr %>%
#' @return A vector of strings.
#' @export
getAnalytesName <- function(runs, dataPath = ".",
                            query = NULL, oswMerged = TRUE, nameCutPattern = "(.*)(/)(.*)",
                            maxFdrQuery = 0.05){
  filenames <- getRunNames(dataPath, oswMerged, nameCutPattern)
  filenames <- filenames[filenames$runs %in% runs,]
  rownames(filenames) <- paste0("run", 0:(length(runs)-1), "")

  # Get Precursors from the query and respective chromatogram indices.
  oswFiles <- getMetaboswFiles(filenames, dataPath, peptides = NULL,  query = NULL,
                               oswMerged = oswMerged, maxFdrQuery = maxFdrQuery)
  AnalytesFound <- oswFiles[[1]] %>% .$transition_group_id
  for(x in oswFiles){
    AnalytesFound <- x %>% .$transition_group_id %>% dplyr::intersect(AnalytesFound)
  }
  return(AnalytesFound)
}
