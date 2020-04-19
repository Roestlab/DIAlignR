#' Fetch features of analytes
#'
#' Get a data-frame of analytes' transition_group_ids, their OpenSwath features, chromatogram indices and associated FDR-scores.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-13
#' @param oswName (char) path to the osw file.
#' @param maxFdrQuery (numeric) A numeric value between 0 and 1. It is used to filter features from osw file which have SCORE_MS2.QVALUE less than itself.
#' @param oswMerged (logical) TRUE for experiment-wide FDR and FALSE for run-specific FDR by pyprophet.
#' @param analytes (vector of strings) transition_group_ids for which features are to be extracted. analyteInGroupLabel must be set according the pattern used here.
#' @param filename (data-frame) Should be from the RUN.FILENAME column from osw files.
#' @param runType (char) This must be one of the strings "DIA_proteomics", "DIA_Metabolomics".
#' @param analyteInGroupLabel (logical) TRUE for getting analytes as PRECURSOR.GROUP_LABEL from osw file.
#'  FALSE for fetching analytes as PEPTIDE.MODIFIED_SEQUENCE and PRECURSOR.CHARGE from osw file.
#' @return (data-frames) Data-frame has following columns:
#' \item{transition_group_id}{(string) it is either fetched from PRECURSOR.GROUP_LABEL or a combination of PEPTIDE.MODIFIED_SEQUENCE and PRECURSOR.CHARGE from osw file.}
#' \item{filename}{(string) as mentioned in RUN table of osw files.}
#' \item{RT}{(numeric) retention time as in FEATURE.EXP_RT of osw files.}
#' \item{delta_rt}{(numeric) as in FEATURE.DELTA_RT of osw files.}
#' \item{assay_RT}{(numeric) library retention time as in PRECURSOR.LIBRARY_RT of osw files.}
#' \item{Intensity}{(numeric) peak intensity as in FEATURE_MS2.AREA_INTENSITY of osw files.}
#' \item{leftWidth}{(numeric) as in FEATURE.LEFT_WIDTH of osw files.}
#' \item{rightWidth}{(numeric) as in FEATURE.RIGHT_WIDTH of osw files.}
#' \item{peak_group_rank}{(integer) rank of each feature associated with transition_group_id.}
#' \item{m_score}{(numeric) q-value of each feature associated with transition_group_id.}
#' \item{transition_id}{(integer) fragment-ion ID associated with transition_group_id. This is matched with chromatogram ID in mzML file.}
#'
#' @seealso \code{\link{getRunNames}}
#' @keywords internal
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' filenames <- DIAlignR::getRunNames(dataPath = dataPath)
#' oswName <- paste0(dataPath,"/osw/merged.osw")
#' \dontrun{
#' analytesInfo <- fetchAnalytesInfo(oswName, maxFdrQuery = 0.05, oswMerged = TRUE,
#'  analytes = c("19051_KLIVTSEGC[160]FK/2"), filename = filenames$filename[2],
#'   runType = "DIA_proteomics", analyteInGroupLabel = TRUE)
#' analytesInfo <- fetchAnalytesInfo(oswName, maxFdrQuery = 0.05, oswMerged = TRUE,
#'  analytes = c("IHFLSPVRPFTLTPGDEEESFIQLITPVR_3"), filename = filenames$filename[3],
#'   runType = "DIA_proteomics", analyteInGroupLabel = FALSE)
#' }
fetchAnalytesInfo <- function(oswName, maxFdrQuery, oswMerged,
                              analytes, filename, runType, analyteInGroupLabel = FALSE){
  # Establish a connection of SQLite file.
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = oswName)
  # Generate a query.
  query <- getQuery(maxFdrQuery, oswMerged, analytes = analytes,
                    filename = filename, runType = runType,
                    analyteInGroupLabel = analyteInGroupLabel)
  # Run query to get peptides, their coordinates and scores.
  analytesInfo <- tryCatch(expr = DBI::dbGetQuery(con, statement = query),
                           finally = DBI::dbDisconnect(con))
  analytesInfo
}


#' Fetch analytes from OSW file
#'
#' Get a data-frame of analytes, their chromatogram indices and associated FDR-scores.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-13
#' @param dataPath (char) path to mzml and osw directory.
#' @param filenames (data-frame) column "filename" contains RUN table from osw files. column "runs" contain respective mzML names without extension.
#' To get filenames use DIAlignR::getRunNames function.
#' @param oswMerged (logical) TRUE for experiment-wide FDR and FALSE for run-specific FDR by pyprophet.
#' @param analyteInGroupLabel (logical) TRUE for getting analytes as PRECURSOR.GROUP_LABEL from osw file.
#'  FALSE for fetching analytes as PEPTIDE.MODIFIED_SEQUENCE and PRECURSOR.CHARGE from osw file.
#' @param maxFdrQuery (numeric) A numeric value between 0 and 1. It is used to filter features from osw file which have SCORE_MS2.QVALUE less than itself.
#' @param runType (char) This must be one of the strings "DIA_proteomics", "DIA_Metabolomics".
#' @return (A list of data-frames) Each data-frame has following columns:
#' \item{transition_group_id}{(string) it is either fetched from PRECURSOR.GROUP_LABEL or a combination of PEPTIDE.MODIFIED_SEQUENCE and PRECURSOR.CHARGE from osw file.}
#' \item{filename}{(string) as mentioned in RUN table of osw files.}
#' \item{peak_group_rank}{(integer) rank of each feature associated with transition_group_id.}
#' \item{m_score}{(numeric) q-value of each feature associated with transition_group_id.}
#' \item{transition_id}{(integer) fragment-ion ID associated with transition_group_id. This is matched with chromatogram ID in mzML file.}
#'
#' @keywords internal
#' @seealso \code{\link{getRunNames}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' filenames <- getRunNames(dataPath = dataPath)
#' \dontrun{
#' oswFiles <- getOswAnalytes(dataPath = dataPath, filenames = filenames,
#'  analyteInGroupLabel = TRUE)
#' oswFiles[["run0"]][1,]
#' oswFiles <- getOswAnalytes(dataPath = dataPath, filenames = filenames,
#'  analyteInGroupLabel = FALSE)
#' oswFiles[["run0"]][1,]
#' }
getOswAnalytes <- function(fileInfo, oswMerged = TRUE, analyteInGroupLabel = FALSE,
                           maxFdrQuery = 0.05, runType  = "DIA_proteomics"){
  oswFiles <- list()
  for(i in 1:nrow(fileInfo)){
    # Get a query to search against the osw files.
    oswName <- as.character(fileInfo[["featureFile"]][[i]])
    # Establish a connection of SQLite file.
    con <- DBI::dbConnect(RSQLite::SQLite(), dbname = oswName)
    # Generate a query.
    query <- getAnalytesQuery(maxFdrQuery = maxFdrQuery, oswMerged = oswMerged,
                      filename = fileInfo$spectraFile[i], runType = runType,
                      analyteInGroupLabel = analyteInGroupLabel)
    # Run query to get peptides, their coordinates and scores.
          oswAnalytes <- tryCatch(expr = DBI::dbGetQuery(con, statement = query),
                             finally = DBI::dbDisconnect(con))
    oswFiles[[i]] <- oswAnalytes
  }
  names(oswFiles) <- rownames(fileInfo)
  oswFiles
}



#' Get precursors from a feature file
#'
#' Get a data-frame of analytes' transition_group_id, transition_ids, peptide_id and amino-acid sequences.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-04-04
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @param filename (string) Should be from the RUN.FILENAME column from osw files.
#' @param runType (string) This must be one of the strings "DIA_proteomics", "DIA_Metabolomics".
#' @return (data-frames) Data-frame has following columns:
#' \item{transition_group_id}{(integer) a unique id for each precursor.}
#' \item{transition_id}{(list) fragment-ion ID associated with transition_group_id. This is matched with chromatogram ID in mzML file.}
#' \item{peptide_id}{(integer) a unique id for each peptide. A peptide can have multiple precursors.}
#' \item{sequence}{(string) amino-acid sequence of the precursor with possible modifications.}
#' \item{charge}{(integer) charge on the precursor.}
#' \item{group_label}{(string) TODO Figure it out.}
#'
#' @seealso \code{\link{getRunNames}, \link{getPrecursors}, \link{getPrecursorsQuery}}
#' @keywords internal
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' filename <- paste0(dataPath,"/osw/merged.osw")
#' \dontrun{
#' precursorsInfo <- fetchPrecursorsInfo(filename, runType = "DIA_proteomics")
#' dim(precursorsInfo) # 303  6
#' }
fetchPrecursorsInfo <- function(filename, runType, selectIDs = NULL){
  # Establish a connection of SQLite file.
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = as.character(filename))
  # Generate a query.

  if(is.null(selectIDs)){
    query <- getPrecursorsQuery(runType)
  } else{
    query <- getPrecursorsQueryID(selectIDs, runType)
  }
  # Run query to get peptides, their coordinates and scores.
  precursorsInfo <- tryCatch(expr = { output <- DBI::dbSendQuery(con, statement = query)
                                        DBI::dbFetch(output)},
                           finally = {DBI::dbClearResult(output)
                             DBI::dbDisconnect(con)})
  # Each precursor has only one row. tidyr::nest creates a tibble object that is twice as heavy to regular list.
  precursorsInfo <- dplyr::group_by(precursorsInfo, .data$transition_group_id, .data$peptide_id, .data$sequence, .data$charge, .data$group_label) %>%
    dplyr::summarise(transition_ids = base::list(.data$transition_id)) %>% dplyr::ungroup() %>% as.data.frame()
  precursorsInfo
}

#' Get precursors from all feature files
#'
#' Get a data-frame of analytes' transition_group_id, transition_ids, peptide_id and amino-acid sequences.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-04-06
#' @importFrom dplyr %>%
#' @param fileInfo (data-frame) Output of DIAlignR::getRunNames function.
#' @param oswMerged (logical) TRUE for experiment-wide FDR and FALSE for run-specific FDR by pyprophet.
#' @param runType (char) This must be one of the strings "DIA_proteomics", "DIA_Metabolomics".
#' @return (data-frames) A data-frame having following columns:
#' \item{transition_group_id}{(integer) a unique id for each precursor.}
#' \item{transition_id}{(list) fragment-ion ID associated with transition_group_id. This is matched with chromatogram ID in mzML file.}
#' \item{peptide_id}{(integer) a unique id for each peptide. A peptide can have multiple precursors.}
#' \item{sequence}{(string) amino-acid sequence of the precursor with possible modifications.}
#' \item{charge}{(integer) charge on the precursor.}
#' \item{group_label}{(string) TODO Figure it out.}
#'
#' @seealso \code{\link{getRunNames}, \link{fetchPrecursorsInfo}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' fileInfo <- DIAlignR::getRunNames(dataPath = dataPath)
#' \dontrun{
#' precursorsInfo <- getPrecursors(fileInfo, oswMerged = TRUE, runType = "DIA_proteomics")
#' dim(precursorsInfo) # 322  6
#' }
#' @export
getPrecursors <- function(fileInfo, oswMerged = TRUE, runType = "DIA_proteomics"){
  if(oswMerged == TRUE){
    # Get precursor information from merged.osw file
    oswName <- unique(fileInfo[["featureFile"]])
    precursors <- fetchPrecursorsInfo(oswName, runType)
  } else {
    # Iterate over each file and collect precursor information
    precursors <- data.frame("transition_group_id" = integer())
    for(i in 1:nrow(fileInfo)){
      oswName <- fileInfo[["featureFile"]][[i]]
      temp <- fetchPrecursorsInfo(oswName, runType)
      precursors <- merge(precursors, temp, by = c("transition_group_id", all.x = TRUE, all.y = TRUE))
    }
  }
  message(paste0(nrow(precursors), " precursors are found."))
  precursors
}


#' Find precursors given their IDs
#'
#' Get a data-frame of analytes' transition_group_id, transition_ids, peptide_id and amino-acid sequences.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2019-04-06
#' @importFrom dplyr %>%
#' @param analytes (integer) a vector of integers.
#' @param fileInfo (data-frame) output of getRunNames function.
#' @param oswMerged (logical) TRUE for experiment-wide FDR and FALSE for run-specific FDR by pyprophet.
#' @param runType (string) This must be one of the strings "DIA_proteomics", "DIA_Metabolomics".
#' @return (data-frames) A data-frame having following columns:
#' \item{transition_group_id}{(integer) a unique id for each precursor.}
#' \item{transition_id}{(list) fragment-ion ID associated with transition_group_id. This is matched with chromatogram ID in mzML file.}
#' \item{peptide_id}{(integer) a unique id for each peptide. A peptide can have multiple precursors.}
#' \item{sequence}{(string) amino-acid sequence of the precursor with possible modifications.}
#' \item{charge}{(integer) charge on the precursor.}
#' \item{group_label}{(string) TODO Figure it out.}
#'
#' @seealso \code{\link{getRunNames}, \link{fetchPrecursorsInfo}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' fileInfo <- getRunNames(dataPath = dataPath, oswMerged = TRUE)
#' precursors <- getPrecursorByID(c(32L, 2474L), fileInfo, oswMerged = TRUE)
#' @export
getPrecursorByID <- function(analytes, fileInfo, oswMerged = TRUE, runType = "DIA_proteomics"){
  if(oswMerged == TRUE){
    # Get precursor information from merged.osw file
    oswName <- unique(fileInfo[["featureFile"]])
    precursors <- fetchPrecursorsInfo(oswName, runType, analytes)
  } else {
    # Iterate over each file and collect precursor information
    precursors <- data.frame("transition_group_id" = integer())
    for(i in 1:nrow(fileInfo)){
      oswName <- fileInfo[["featureFile"]][[i]]
      temp <- fetchPrecursorsInfo(oswName, runType, analytes)
      precursors <- merge(precursors, temp, by = c("transition_group_id", all.x = TRUE, all.y = TRUE))
    }
  }
  message(paste0(nrow(precursors), " precursors are found."))
  precursors
}


#' Get features from a feature file.
#'
#' Get a data-frame of OpenSwath features that contains retention time, intensities, boundaries etc.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-04-04
#' @importFrom dplyr %>%
#' @param filename (string) Path to the feature file.
#' @param runID (string) id in RUN.ID column of the feature file.
#' @param maxFdrQuery (numeric) A numeric value between 0 and 1. It is used to filter features from osw file which have SCORE_MS2.QVALUE less than itself.
#' @param runType (char) This must be one of the strings "DIA_proteomics", "DIA_Metabolomics".
#' @return (data-frames) Data-frame has following columns:
#' \item{transition_group_id}{(integer) a unique id for each precursor.}
#' \item{RT}{(numeric) retention time as in FEATURE.EXP_RT of osw files.}
#' \item{Intensity}{(numeric) peak intensity as in FEATURE_MS2.AREA_INTENSITY of osw files.}
#' \item{leftWidth}{(numeric) as in FEATURE.LEFT_WIDTH of osw files.}
#' \item{rightWidth}{(numeric) as in FEATURE.RIGHT_WIDTH of osw files.}
#' \item{peak_group_rank}{(integer) rank of each feature associated with transition_group_id.}
#' \item{m_score}{(numeric) q-value of each feature associated with transition_group_id.}
#'
#' @seealso \code{\link{getRunNames}, \link{getFeatures}, \link{getFeaturesQuery}}
#' @keywords internal
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' fileInfo <- DIAlignR::getRunNames(dataPath = dataPath)
#' \dontrun{
#' featuresInfo <- fetchFeaturesFromRun(fileInfo$featureFile[1], fileInfo$spectraFileID[1],
#'  maxFdrQuery = 0.05)
#' dim(featuresInfo) # 211  7
#' }
fetchFeaturesFromRun <- function(filename, runID, maxFdrQuery = 1.00, runType = "DIA_proteomics"){
  # Establish a connection of SQLite file.
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = as.character(filename))
  # Generate a query.
  query <- getFeaturesQuery(runType)
  # Run query to get peptides, their coordinates and scores.
  featuresInfo <- tryCatch(expr = {output <- DBI::dbSendQuery(con, statement = query)
                          DBI::dbBind(output, list("FDR"=maxFdrQuery, "runID" = runID))
                          DBI::dbFetch(output)},
                          finally = {DBI::dbClearResult(output)
                            DBI::dbDisconnect(con)})
  featuresInfo
}


#' Get features from all feature files
#'
#' Get a list of data-frame of OpenSwath features that contains retention time, intensities, boundaries etc.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-04-06
#' @importFrom dplyr %>%
#' @param fileInfo (data-frame) Output of DIAlignR::getRunNames function.
#' @param maxFdrQuery (numeric) A numeric value between 0 and 1. It is used to filter features from osw file which have SCORE_MS2.QVALUE less than itself.
#' @param runType (char) This must be one of the strings "DIA_proteomics", "DIA_Metabolomics".
#' @return (data-frames) Data-frame has following columns:
#' \item{transition_group_id}{(integer) a unique id for each precursor.}
#' \item{RT}{(numeric) retention time as in FEATURE.EXP_RT of osw files.}
#' \item{Intensity}{(numeric) peak intensity as in FEATURE_MS2.AREA_INTENSITY of osw files.}
#' \item{leftWidth}{(numeric) as in FEATURE.LEFT_WIDTH of osw files.}
#' \item{rightWidth}{(numeric) as in FEATURE.RIGHT_WIDTH of osw files.}
#' \item{peak_group_rank}{(integer) rank of each feature associated with transition_group_id.}
#' \item{m_score}{(numeric) q-value of each feature associated with transition_group_id.}
#'
#' @seealso \code{\link{getRunNames}, \link{fetchPrecursorsInfo}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' fileInfo <- DIAlignR::getRunNames(dataPath = dataPath)
#' \dontrun{
#' features <- getFeatures(fileInfo, maxFdrQuery = 1.00, runType = "DIA_proteomics")
#' dim(features[[2]]) # 227  7
#' }
#' @export
getFeatures <- function(fileInfo, maxFdrQuery = 0.05, runType = "DIA_proteomics"){
  features <- vector(mode = "list", length = nrow(fileInfo))
  for(i in 1:nrow(fileInfo)){
    run <- rownames(fileInfo)[i]
    oswName <- fileInfo[["featureFile"]][[i]]
    runID <- fileInfo[["spectraFileID"]][[i]]
    names(runID) <- rownames(fileInfo)[[i]]
    features[[i]] <- fetchFeaturesFromRun(oswName, runID, maxFdrQuery, runType)
    message(paste0(nrow(features[[i]]), " peakgroups are founds below ", maxFdrQuery,
                   " FDR in run ", fileInfo[["runName"]][[i]], ", ID = ", runID))
  }
  names(features) <- rownames(fileInfo)
  features
}
