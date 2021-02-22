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
#' @param runType (char) This must be one of the strings "DIA_Proteomics", "DIA_Metabolomics".
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
#' filenames <- getRunNames(dataPath = dataPath)
#' oswName <- paste0(dataPath,"/osw/merged.osw")
#' \dontrun{
#' analytesInfo <- fetchAnalytesInfo(oswName, maxFdrQuery = 0.05, oswMerged = TRUE,
#'  analytes = c("19051_KLIVTSEGC[160]FK/2"), filename = filenames$filename[2],
#'   runType = "DIA_Proteomics", analyteInGroupLabel = TRUE)
#' analytesInfo <- fetchAnalytesInfo(oswName, maxFdrQuery = 0.05, oswMerged = TRUE,
#'  analytes = c("IHFLSPVRPFTLTPGDEEESFIQLITPVR_3"), filename = filenames$filename[3],
#'   runType = "DIA_Proteomics", analyteInGroupLabel = FALSE)
#' }
fetchAnalytesInfo <- function(oswName, maxFdrQuery, oswMerged,
                              analytes, filename, runType = "DIA_Proteomics", analyteInGroupLabel = FALSE){
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
#' @param dataPath (char) path to xics and osw directory.
#' @param filenames (data-frame) column "filename" contains RUN table from osw files. column "runs" contain respective mzML names without extension.
#' To get filenames use \code{\link{getRunNames}} function.
#' @param oswMerged (logical) TRUE for experiment-wide FDR and FALSE for run-specific FDR by pyprophet.
#' @param analyteInGroupLabel (logical) TRUE for getting analytes as PRECURSOR.GROUP_LABEL from osw file.
#'  FALSE for fetching analytes as PEPTIDE.MODIFIED_SEQUENCE and PRECURSOR.CHARGE from osw file.
#' @param maxFdrQuery (numeric) A numeric value between 0 and 1. It is used to filter features from osw file which have SCORE_MS2.QVALUE less than itself.
#' @param runType (char) This must be one of the strings "DIA_Proteomics", "DIA_Metabolomics".
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
                           maxFdrQuery = 0.05, runType = "DIA_Proteomics"){
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
#' @importFrom data.table setDT
#' @inheritParams getPrecursors
#' @param filename (string) Should be from the RUN.FILENAME column from osw files.
#' @param runType (string) This must be one of the strings "DIA_Proteomics", "DIA_Metabolomics".
#' @param selectIDs (integer) a vector of integers.
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
#' precursorsInfo <- fetchPrecursorsInfo(filename, runType = "DIA_Proteomics", context = "experiment-wide")
#' dim(precursorsInfo) # 234  6
#' }
fetchPrecursorsInfo <- function(filename, runType = "DIA_Proteomics", selectIDs = NULL,
                                context = "global", maxPeptideFdr = 0.05, level = "Peptide"){
  # Establish a connection of SQLite file.
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = as.character(filename))
  # Generate a query.
  all = FALSE
  if(is.null(selectIDs)) all = TRUE
  if(all){
    query <- getPrecursorsQuery(runType, level)
  } else{
    query <- getPrecursorsQueryID(selectIDs, runType)
  }

  # Run query to get peptides, their coordinates and scores.
  precursorsInfo <- tryCatch(expr = { output <- DBI::dbSendQuery(con, statement = query)
                     if(all && runType != "DIA_Metabolomics") {DBI::dbBind(output, list("CONTEXT"=context, "FDR"=maxPeptideFdr))}
                                        DBI::dbFetch(output)},
                             finally = {DBI::dbClearResult(output)
                             DBI::dbDisconnect(con)})
  # Each precursor has only one row.
  setDT(precursorsInfo)
  precursorsInfo[, `:=`(transition_ids = list(transition_id)),
                 by = .(transition_group_id)][, transition_id := NULL]
  precursorsInfo <- unique(precursorsInfo, by = c("transition_group_id"))
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
#' @importFrom data.table data.table setkeyv
#' @param fileInfo (data-frame) Output of \code{\link{getRunNames}} function.
#' @param oswMerged (logical) TRUE for experiment-wide FDR and FALSE for run-specific FDR by pyprophet.
#' @param runType (char) This must be one of the strings "DIA_Proteomics", "DIA_Metabolomics".
#' @param context (string) Context used in pyprophet peptide. Must be either "run-specific", "experiment-wide", or "global".
#' @param maxPeptideFdr (numeric) A numeric value between 0 and 1. It is used to filter peptides from osw file which have SCORE_PEPTIDE.QVALUE less than itself.
#' @param level (string) Apply maxPeptideFDR on Protein as well if specified as "Protein". Default: "Peptide".
#' @return (data-frames) A data-frame having following columns:
#' \item{transition_group_id}{(integer) a unique id for each precursor.}
#' \item{peptide_id}{(integer) a unique id for each peptide. A peptide can have multiple precursors.}
#' \item{sequence}{(string) amino-acid sequence of the precursor with possible modifications.}
#' \item{charge}{(integer) charge on the precursor.}
#' \item{group_label}{(string) TODO Figure it out.}
#' \item{transition_ids}{(list) fragment-ion ID associated with transition_group_id. This is matched with chromatogram ID in mzML file.}
#'
#' @seealso \code{\link{getRunNames}, \link{fetchPrecursorsInfo}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' fileInfo <- getRunNames(dataPath = dataPath)
#' precursorsInfo <- getPrecursors(fileInfo, oswMerged = TRUE, runType = "DIA_Proteomics",
#' context = "experiment-wide", maxPeptideFdr = 0.05)
#' dim(precursorsInfo) # 234  6
#' @export
getPrecursors <- function(fileInfo, oswMerged = TRUE, runType = "DIA_Proteomics",
                          context = "global", maxPeptideFdr = 0.05, level = "Peptide"){
  if(oswMerged == TRUE){
    # Get precursor information from merged.osw file
    oswName <- unique(fileInfo[["featureFile"]])
    precursors <- fetchPrecursorsInfo(oswName, runType, NULL, context, maxPeptideFdr, level)
  } else {
    # Iterate over each file and collect precursor information
    oswName <- fileInfo[["featureFile"]][[1]]
    precursors <- fetchPrecursorsInfo(oswName, runType, NULL, context, maxPeptideFdr, level)
    #for(i in 1:nrow(fileInfo)){
      #oswName <- fileInfo[["featureFile"]][[i]]
      #temp <- fetchPrecursorsInfo(oswName, runType, NULL, context, maxPeptideFdr, level)
      #temp <- temp[!precursors, on = .(transition_group_id)] # Anti-join
      #precursors <- rbind(precursors, temp, use.names = FALSE)
    #}
  }
  setkeyv(precursors, c("peptide_id", "transition_group_id"))
  message(precursors[,.N], " precursors are found.")
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
#' @importFrom data.table .N
#' @inheritParams getPrecursors
#' @param analytes (integer) a vector of integers.
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
getPrecursorByID <- function(analytes, fileInfo, oswMerged = TRUE, runType = "DIA_Proteomics"){
  if(oswMerged == TRUE){
    # Get precursor information from merged.osw file
    oswName <- unique(fileInfo[["featureFile"]])
    precursors <- fetchPrecursorsInfo(oswName, runType, analytes, maxPeptideFdr = 1.00)
  } else {
    # Iterate over each file and collect precursor information
    precursors <- data.table("transition_group_id" = integer(), "transition_id"= integer(), "peptide_id" = integer(),
                             "sequence" = character(), "charge" = integer(), "group_label" = character(),
                             "transition_ids" = list())
    for(i in 1:nrow(fileInfo)){
      oswName <- fileInfo[["featureFile"]][[i]]
      temp <- fetchPrecursorsInfo2(oswName, runType, analytes, maxPeptideFdr = 1.00)
      temp <- temp[!precursors, on = .(transition_group_id)] # Anti-join
      precursors <- rbind(precursors, temp, use.names = FALSE)
    }
  }
  setkeyv(precursors, c("peptide_id", "transition_group_id"))
  message(precursors[,.N], " precursors are found.")
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
#' @importFrom data.table setDT setkey
#' @param filename (string) Path to the feature file.
#' @param runID (string) id in RUN.ID column of the feature file.
#' @param maxFdrQuery (numeric) A numeric value between 0 and 1. It is used to filter features from osw file which have SCORE_MS2.QVALUE less than itself.
#' @param runType (char) This must be one of the strings "DIA_Proteomics", "DIA_Metabolomics".
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
#' fileInfo <- getRunNames(dataPath = dataPath)
#' \dontrun{
#' featuresInfo <- fetchFeaturesFromRun(fileInfo$featureFile[1], fileInfo$spectraFileID[1],
#'  maxFdrQuery = 0.05)
#' dim(featuresInfo) # 211  8
#' }
fetchFeaturesFromRun <- function(filename, runID, maxFdrQuery = 1.00, runType = "DIA_Proteomics"){
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
  setDT(featuresInfo)
  setkey(featuresInfo, "transition_group_id")
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
#' @inheritParams alignTargetedRuns
#' @param fileInfo (data-frame) output of \code{\link{getRunNames}} function.
#' @param maxFdrQuery (numeric) a numeric value between 0 and 1. It is used to filter features from osw file which have SCORE_MS2.QVALUE less than itself.
#' @param runType (char) yhis must be one of the strings "DIA_Proteomics", "DIA_Metabolomics".
#' @return (list of dataframes) each dataframe has following columns:
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
#' fileInfo <- getRunNames(dataPath = dataPath)
#' \dontrun{
#' features <- getFeatures(fileInfo, maxFdrQuery = 1.00, runType = "DIA_Proteomics")
#' dim(features[[2]]) # 938  8
#' }
#' @export
getFeatures <- function(fileInfo, maxFdrQuery = 0.05, runType = "DIA_Proteomics", applyFun = lapply){
  features <- applyFun(1:nrow(fileInfo), function(i){
    run <- rownames(fileInfo)[i]
    oswName <- fileInfo[["featureFile"]][[i]]
    runID <- fileInfo[["spectraFileID"]][[i]]
    names(runID) <- rownames(fileInfo)[[i]]
    df <- fetchFeaturesFromRun(oswName, runID, maxFdrQuery, runType)
    message(paste0(nrow(df), " peakgroups are founds below ", maxFdrQuery,
                   " FDR in run ", fileInfo[["runName"]][[i]], ", ID = ", runID))
    df
  })
  names(features) <- rownames(fileInfo)
  features
}

dummyFeatures <- function(precursors, numMerge = 0L, startIdx = 1L){
  stpIdx <- startIdx + numMerge - 1
  masters <- paste0("master", startIdx:stpIdx)
  transition_group_ids <- .subset2(precursors, "transition_group_id")

  features <- lapply(masters, function(run) {
    data.table("transition_group_id" = rep(transition_group_ids, each = 5L),
               "feature_id" = bit64::NA_integer64_, "RT" = NA_real_, "intensity" = NA_real_,
               "leftWidth" = NA_real_, "rightWidth" = NA_real_, "peak_group_rank" = NA_integer_,
               "m_score" = NA_real_, key = "transition_group_id")
  })
  names(features) <- masters
  features
}

#' Get scores of all peptides
#'
#' Return a scores, pvalues, and qvalues for all peptides from the osw file.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-07-01
#' @keywords internal
#' @importFrom data.table setDT
#' @inheritParams getPrecursors
#' @param oswName (char) path to the osw file.
#' @return (dataframe) with following columns:
#' \item{peptide_id}{(integer) a unique id for each precursor.}
#' \item{run}{(character) as in SCORE_PEPTIDE.RUN_ID of osw files.}
#' \item{score}{(numeric) as in SCORE_PEPTIDE.SCORE of osw files.}
#' \item{pvalue}{(numeric) as in SCORE_PEPTIDE.PVALUE of osw files.}
#' \item{qvalue}{(numeric) as in SCORE_PEPTIDE.QVALUE of osw files.}
#'
#' @seealso \code{\link{getPeptideQuery}, \link{getPeptideScores}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' fileInfo <- getRunNames(dataPath = dataPath)
#' oswName <- fileInfo[["featureFile"]][1]
#' \dontrun{
#' precursorsInfo <- fetchPeptidesInfo(fileInfo, runType = "DIA_Proteomics", context = "experiment-wide")
#' }
fetchPeptidesInfo <- function(oswName, runType, context){
  # Establish a connection of SQLite file.
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = oswName)
  # Generate a query.
  query <- getPeptideQuery(runType)

  # Run query to get peptides, their scores and pvalues.
  peptidesInfo <- tryCatch(expr = {output <- DBI::dbSendQuery(con, statement = query)
                  if(runType != "DIA_Metabolomics") DBI::dbBind(output, list("CONTEXT"=context))
                     DBI::dbFetch(output)},
              finally = {DBI::dbClearResult(output)
                         DBI::dbDisconnect(con)})

  peptidesInfo <- peptidesInfo[,c("peptide_id", "run", "score", "pvalue", "qvalue")]
  setDT(peptidesInfo)
  peptidesInfo
}

#' Get scores of all peptides
#'
#' Return a scores, pvalues, and qvalues for all peptides from the osw file.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-11-18
#' @keywords internal
#' @inheritParams getPrecursors
#' @param oswName (char) path to the osw file.
#' @return (dataframe) with following columns:
#' \item{peptide_id}{(integer) a unique id for each precursor.}
#' \item{run}{(character) as in SCORE_PEPTIDE.RUN_ID of osw files.}
#' \item{score}{(numeric) as in SCORE_PEPTIDE.SCORE of osw files.}
#' \item{pvalue}{(numeric) as in SCORE_PEPTIDE.PVALUE of osw files.}
#' \item{qvalue}{(numeric) as in SCORE_PEPTIDE.QVALUE of osw files.}
#'
#' @seealso \code{\link{getPeptideQuery}, \link{getPeptideScores}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' fileInfo <- getRunNames(dataPath = dataPath)
#' oswName <- fileInfo[["featureFile"]][1]
#' \dontrun{
#' precursorsInfo <- fetchPeptidesInfo(fileInfo, runType = "DIA_Proteomics", context = "experiment-wide")
#' }
fetchPeptidesInfo2 <- function(oswName, runType, context, runID){
  # Establish a connection of SQLite file.
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = oswName)
  # Generate a query.
  query <- getPeptideQuery2(runType)

  # Run query to get peptides, their scores and pvalues.
  peptidesInfo <- tryCatch(expr = {output <- DBI::dbSendQuery(con, statement = query)
  DBI::dbBind(output, list("CONTEXT"=context, "runID" = runID ))
  DBI::dbFetch(output)},
  finally = {DBI::dbClearResult(output)
    DBI::dbDisconnect(con)})

  setDT(peptidesInfo)
  peptidesInfo
}


#' Get scores of peptide
#'
#' Get a list of dataframes that contains peptide scores, pvalues, and qvalues across all runs.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-07-01
#' @importFrom data.table setnames setcolorder setkey
#' @inheritParams getPrecursors
#' @param peptides (integer) Ids of peptides for which scores are required.
#' @return (list of dataframes) dataframe has following columns:
#' \item{peptide_id}{(integer) a unique id for each precursor.}
#' \item{run}{(character) as in SCORE_PEPTIDE.RUN_ID of osw files.}
#' \item{score}{(numeric) as in SCORE_PEPTIDE.SCORE of osw files.}
#' \item{pvalue}{(numeric) as in SCORE_PEPTIDE.PVALUE of osw files.}
#' \item{qvalue}{(numeric) as in SCORE_PEPTIDE.QVALUE of osw files.}
#'
#' @seealso \code{\link{getRunNames}, \link{fetchPeptidesInfo}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' fileInfo <- getRunNames(dataPath = dataPath)
#' precursorsInfo <- getPrecursors(fileInfo, oswMerged = TRUE, runType = "DIA_Proteomics",
#' context = "experiment-wide", maxPeptideFdr = 0.05)
#' peptidesInfo <- getPeptideScores(fileInfo, unique(precursorsInfo$peptide_id))
#' dim(peptidesInfo) # 684 5
#' @export
getPeptideScores <- function(fileInfo, peptides, oswMerged = TRUE, runType = "DIA_Proteomics", context = "global"){
  if(context == "global") context <- "experiment-wide"

  if(oswMerged == TRUE){
    # Get precursor information from merged.osw file
    oswName <- unique(fileInfo[["featureFile"]])
    peptidesInfo <- fetchPeptidesInfo(oswName, runType, context)
  } else {
    oswName <- fileInfo[["featureFile"]][[1]]
    peptidesInfo <- fetchPeptidesInfo(oswName, runType, context)
  }

  ids <- bit64::as.integer64(fileInfo$spectraFileID)
  peptidesInfo <- peptidesInfo[list(ids), on = "run"]
  peptidesInfo[, col2 := rownames(fileInfo)[match(run, ids)]][,run:= NULL]
  setnames(peptidesInfo, "col2", "run")
  setcolorder(peptidesInfo, c("peptide_id","run"))
  peptidesInfo <- peptidesInfo[list(peptides), on = "peptide_id"]

  setkey(peptidesInfo, peptide_id)
  message(peptidesInfo[,.N], " peptides scores are fetched.")
  peptidesInfo
}


#' Get transitions from a feature file
#'
#' Get a data-frame of OpenSwath features that contains retention time, transition intensities, boundaries etc.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-11-15
#' @importFrom data.table setnames setDT setcolorder .SD ":="
#' @param filename (string) Path to the feature file.
#' @param runID (string) id in RUN.ID column of the feature file.
#' @param maxFdrQuery (numeric) A numeric value between 0 and 1. It is used to filter features from osw file which have SCORE_MS2.QVALUE less than itself.
#' @param runType (char) This must be one of the strings "DIA_Proteomics", "DIA_Metabolomics".
#' @return (data-frames) Data-frame has following columns:
#' \item{transition_group_id}{(integer) a unique id for each precursor.}
#' \item{RT}{(numeric) retention time as in FEATURE.EXP_RT of osw files.}
#' \item{intensity}{(list) of peak intensities as in FEATURE_TRANSITION.AREA_INTENSITY of osw files.}
#' \item{leftWidth}{(numeric) as in FEATURE.LEFT_WIDTH of osw files.}
#' \item{rightWidth}{(numeric) as in FEATURE.RIGHT_WIDTH of osw files.}
#' \item{peak_group_rank}{(integer) rank of each feature associated with transition_group_id.}
#' \item{m_score}{(numeric) q-value of each feature associated with transition_group_id.}
#'
#' @seealso \code{\link{getRunNames}, \link{getTransitions}, \link{getTransitionsQuery}}
#' @keywords internal
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' fileInfo <- getRunNames(dataPath = dataPath)
#' \dontrun{
#' transitionsInfo <- fetchTransitionsFromRun(fileInfo$featureFile[1], fileInfo$spectraFileID[1],
#'  maxFdrQuery = 0.05)
#' dim(transitionsInfo) # 211  8
#' }
fetchTransitionsFromRun <- function(filename, runID, maxFdrQuery = 1.00, runType = "DIA_Proteomics"){
  # Establish a connection of SQLite file.
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = filename)
  # Generate a query.
  query <- getTransitionsQuery(runType)
  # Run query to get peptides, their coordinates and scores.
  transitionInfo <- tryCatch(expr = {output <- DBI::dbSendQuery(con, statement = query)
                  DBI::dbBind(output, list("FDR"=maxFdrQuery, "runID" = runID))
                  DBI::dbFetch(output)},
                  finally = {DBI::dbClearResult(output)
                    DBI::dbDisconnect(con)})
  setDT(transitionInfo)
  transitionInfo <- transitionInfo[, `:=`(intensity2 = list(intensity)),
                                   keyby = .(transition_group_id, peak_group_rank)][
       ,intensity := NULL][
       , head(.SD, 1), by=.(transition_group_id, peak_group_rank),
       .SDcols = c("feature_id", "RT", "intensity2", "leftWidth", "rightWidth", "m_score")]
  setnames(transitionInfo, "intensity2", "intensity")
  setcolorder(transitionInfo, c("transition_group_id", "feature_id", "RT", "intensity", "leftWidth",
                                "rightWidth", "peak_group_rank", "m_score"))
  transitionInfo
}


#' Get transitions from all feature files
#'
#' Get a list of data-frame of OpenSwath features that contains retention time, intensities, boundaries etc.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-11-15
#' @inheritParams alignTargetedRuns
#' @param fileInfo (data-frame) output of \code{\link{getRunNames}} function.
#' @param maxFdrQuery (numeric) a numeric value between 0 and 1. It is used to filter features from osw file which have SCORE_MS2.QVALUE less than itself.
#' @param runType (char) yhis must be one of the strings "DIA_Proteomics", "DIA_Metabolomics".
#' @return (list of dataframes) each dataframe has following columns:
#' \item{transition_group_id}{(integer) a unique id for each precursor.}
#' \item{RT}{(numeric) retention time as in FEATURE.EXP_RT of osw files.}
#' \item{intensity}{(list) of peak intensities as in FEATURE_TRANSITION.AREA_INTENSITY of osw files.}
#' \item{leftWidth}{(numeric) as in FEATURE.LEFT_WIDTH of osw files.}
#' \item{rightWidth}{(numeric) as in FEATURE.RIGHT_WIDTH of osw files.}
#' \item{peak_group_rank}{(integer) rank of each feature associated with transition_group_id.}
#' \item{m_score}{(numeric) q-value of each feature associated with transition_group_id.}
#'
#' @seealso \code{\link{getRunNames}, \link{fetchTransitionsFromRun}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' fileInfo <- getRunNames(dataPath = dataPath)
#' transitions <- getTransitions(fileInfo, maxFdrQuery = 1.00, runType = "DIA_Proteomics")
#' dim(transitions[[2]]) # 938  8
#' @export
getTransitions <- function(fileInfo, maxFdrQuery = 0.05, runType = "DIA_Proteomics", applyFun = lapply){
  transitions <- applyFun(1:nrow(fileInfo), function(i){
    run <- rownames(fileInfo)[i]
    oswName <- fileInfo[["featureFile"]][[i]]
    runID <- fileInfo[["spectraFileID"]][[i]]
    names(runID) <- rownames(fileInfo)[[i]]
    df <- fetchTransitionsFromRun(oswName, runID, maxFdrQuery, runType)
    message(paste0(nrow(df), " peakgroups are founds below ", maxFdrQuery,
                   " FDR in run ", fileInfo[["runName"]][[i]], ", ID = ", runID))
    df
  })
  names(transitions) <- rownames(fileInfo)
  transitions
}
