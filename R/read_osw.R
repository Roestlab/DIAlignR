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
#' filenames <- getRunNames(dataPath = dataPath)
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
#' To get filenames use \code{\link{getRunNames}} function.
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
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @inheritParams getPrecursors
#' @param filename (string) Should be from the RUN.FILENAME column from osw files.
#' @param runType (string) This must be one of the strings "DIA_proteomics", "DIA_Metabolomics".
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
#' precursorsInfo <- fetchPrecursorsInfo(filename, runType = "DIA_proteomics", context = "experiment-wide")
#' dim(precursorsInfo) # 234  6
#' }
fetchPrecursorsInfo <- function(filename, runType, selectIDs = NULL,
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
                                        if(all) {DBI::dbBind(output, list("CONTEXT"=context, "FDR"=maxPeptideFdr))}
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
#' @importFrom magrittr %>%
#' @param fileInfo (data-frame) Output of \code{\link{getRunNames}} function.
#' @param oswMerged (logical) TRUE for experiment-wide FDR and FALSE for run-specific FDR by pyprophet.
#' @param runType (string) This must be one of the strings "DIA_proteomics", "DIA_Metabolomics".
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
#' precursorsInfo <- getPrecursors(fileInfo, oswMerged = TRUE, runType = "DIA_proteomics",
#' context = "experiment-wide", maxPeptideFdr = 0.05)
#' dim(precursorsInfo) # 234  6
#' @export
getPrecursors <- function(fileInfo, oswMerged = TRUE, runType = "DIA_proteomics",
                          context = "global", maxPeptideFdr = 0.05, level = "Peptide"){
  if(oswMerged == TRUE){
    # Get precursor information from merged.osw file
    oswName <- unique(fileInfo[["featureFile"]])
    precursors <- fetchPrecursorsInfo(oswName, runType, NULL, context, maxPeptideFdr, level)
  } else {
    # Iterate over each file and collect precursor information
    oswName <- fileInfo[["featureFile"]][[1]]
    precursors <- fetchPrecursorsInfo(oswName, runType, NULL, context, maxPeptideFdr, level)
    #precursors <- data.frame()
    #ids <- integer()
    #for(i in 1:nrow(fileInfo)){
    #  oswName <- fileInfo[["featureFile"]][[i]]
    #  temp <- fetchPrecursorsInfo(oswName, runType, NULL, context, maxPeptideFdr, level)
    #  idx <- !(temp$transition_group_id %in% ids)
    #  precursors <- rbind(precursors, temp[idx,])
    #  ids <- precursors$transition_group_id
    # }
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
#' @importFrom magrittr %>%
#' @import dplyr
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
getPrecursorByID <- function(analytes, fileInfo, oswMerged = TRUE, runType = "DIA_proteomics"){
  if(oswMerged == TRUE){
    # Get precursor information from merged.osw file
    oswName <- unique(fileInfo[["featureFile"]])
    precursors <- fetchPrecursorsInfo(oswName, runType, analytes, maxPeptideFdr = 1.00)
  } else {
    # Iterate over each file and collect precursor information
    precursors <- data.frame("transition_group_id" = integer())
    for(i in 1:nrow(fileInfo)){
      oswName <- fileInfo[["featureFile"]][[i]]
      temp <- fetchPrecursorsInfo(oswName, runType, analytes, maxPeptideFdr = 1.00)
      idx <- !(temp$transition_group_id %in% precursors$transition_group_id)
      precursors <- dplyr::bind_rows(precursors, temp[idx,])
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
#' @importFrom magrittr %>%
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
#' fileInfo <- getRunNames(dataPath = dataPath)
#' \dontrun{
#' featuresInfo <- fetchFeaturesFromRun(fileInfo$featureFile[1], fileInfo$spectraFileID[1],
#'  maxFdrQuery = 0.05)
#' dim(featuresInfo) # 211  8
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
#' @importFrom magrittr %>%
#' @inheritParams alignTargetedRuns
#' @param fileInfo (data-frame) output of \code{\link{getRunNames}} function.
#' @param maxFdrQuery (numeric) a numeric value between 0 and 1. It is used to filter features from osw file which have SCORE_MS2.QVALUE less than itself.
#' @param runType (char) yhis must be one of the strings "DIA_proteomics", "DIA_Metabolomics".
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
#' features <- getFeatures(fileInfo, maxFdrQuery = 1.00, runType = "DIA_proteomics")
#' dim(features[[2]]) # 938  8
#' @export
getFeatures <- function(fileInfo, maxFdrQuery = 0.05, runType = "DIA_proteomics", applyFun = lapply){
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
#' precursorsInfo <- fetchPeptidesInfo(fileInfo, runType = "DIA_proteomics", context = "experiment-wide")
#' }
fetchPeptidesInfo <- function(oswName, runType, context){
  # Establish a connection of SQLite file.
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = oswName)
  # Generate a query.
  query <- getPeptideQuery(runType)

  # Run query to get peptides, their scores and pvalues.
  peptidesInfo <- tryCatch(expr = {output <- DBI::dbSendQuery(con, statement = query)
                     DBI::dbBind(output, list("CONTEXT"=context))
                     DBI::dbFetch(output)},
              finally = {DBI::dbClearResult(output)
                         DBI::dbDisconnect(con)})

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
#' precursorsInfo <- fetchPeptidesInfo(fileInfo, runType = "DIA_proteomics", context = "experiment-wide")
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
#' @import dplyr bit64
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
#' precursorsInfo <- getPrecursors(fileInfo, oswMerged = TRUE, runType = "DIA_proteomics",
#' context = "experiment-wide", maxPeptideFdr = 0.05)
#' peptidesInfo <- getPeptideScores(fileInfo, unique(precursorsInfo$peptide_id))
#' dim(peptidesInfo) # 684 5
#' @export
getPeptideScores <- function(fileInfo, peptides, oswMerged = TRUE, runType = "DIA_proteomics", context = "global"){
  if(context == "global") context <- "experiment-wide"

  if(oswMerged == TRUE){
    # Get precursor information from merged.osw file
    oswName <- unique(fileInfo[["featureFile"]])
    peptidesInfo <- fetchPeptidesInfo(oswName, runType, context)
  } else {
    oswName <- fileInfo[["featureFile"]][[1]]
    peptidesInfo <- fetchPeptidesInfo(oswName, runType, context)
    # Iterate over each file and collect peptides information
    #peptidesInfo <- data.frame()
    #for(i in 1:nrow(fileInfo)){
    #  oswName <- fileInfo[["featureFile"]][[i]]
    #  runID <- fileInfo[["spectraFileID"]][[i]]
    #  names(runID) <- rownames(fileInfo)[[i]]
    #  temp <- fetchPeptidesInfo2(oswName, runType, context, runID)
    #  peptidesInfo <- dplyr::bind_rows(peptidesInfo, temp)
    #}
  }
  peptidesInfo <- dplyr::filter(peptidesInfo, .data$peptide_id %in% peptides)
  runs <- bit64::as.integer64(fileInfo$spectraFileID)
  peptidesInfo <- dplyr::filter(peptidesInfo, .data$run %in% runs)
  peptidesInfo$run <- rownames(fileInfo)[match(peptidesInfo$run, runs)]
  if(length(unique(peptidesInfo$peptide_id)) != length(peptides)) {
    warning("Unable to find scores for few peptides. Appending NAs.")
    temp <- data.frame(peptide_id = setdiff(peptides, unique(peptidesInfo$peptide_id)),
                       run = NA_character_, score = NA_real_, pvalue = NA_real_, qvalue = NA_real_)
    peptidesInfo <- rbind(peptidesInfo, temp)
  }
  message(nrow(peptidesInfo), " peptides scores are fetched.")
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
#' @importFrom magrittr %>%
#' @param filename (string) Path to the feature file.
#' @param runID (string) id in RUN.ID column of the feature file.
#' @param maxFdrQuery (numeric) A numeric value between 0 and 1. It is used to filter features from osw file which have SCORE_MS2.QVALUE less than itself.
#' @param runType (char) This must be one of the strings "DIA_proteomics", "DIA_Metabolomics".
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
fetchTransitionsFromRun <- function(filename, runID, maxFdrQuery = 1.00, runType = "DIA_proteomics"){
  # Establish a connection of SQLite file.
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = as.character(filename))
  # Generate a query.
  query <- getTransitionsQuery(runType)
  # Run query to get peptides, their coordinates and scores.
  transitionInfo <- tryCatch(expr = {output <- DBI::dbSendQuery(con, statement = query)
  DBI::dbBind(output, list("FDR"=maxFdrQuery, "runID" = runID))
  DBI::dbFetch(output)},
  finally = {DBI::dbClearResult(output)
    DBI::dbDisconnect(con)})
  df1 <- dplyr::select(transitionInfo, -.data$intensity) %>%
    dplyr::group_by(.data$transition_group_id, .data$peak_group_rank) %>%
    filter(row_number()==1) %>% dplyr::ungroup() %>% as.data.frame()
  df2 <- dplyr::group_by(transitionInfo, .data$transition_group_id, .data$peak_group_rank) %>%
    dplyr::summarise(intensity = base::list(.data$intensity)) %>% dplyr::ungroup() %>% as.data.frame()
  transitionInfo <- merge(df1, df2, by = c("transition_group_id", "peak_group_rank"))
  transitionInfo <- dplyr::select(transitionInfo, .data$transition_group_id, .data$feature_id, .data$RT,
                                  .data$intensity, .data$leftWidth, .data$rightWidth, .data$peak_group_rank, .data$m_score) %>%
    dplyr::arrange(transitionInfo, .data$transition_group_id, .data$peak_group_rank)
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
#' @importFrom magrittr %>%
#' @inheritParams alignTargetedRuns
#' @param fileInfo (data-frame) output of \code{\link{getRunNames}} function.
#' @param maxFdrQuery (numeric) a numeric value between 0 and 1. It is used to filter features from osw file which have SCORE_MS2.QVALUE less than itself.
#' @param runType (char) yhis must be one of the strings "DIA_proteomics", "DIA_Metabolomics".
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
#' transitions <- getTransitions(fileInfo, maxFdrQuery = 1.00, runType = "DIA_proteomics")
#' dim(transitions[[2]]) # 938  8
#' @export
getTransitions <- function(fileInfo, maxFdrQuery = 0.05, runType = "DIA_proteomics", applyFun = lapply){
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
