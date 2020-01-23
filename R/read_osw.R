#' Fetch features of analytes
#'
#' Get a data-frame of analytes' transition_group_ids, their OpenSwath features, chromatogram indices and associated FDR-scores.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + MIT
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
#' License: (c) Author (2019) + MIT
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
#' @seealso \code{\link{getRunNames}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' filenames <- DIAlignR::getRunNames(dataPath = dataPath)
#' oswFiles <- getOswAnalytes(dataPath = dataPath, filenames = filenames,
#'  analyteInGroupLabel = TRUE)
#' oswFiles[["run0"]][1,]
#' oswFiles <- getOswAnalytes(dataPath = dataPath, filenames = filenames,
#'  analyteInGroupLabel = FALSE)
#' oswFiles[["run0"]][1,]
#' @export
getOswAnalytes <- function(dataPath, filenames, oswMerged = TRUE, analyteInGroupLabel = FALSE,
                           maxFdrQuery = 0.05, runType  = "DIA_proteomics"){
  oswFiles <- list()
  for(i in 1:nrow(filenames)){
    # Get a query to search against the osw files.
    if(oswMerged == TRUE){
      oswName <- list.files(path = file.path(dataPath, "osw"), pattern="*merged.osw")
      oswName <- file.path(dataPath, "osw", oswName[1])
    } else{
      oswName <- paste0(file.path(dataPath, "osw", filenames$runs[i]), ".osw")
    }
    # Establish a connection of SQLite file.
    con <- DBI::dbConnect(RSQLite::SQLite(), dbname = oswName)
    # Generate a query.
    query <- getAnalytesQuery(maxFdrQuery = maxFdrQuery, oswMerged = oswMerged,
                      filename = filenames$filename[i], runType = runType,
                      analyteInGroupLabel = analyteInGroupLabel)
    # Run query to get peptides, their coordinates and scores.
    oswAnalytes <- tryCatch(expr = DBI::dbGetQuery(con, statement = query),
                             finally = DBI::dbDisconnect(con))
    oswFiles[[i]] <- oswAnalytes
  }
  names(oswFiles) <- rownames(filenames)
  oswFiles
}

