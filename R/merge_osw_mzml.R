#' Coerce chromatogram ids as integer
#'
#' chromatogramHeader has 10 columns. The two important columns are:
#' "chromatogramId" which has fragment-ion ID that matches with transition ID in osw file.
#' "chromatogramIndex" that lists indices of chromatograms in mzML file.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + MIT
#' Date: 2019-12-13
#' @param chromatogramHeader (dataframe)
#' @return Invisible NULL
chromatogramIdAsInteger <- function(chromatogramHeader){
  assign("chromHead", dplyr::mutate(dplyr::select(chromatogramHeader, chromatogramId, chromatogramIndex),
                                    chromatogramId = as.integer(chromatogramId)),
         envir = parent.frame(n = 1))
  invisible(NULL)
}


#' Merge dataframes from OSW and mzML files
#'
#' Merges dataframes on transition_id(OSW) = chromatogramId(mzML).
#' @importFrom dplyr %>%
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + MIT
#' Date: 2019-12-13
#'
#' @param oswAnalytes (dataframe) This is an output of getOswFiles.
#' @param chromHead (dataframe) This has two columns: chromatogramId and chromatogramIndex with integer values.
#' @param analyteFDR (numeric) Not used.
#' @param runType (char) This must be one of the strings "DIA_proteomics", "DIA_Metabolomics".
#' @return Invisible NULL
#' @seealso \code{\link{getOswFiles}}
mergeOswAnalytes_ChromHeader <- function(oswAnalytes, chromHead, analyteFDR =  1.00, runType = "DIA_proteomics"){
  # TODO: Make sure that transition_id has same order across runs. IMO should be specified in query.
  assign("oswAnalytes", dplyr::left_join(oswAnalytes, chromHead,
                                  by = c("transition_id" = "chromatogramId")) %>%
    dplyr::group_by(transition_group_id, peak_group_rank) %>%
    dplyr::mutate(transition_ids = paste0(transition_id, collapse = ","),
                  chromatogramIndex = paste0(chromatogramIndex, collapse = ",")) %>%
    dplyr::ungroup() %>% dplyr::select(-transition_id) %>% dplyr::distinct(),
    envir = parent.frame(n = 1))
  invisible(NULL)
}

#' Get list of peptides and their chromatogram indices.
#'
#' This function reads all osw and mzml files in the directories at dataPath. It selects analytes which has associated features with m-score < maxFdrQuery.
#' For these analytes it fetches chromatogram indices by matching transition_id(osw) with chromatogramID(mzml).
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + MIT
#' Date: 2019-12-13
#' @param dataPath (char) path to mzml and osw directory.
#' @param filenames (data-frame) column "filename" contains RUN table from osw files. column "runs" contain respective mzML names without extension.
#' To get filenames use DIAlignR::getRunNames function.
#' @param maxFdrQuery (numeric) A numeric value between 0 and 1. It is used to filter features from osw file which have SCORE_MS2.QVALUE less than itself.
#' @param analyteFDR (numeric) Not used.
#' @param oswMerged (logical) TRUE for experiment-wide FDR and FALSE for run-specific FDR by pyprophet.
#' @param analytes (string) analyte is as PRECURSOR.GROUP_LABEL or as PEPTIDE.MODIFIED_SEQUENCE and PRECURSOR.CHARGE from osw file.
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
#' \item{chromatogramIndex}{(integer) Index of chromatogram in mzML file.}
#' \item{transition_ids}{(integer) fragment-ion ID associated with transition_group_id. This is matched with chromatogram ID in mzML file.}
#'
#' @seealso \code{\link{getRunNames}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' filenames <- DIAlignR::getRunNames(dataPath = dataPath)
#' oswFiles <- getOswFiles(dataPath = dataPath, filenames =filenames, analyteInGroupLabel = TRUE)
#' @export
getOswFiles <- function(dataPath, filenames, maxFdrQuery = 0.05, analyteFDR = 0.01, oswMerged = TRUE,
                         analytes = NULL, runType = "DIA_proteomics", analyteInGroupLabel = FALSE){
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
    # Get transition indices for MS2 fragment-ions.
    oswAnalytes <- fetchAnalytesInfo(oswName, maxFdrQuery, oswMerged, analytes = analytes,
                                     filename = filenames$filename[i], runType, analyteInGroupLabel)

    # Get chromatogram indices from the header file.
    mzmlName <- file.path(dataPath, "mzml", paste0(filenames$runs[i], ".chrom.mzML"))
    chromHead <- readChromatogramHeader(mzmlName)
    chromatogramIdAsInteger(chromHead)
    # Merge chromatogram indices with transition indices and save them.
    # Following function merges analytesInfo dataframe with the chromatogram Header.
    mergeOswAnalytes_ChromHeader(oswAnalytes, chromHead, analyteFDR, runType)
    oswFiles[[i]] <- oswAnalytes
    message("Fetched chromatogram indices from ", filenames$filename[i])
  }
  # Assign rownames to the each element of list
  names(oswFiles) <- rownames(filenames)
  return(oswFiles)
}
