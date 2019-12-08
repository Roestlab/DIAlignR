#' Get names of analytes found in all runs.
#'
#' @importFrom dplyr %>%
#' @return A vector of strings.
getAnalytesName <- function(oswFiles, analyteFDR = 1.00, commonAnalytes = TRUE){
  analytes <- c()
  analytes <- dplyr::filter(oswFiles[[1]], m_score < analyteFDR) %>%
    .$transition_group_id %>% union(analytes)
  if(commonAnalytes){
    # Get intersect
    for(oswAnalytes in oswFiles){
      analytes <- dplyr::filter(oswAnalytes, m_score < analyteFDR) %>%
        .$transition_group_id %>% dplyr::intersect(analytes)
    }
  } else {
    # Get union
    for(oswAnalytes in oswFiles){
        analytes <- dplyr::filter(oswAnalytes, m_score < analyteFDR) %>%
        .$transition_group_id %>% dplyr::union(analytes)
    }
  }
  if(length(analytes) == 0){
    return(NULL)
  }
  analytes
}

#' Get names of analytes found in all runs.
#'
#' @importFrom dplyr %>%
#' @return A vector of strings.
#' @export
getAnalytes <- function(dataPath, runs = NULL, oswMerged = TRUE, runType = "DIA_Proteomics",
                        commonAnalytes = FALSE, maxFdrQuery = 0.05, nameCutPattern = "(.*)(/)(.*)",
                        analyteInGroupLabel = FALSE){
  # Get filenames from .merged.osw file and check if names are consistent between osw and mzML files.
  filenames <- getRunNames(dataPath, oswMerged, nameCutPattern)
  if(!is.null(runs)){
    filenames <- filenames[filenames$runs %in% runs,]
    missingRun <- setdiff(runs, filenames$runs)
    if(length(missingRun) != 0){
      return(stop(missingRun, " runs are not found."))
    }
  }

  # Get Precursors from the query and respective chromatogram indices.
  oswFiles <- getOswAnalytes(dataPath = dataPath, filenames = filenames, oswMerged = oswMerged,
                             maxFdrQuery = maxFdrQuery, runType = runType,
                             analyteInGroupLabel = analyteInGroupLabel)
  analytes <- getAnalytesName(oswFiles, analyteFDR = maxFdrQuery, commonAnalytes = commonAnalytes)
  message(length(analytes), " analytes are found.")
  analytes
}
