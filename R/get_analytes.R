#' Outputs analytes below FDR
#'
#' Provides all found analytes or only common analytes that have m-score less than analyteFDR.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-04-14
#' @importFrom dplyr %>% pull
#' @importFrom rlang .data
#' @param features (list of data-frames) contains features and their properties identified in each run.
#' @param analyteFDR (numeric) only analytes that have m-score less than this, will be included in the output.
#' @param commonAnalytes (logical) TRUE: intersect across all runs, FASLE: union across all runs.
#' @return a vector of integer.
#' @seealso \code{\link{getFeatures}}
#' @keywords internal
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' fileInfo <- getRunNames(dataPath = dataPath)
#' features <- getFeatures(fileInfo, maxFdrQuery = 1.00, runType = "DIA_proteomics")
#' \dontrun{
#' commonAnalytes <- analytesFromFeatures(features, analyteFDR = 0.01, commonAnalytes = TRUE)
#' }
analytesFromFeatures <- function(features, analyteFDR = 1.00, commonAnalytes = TRUE){
  analytes <- c()
  analytes <- dplyr::filter(features[[1]], .data$m_score < analyteFDR) %>%
    pull(.data$transition_group_id) %>% union(analytes)
  if(commonAnalytes){
    # Get intersect
    for(df in features){
      analytes <- dplyr::filter(df, .data$m_score < analyteFDR) %>%
        pull(.data$transition_group_id) %>% dplyr::intersect(analytes)
    }
  } else {
    # Get union
    for(df in features){
        analytes <- dplyr::filter(df, .data$m_score < analyteFDR) %>%
          pull(.data$transition_group_id) %>% dplyr::union(analytes)
    }
  }
  if(length(analytes) == 0){
    return(NULL)
  }
  unique(analytes)
}
