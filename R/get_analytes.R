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
