#' Pick feature closest to reference peak
#'
#' @return Invisible NULL
pickNearestFeature <- function(eXpRT, analyte, oswFiles, runname, adaptiveRT, featureFDR){
  # Fetch features for an analyte.
  df <- oswFiles[[runname]] %>% dplyr::filter(transition_group_id == analyte) %>%
    dplyr::select(leftWidth, rightWidth, RT, Intensity, peak_group_rank, m_score)
  # select features which are within adaptive RT of mapped retention time.
  df <- df %>% dplyr::filter(abs(RT - eXpRT) <= adaptiveRT)
  # select highest peak_group_rank feature.
  df <- df %>% dplyr::filter(m_score < featureFDR & peak_group_rank == min(peak_group_rank))
  if(df %>% nrow() == 0){
    return(NULL)
  }
  df <-  df %>% as.list()
  df
}
