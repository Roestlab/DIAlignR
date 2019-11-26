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


#' Establishes mapping back from index to time
#'
#' Takes a time vector and index vector of same length. This function create a
#' new time vector given indices specified in \code{idx}. For \code{NA} indices
#' it uses last index to fill time value. For \code{NA} appearing at the start
#' of \code{idx}, it uses next index to get time value.
#' @param timeVec A numeric vector
#' @param idx An integer vector
#' @return A mutated time vector
#' @export
mapIdxToTime <- function(timeVec, idx){
  mutateT <- zoo::na.locf(na.locf(timeVec[idx], na.rm = FALSE), fromLast = TRUE)
  return(mutateT)
}

#' Use AlignObj to map reference retention time on experiment run.
#'
#' @return A numeric
mappedRTfromAlignObj <- function(refRT, tVec.ref, tVec.eXp, AlignObj){
  AlignedIndices <- cbind(AlignObj@indexA_aligned,
                          AlignObj@indexB_aligned,
                          AlignObj@score)
  colnames(AlignedIndices) <- c("indexAligned.ref", "indexAligned.eXp", "score")
  AlignedIndices <- AlignedIndices[(AlignedIndices[,"indexAligned.ref"] != 0L), ]
  AlignedIndices[, 1:2][AlignedIndices[, 1:2] == 0] <- NA
  tAligned.ref <- mapIdxToTime(tVec.ref, AlignedIndices[,"indexAligned.ref"])
  tAligned.eXp <- mapIdxToTime(tVec.eXp, AlignedIndices[,"indexAligned.eXp"])
  # Map retention time from reference to eXp.
  eXpRT <- tAligned.eXp[which.min(abs(tAligned.ref - refRT))]
  eXpRT
}
