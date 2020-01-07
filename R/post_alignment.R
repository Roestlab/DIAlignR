#' Pick feature closest to reference peak
#'
#' It picks a feature that is within adaptiveRT window across eXpRT and has lowest m-score. Feature's m-score also has to be smaller than featureFDR.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + MIT
#' Date: 2019-12-13
#' @param eXpRT (numeric) retention time in experiment run.
#' @param analyte (string) analyte is as PRECURSOR.GROUP_LABEL or as PEPTIDE.MODIFIED_SEQUENCE and PRECURSOR.CHARGE from osw file.
#' @param oswFiles (list of data-frames) it is output from getOswFiles function.
#' @param runname (string) must be a combination of "run" and an iteger e.g. "run2".
#' @param adaptiveRT (numeric) half-width of retention time window. Feature, if found, is picked from within this window.
#' @param featureFDR (numeric) upper m-score cut-off for a feature to be picked.
#' @return (list) Following elements are present in the list:
#' \item{leftWidth}{(numeric) as in FEATURE.LEFT_WIDTH of osw files.}
#' \item{rightWidth}{(numeric) as in FEATURE.RIGHT_WIDTH of osw files.}
#' \item{RT}{(numeric) retention time as in FEATURE.EXP_RT of osw files.}
#' \item{Intensity}{(numeric) peak intensity as in FEATURE_MS2.AREA_INTENSITY of osw files.}
#' \item{peak_group_rank}{(integer) rank of each feature associated with transition_group_id.}
#' \item{m_score}{(numeric) q-value of each feature associated with transition_group_id.}
#'
#' @seealso \code{\link{getOswFiles}, \link{getOswAnalytes}}
#' @examples
#' data(oswFiles_DIAlignR, package="DIAlignR")
#' \dontrun{
#' pickNearestFeature(eXpRT = 5237.8, analyte = "14299_QFNNTDIVLLEDFQK/3",
#'  oswFiles = oswFiles_DIAlignR, runname = "run2", adaptiveRT = 77.82315, featureFDR = 0.05)
#' }
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


#' Establishes mapping from index to time
#'
#' Takes a time vector and index vector of same length. This function create a
#' new time vector given indices specified in \code{idx}. For \code{NA} indices
#' it uses last index to fill time value. For \code{NA} appearing at the start
#' of \code{idx}, it uses next index to get time value.
#' @importFrom dplyr %>%
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + MIT
#' Date: 2019-12-13
#' @param timeVec A numeric vector
#' @param idx An integer vector
#' @return A mutated time vector
#' @examples
#' timeVec <- c(1.3,5.6,7.8)
#' idx <- c(NA, NA, 1L, 2L, NA, NA, 3L, NA)
#' mapIdxToTime(timeVec, idx) # c(1.3, 1.3, 1.3, 5.6, 5.6, 5.6, 7.8, 7.8)
#'
#' @importFrom zoo na.locf
#' @export
mapIdxToTime <- function(timeVec, idx){
  mutateT <- na.locf(na.locf(timeVec[idx], na.rm = FALSE), fromLast = TRUE)
  return(mutateT)
}


#' Map reference run time on the experiment run.
#'
#' Retention time from reference run is mapped to experiment run using AlignObj.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + MIT
#' Date: 2019-12-13
#'
#' @param refRT (numeric) retention time in reference run.
#' @param tVec.ref (numeric) time vector of refernce run.
#' @param tVec.eXp (numeric) time vector of experiment run.
#' @param AlignObj (S4 object)
#' @return (numeric)
#' @examples
#' data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
#' tVec.ref <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["run1"]][["14299_QFNNTDIVLLEDFQK/3"]][[1]][, "time"]
#' tVec.eXp <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["run2"]][["14299_QFNNTDIVLLEDFQK/3"]][[1]][, "time"]
#' \dontrun{
#' AlignObj <- testAlignObj()
#' mappedRTfromAlignObj(refRT= 5238.35, tVec.ref, tVec.eXp, AlignObj)
#' }
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
