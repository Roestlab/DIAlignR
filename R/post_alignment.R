#' Pick feature closest to reference peak
#'
#' It picks a feature that is within adaptiveRT window across eXpRT and has lowest m-score.
#' Feature's m-score also has to be smaller than featureFDR.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-13
#' @importFrom rlang .data
#' @param eXpRT (numeric) retention time in experiment run.
#' @param analyte (integer) vector of precursor IDs.
#' @param oswFiles (list of data-frames) it is output from getFeatures function.
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
#' @seealso \code{\link{getFeatures}}
#' @keywords internal
#' @examples
#' data(oswFiles_DIAlignR, package="DIAlignR")
#' \dontrun{
#' pickNearestFeature(eXpRT = 5237.8, analyte = 4618L, oswFiles = oswFiles_DIAlignR,
#'  runname = "run2", adaptiveRT = 77.82315, featureFDR = 0.05)
#' }
pickNearestFeature <- function(eXpRT, analyte, oswFiles, runname, adaptiveRT, featureFDR){
  # Fetch features for an analyte.
  df <- oswFiles[[runname]] %>% dplyr::filter(.data$transition_group_id == analyte) %>%
    dplyr::select(.data$leftWidth, .data$rightWidth, .data$RT, .data$intensity, .data$peak_group_rank, .data$m_score)
  # select features which are within adaptive RT of mapped retention time.
  df <- df %>% dplyr::filter(abs(.data$RT - eXpRT) <= adaptiveRT)
  # select highest peak_group_rank feature.
  df <- df %>% dplyr::filter(.data$m_score < featureFDR & .data$peak_group_rank == min(.data$peak_group_rank))
  if(df %>% nrow() == 0){
    return(NULL)
  }
  df <-  df %>% as.list()
  df
}


#' Establishes mapping from index to time
#'
#' Takes a time vector and index vector of same length. This function create a
#' new time vector given indices specified in \code{idx}.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-13
#' @param timeVec A numeric vector
#' @param idx An integer vector
#' @return A mutated time vector
#' @examples
#' timeVec <- c(1.3,5.6,7.8)
#' idx <- c(NA, NA, 1L, 2L, NA, NA, 3L, NA)
#' mapIdxToTime(timeVec, idx) # c(NA, NA, 1.3, 5.6, 6.333, 7.067, 7.8, NA)
#'
#' @importFrom zoo na.approx
#' @export
mapIdxToTime <- function(timeVec, idx){
  mutateT <- na.approx(timeVec[idx], na.rm = FALSE)
  return(mutateT)
}


#' Map reference run time on the experiment run.
#'
#' Retention time from reference run is mapped to experiment run using AlignObj.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-13
#'
#' @param refRT (numeric) retention time in reference run.
#' @param tVec.ref (numeric) time vector of refernce run.
#' @param tVec.eXp (numeric) time vector of experiment run.
#' @param AlignObj (S4 object) must have indexA_aligned, indexB_aligned and score as slots.
#' @return (numeric)
#' @keywords internal
#' @examples
#' data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
#' tVec.ref <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["run1"]][["14299_QFNNTDIVLLEDFQK/3"]][[1]][, "time"]
#' tVec.eXp <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["run2"]][["14299_QFNNTDIVLLEDFQK/3"]][[1]][, "time"]
#' \dontrun{
#' AlignObj <- testAlignObj()
#' mappedRTfromAlignObj(refRT= 5238.35, tVec.ref, tVec.eXp, AlignObj)
#' }
mappedRTfromAlignObj <- function(refRT, tVec.ref, tVec.eXp, AlignObj){
  AlignedIndices <- cbind(slot(AlignObj, "indexA_aligned"),
                          slot(AlignObj, "indexB_aligned"),
                          slot(AlignObj, "score"))
  colnames(AlignedIndices) <- c("indexAligned.ref", "indexAligned.eXp", "score")
  AlignedIndices <- AlignedIndices[(AlignedIndices[,"indexAligned.ref"] != 0L), ]
  AlignedIndices[, 1:2][AlignedIndices[, 1:2] == 0] <- NA
  tAligned.ref <- mapIdxToTime(tVec.ref, AlignedIndices[,"indexAligned.ref"])
  tAligned.eXp <- mapIdxToTime(tVec.eXp, AlignedIndices[,"indexAligned.eXp"])
  # Map retention time from reference to eXp.
  eXpRT <- tAligned.eXp[which.min(abs(tAligned.ref - refRT))]
  eXpRT
}


#' Set Alignment rank to the aligned feature
#'
#' Picks the top feature in run by comparing m-score to unaligned FDR and aligned FDR.
#' If no satisfactory feature is found, peak-integration is carried out by mapping left and right peak
#' boundaries from the reference feature and integrating area under the curve.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-04-13
#'
#' @param multipeptide (list of data-frames) each element of the list is collection of features associated with a precursor.
#' @param ref (numeric) name of the refernce run.
#' @param eXp (numeric) name of the experiment run.
#' @param analyte_chr (string) coerced transition_group_id.
#' @param unalignedFDR (numeric) must be between 0 and maxFdrQuery. Features below unalignedFDR are
#'  considered for quantification even without the RT alignment.
#' @param alignedFDR (numeric) must be between unalignedFDR and 1. Features below alignedFDR are
#'  considered for quantification after the alignment.
#' @param adaptiveRT (numeric) defines the window around the aligned retention time, within which
#'  features with m-score below aligned FDR are considered for quantification.
#' @param tAligned (list) the first element corresponds to the aligned reference time,
#'  the second element is the aligned experiment time.
#' @param XICs.ref (list) list of extracted ion chromatograms from reference run.
#' @param XICs.eXp (list) list of extracted ion chromatograms from experiment run.
#' @param integrationType (string) method to ompute the area of a peak contained in XICs. Must be
#'  from "intensity_sum", "trapezoid", "simpson".
#' @param baselineType (string) method to estimate the background of a peak contained in XICs. Must be
#'  from "base_to_base", "vertical_division_min", "vertical_division_max".
#' @param fitEMG (logical) enable/disable exponentially modified gaussian peak model fitting.
#' @param recalIntensity (logical) recalculate intensity for all analytes.
#' @param fillMissing (logical) calculate intensity for ananlytes for which features are not found.
#' @return (NULL)
#' @seealso \code{\link{getMultipeptide}, \link{calculateIntensity}}
#' @keywords internal
#'
#' @examples
#' data(multipeptide_DIAlignR, package="DIAlignR")
#' data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
#' XICs.ref <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["run1"]][["14299_QFNNTDIVLLEDFQK/3"]]
#' XICs.eXp <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["run2"]][["14299_QFNNTDIVLLEDFQK/3"]]
#' \dontrun{
#' # Use getAlignedIndices() to get tAligned.
#' setAlignmentRank(multipeptide_DIAlignR, ref = "run1", eXp = "run2", analyte_chr = "4618",
#'  unalignedFDR = 0.01, alignedFDR = 0.05, adaptiveRT = 30, tAligned, XICs.ref, XICs.eXp,
#'  integrationType = "intensity_sum", baselineType = "base_to_base", fitEMG = FALSE,
#'  recalIntensity = FALSE, fillMissing = TRUE)
#' multipeptide[["4618"]]
#' }
setAlignmentRank <- function(multipeptide, ref, eXp, analyte_chr, unalignedFDR, alignedFDR, adaptiveRT,
                             tAligned, XICs.ref, XICs.eXp, integrationType, baselineType, fitEMG,
                             recalIntensity, fillMissing = TRUE){
  # reference run.
  df <- multipeptide[[analyte_chr]]
  refIdx <- which(df[["run"]] == ref & df[["peak_group_rank"]] == 1)
  df[["alignment_rank"]][refIdx] <- 1L
  refRT <- df[["RT"]][refIdx]
  leftRef <- df[["leftWidth"]][refIdx]; rightRef <- df[["rightWidth"]][refIdx]
  if(recalIntensity){
    df[["intensity"]][refIdx] <-calculateIntensity(XICs.ref, leftRef, rightRef,
                                                   integrationType, baselineType, fitEMG)
    }

  # Experiment run.
  idx <- which(df[["run"]] == eXp)
  eXpRT <- tAligned[[2]][which.min(abs(tAligned[[1]] - refRT))]
  left <- tAligned[[2]][which.min(abs(tAligned[[1]] - leftRef))]
  right <- tAligned[[2]][which.min(abs(tAligned[[1]] - rightRef))]
  # TODO. Save for the edge cases. or use wider chromatogram.
  if(any(length(left)==0, length(right)==0, length(eXpRT)==0 )){
    return(NULL) # Can happen if XICs have all zero intensities.
  }

  featurePresent <- FALSE
  if(any(df[["m_score"]][idx] < unalignedFDR, na.rm = TRUE)){
    # Check if any feature is below unaligned FDR. If present alignment_rank = 1.
    featurePresent <- TRUE
  } else if(any(df[["m_score"]][idx] < alignedFDR, na.rm = TRUE)){
    # Check if any feature is below aligned FDR and within adaptive RT. If present alignment_rank = 1.
    idx <- idx[abs(df[["rightWidth"]][idx] - eXpRT) < adaptiveRT |
                     abs(df[["leftWidth"]][idx] - eXpRT) < adaptiveRT |
                     abs(df[["RT"]][idx] - eXpRT) < adaptiveRT]
    idx <- idx[!is.na(idx)]
    if(length(idx!=0)){
      featurePresent <- TRUE
    }
  }

  if(featurePresent){
    idx <- idx[which.min(df[["m_score"]][idx])]
    df[["alignment_rank"]][idx] <- 1L
    if(recalIntensity){
      df[["intensity"]][idx] <- calculateIntensity(XICs.eXp, left, right,
                                                     integrationType, baselineType, fitEMG)}
  } else if(fillMissing){
    # Otherwise create new feature and alignment rank = 1.
    intensity <- calculateIntensity(XICs.eXp, left, right, integrationType, baselineType, fitEMG)
    row <- data.frame("transition_group_id" = df[["transition_group_id"]][1], "RT" = eXpRT,
                      "intensity"= intensity, "leftWidth" = left, "rightWidth" = right,
                      "m_score" = NA_integer_, "peak_group_rank" = NA_real_, "run" = eXp,
                      "alignment_rank" = 1L)
    df <- rbind(df, row)
    }

  # Update the list element in the parent frame.
  assign("temp", df, envir = parent.frame(n = 1))
  with(parent.frame(n = 1), multipeptide[[analyte_chr]] <- temp)
  invisible(NULL)
}



#' Calculates area of a peak in XIC group
#'
#' Retention time from reference run is mapped to experiment run using AlignObj.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-04-13
#'
#' @param XICs (list) list of extracted ion chromatograms of a precursor.
#' @param left (numeric) left boundary of the peak.
#' @param right (numeric) right boundary of the peak.
#' @param integrationType (string) method to ompute the area of a peak contained in XICs. Must be
#'  from "intensity_sum", "trapezoid", "simpson".
#' @param baselineType (string) method to estimate the background of a peak contained in XICs. Must be
#'  from "base_to_base", "vertical_division_min", "vertical_division_max".
#' @param fitEMG (logical) enable/disable exponentially modified gaussian peak model fitting.
#' @return (numeric)
#' @keywords internal
#' @seealso \code{\link{getMultipeptide}, \link{setAlignmentRank}}
#' @examples
#' data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
#' XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["run1"]][["14299_QFNNTDIVLLEDFQK/3"]]
#' \dontrun{
#' calculateIntensity(XICs, 5220, 5261, integrationType = "intensity_sum",
#'  baselineType = "base_to_base", fitEMG = FALSE)
#' }
calculateIntensity <- function(XICs, left, right, integrationType, baselineType,
                               fitEMG){
  time <- lapply(XICs, `[[`, 1)
  intensityList <- lapply(XICs, `[[`, 2)
  intensity <- areaIntegrator(time, intensityList, left, right, integrationType, baselineType,
                              fitEMG)
  intensity
}
