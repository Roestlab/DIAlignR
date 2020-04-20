#' Calculates LOESS fit between RT of two runs
#'
#' This function selects features from oswFiles which has m-score < maxFdrLoess. It fits LOESS on these feature.
#' Loess mapping is established from reference to experiment run.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-14
#' @importFrom rlang .data
#' @param oswFiles (list of data-frames) it is output from getFeatures function.
#' @param ref (string) Must be a combination of "run" and an iteger e.g. "run2".
#' @param eXp (string) Must be a combination of "run" and an iteger e.g. "run2".
#' @param maxFdrGlobal (numeric) A numeric value between 0 and 1. Features should have m-score lower than this value for participation in LOESS fit.
#' @param spanvalue (numeric) Spanvalue for LOESS fit. For targeted proteomics 0.1 could be used.
#' @importFrom dplyr %>%
#' @importFrom stats loess loess.control
#' @return An object of class "loess".
#' @seealso \code{\link{getLinearfit}, \link{getFeatures}}
#' @keywords internal
#' @examples
#' data(oswFiles_DIAlignR, package="DIAlignR")
#' \dontrun{
#' Loess.fit <- getLOESSfit(oswFiles = oswFiles_DIAlignR, ref = "run1", eXp = "run2",
#'  maxFdrGlobal = 0.05, spanvalue = 0.1)
#' }
getLOESSfit <- function(oswFiles, ref, eXp, maxFdrGlobal, spanvalue = 0.1){
  df.ref <-  oswFiles[[ref]] %>% dplyr::filter(.data$m_score <= maxFdrGlobal & .data$peak_group_rank == 1) %>%
    dplyr::select(.data$transition_group_id, .data$RT)
  df.eXp <-  oswFiles[[eXp]] %>% dplyr::filter(.data$m_score <= maxFdrGlobal & .data$peak_group_rank == 1) %>%
    dplyr::select(.data$transition_group_id, .data$RT)
  RUNS_RT <- dplyr::inner_join(df.ref, df.eXp, by = "transition_group_id", suffix = c(".ref", ".eXp"))
  Loess.fit <- loess(RT.eXp ~ RT.ref, data = RUNS_RT,
                     span = spanvalue,
                     control=loess.control(surface="direct"))
  # direct surface allows to extrapolate outside of training data boundary while using predict.
  Loess.fit
}

#' Calculates linear fit between RT of two runs
#'
#' This function selects features from oswFiles which has m-score < maxFdrLoess. It fits Linear model on these feature.
#' Loess mapping is established from reference to experiment run.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-14
#' @param oswFiles (list of data-frames) it is output from getFeatures function.
#' @param ref (string) Must be a combination of "run" and an iteger e.g. "run2".
#' @param eXp (string) Must be a combination of "run" and an iteger e.g. "run2".
#' @param maxFdrGlobal (numeric) A numeric value between 0 and 1. Features should have m-score lower than this value for participation in linear fit.
#' @importFrom dplyr %>%
#' @importFrom stats lm
#' @return An object of class "lm".
#' @seealso \code{\link{getLOESSfit}, \link{getFeatures}}
#' @keywords internal
#' @examples
#' data(oswFiles_DIAlignR, package="DIAlignR")
#' \dontrun{
#' lm.fit <- getLinearfit(oswFiles = oswFiles_DIAlignR, ref = "run1", eXp = "run2",
#'  maxFdrGlobal = 0.05)
#' }
getLinearfit <- function(oswFiles, ref, eXp, maxFdrGlobal){
  df.ref <-  oswFiles[[ref]] %>% dplyr::filter(.data$m_score <= maxFdrGlobal & .data$peak_group_rank == 1) %>%
    dplyr::select(.data$transition_group_id, .data$RT)
  df.eXp <-  oswFiles[[eXp]] %>% dplyr::filter(.data$m_score <= maxFdrGlobal & .data$peak_group_rank == 1) %>%
    dplyr::select(.data$transition_group_id, .data$RT)
  RUNS_RT <- dplyr::inner_join(df.ref, df.eXp, by = "transition_group_id", suffix = c(".ref", ".eXp"))
  # For testing we want to avoid validation peptides getting used in the fit.
  lm.fit <- lm(RT.eXp ~ RT.ref, data = RUNS_RT)
  lm.fit
}


#' Calculates global alignment between RT of two runs
#'
#' This function selects features from oswFiles which has m-score < maxFdrLoess. It fits linear/loess regression on these feature.
#' Retention-time mapping is established from reference to experiment run.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-14
#' @param oswFiles (list of data-frames) it is output from getFeatures function.
#' @param ref (string) Must be a combination of "run" and an iteger e.g. "run2".
#' @param eXp (string) Must be a combination of "run" and an iteger e.g. "run2".
#' @param fitType (string) Must be from "loess" or "linear".
#' @param maxFdrGlobal (numeric) A numeric value between 0 and 1. Features should have m-score lower than this value for participation in global fit.
#' @param spanvalue (numeric) Spanvalue for LOESS fit. For targeted proteomics 0.1 could be used.
#' @importFrom dplyr %>%
#' @return An object of class "loess".
#' @seealso \code{\link{getFeatures}}
#' @examples
#' data(oswFiles_DIAlignR, package="DIAlignR")
#' lm.fit <- getGlobalAlignment(oswFiles = oswFiles_DIAlignR, ref = "run1", eXp = "run2",
#'  fitType = "linear", maxFdrGlobal = 0.05, spanvalue = 0.1)
#' @export
getGlobalAlignment <- function(oswFiles, ref, eXp, fitType = "linear", maxFdrGlobal = 0.01, spanvalue = 0.1){
  if(fitType == "loess"){
    fit <- getLOESSfit(oswFiles, ref, eXp, maxFdrGlobal, spanvalue)
  }
  else{
    fit <- getLinearfit(oswFiles, ref, eXp, maxFdrGlobal)
  }
  fit
}


#' Calculates Residual Standard Error of the fit
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-04-08
#' @importFrom methods is
#' @param fit (lm or loess) Linear or loess fit object between reference and experiment run.
#' @return (numeric)
#' @seealso \code{\link{getLOESSfit}, \link{getLinearfit}, \link{getGlobalAlignment}}
#' @keywords internal
#' @examples
#' data(oswFiles_DIAlignR, package="DIAlignR")
#' \dontrun{
#' Loess.fit <- getGlobalAlignment(oswFiles = oswFiles_DIAlignR, ref = "run1", eXp = "run2",
#' maxFdrGlobal = 0.05, spanvalue = 0.1, fit = "loess")
#' getRSE(Loess.fit)
#' }
getRSE <- function(fit){
  if(is(fit, "loess")){
    RSE <- fit[["s"]]
  } else{
    RSE <- summary(fit)[["sigma"]]
  }
  RSE
}

#' Calculates all global alignment needed in refRun
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-04-19
#' @param refRun (data-frame) Output of getRefRun function. Must have two columsn : transition_group_id and run.
#' @param features (list of data-frames) it is output from getFeatures function.
#' @param fileInfo (data-frame) Output of getRunNames function.
#' @param globalAlignment (string) Must be from "loess" or "linear".
#' @param globalAlignmentFdr (numeric) A numeric value between 0 and 1. Features should have m-score lower than this value for participation in global fit.
#' @param globalAlignmentSpan (numeric) Spanvalue for LOESS fit. For targeted proteomics 0.1 could be used.
#' @return (list) Each element is either of class lm or loess.
#' @seealso \code{\link{getRefRun}, \link{getFeatures}, \link{getGlobalAlignment}}
#' @keywords internal
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' fileInfo <- getRunNames(dataPath, oswMerged = TRUE)
#' features <- getFeatures(fileInfo, maxFdrQuery = 0.05)
#' data("multipeptide_DIAlignR", package = "DIAlignR")
#' \dontrun{
#' refRun <- getRefRun(multipeptide_DIAlignR)
#' fits <- getGlobalFits(refRun, features, fileInfo, "linear", 0.05, 0.1)
#' }
getGlobalFits <- function(refRun, features, fileInfo, globalAlignment,
                          globalAlignmentFdr, globalAlignmentSpan){
  globalFits <- list()
  refs <- unique(refRun[["run"]])
  for(ref in refs){
    exps <- setdiff(rownames(fileInfo), ref)
    for(eXp in exps){
      pair <- paste(ref, eXp, sep = "_")
      globalFit <- getGlobalAlignment(features, ref, eXp,
                                      globalAlignment, globalAlignmentFdr, globalAlignmentSpan)
      globalFits[[pair]] <- globalFit
    }
  }
  globalFits
}
