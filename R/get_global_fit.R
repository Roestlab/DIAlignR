## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("m_score"))

#' Calculates LOESS fit between RT of two runs
#'
#' This function selects features from oswFiles which has m-score < maxFdrLoess. It then fit LOESS on these feature.
#' Loess mapping is established from reference to experiment run.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-14
#' @param oswFiles (list of data-frames) it is output from getOswFiles function.
#' @param ref (string) Must be a combination of "run" and an iteger e.g. "run2".
#' @param eXp (string) Must be a combination of "run" and an iteger e.g. "run2".
#' @param maxFdrGlobal (numeric) A numeric value between 0 and 1. Features should have m-score lower than this value for participation in LOESS fit.
#' @param spanvalue (numeric) Spanvalue for LOESS fit. For targeted proteomics 0.1 could be used.
#' @importFrom dplyr %>%
#' @importFrom stats loess loess.control
#' @return An object of class "loess".
#' @seealso \code{\link{getLinearfit}, \link{getOswFiles}}
#' @keywords internal
#' @examples
#' data(oswFiles_DIAlignR, package="DIAlignR")
#' \dontrun{
#' Loess.fit <- getLOESSfit(oswFiles = oswFiles_DIAlignR, ref = "run1", eXp = "run2",
#'  maxFdrGlobal = 0.05, spanvalue = 0.1)
#' }
getLOESSfit <- function(oswFiles, ref, eXp, maxFdrGlobal, spanvalue = 0.1){
  df.ref <-  oswFiles[[ref]] %>% dplyr::filter(m_score <= maxFdrGlobal & peak_group_rank == 1) %>%
    dplyr::select(transition_group_id, RT)
  df.eXp <-  oswFiles[[eXp]] %>% dplyr::filter(m_score <= maxFdrGlobal & peak_group_rank == 1) %>%
    dplyr::select(transition_group_id, RT)
  RUNS_RT <- dplyr::inner_join(df.ref, df.eXp, by = "transition_group_id", suffix = c(".ref", ".eXp"))
  Loess.fit <- loess(RT.eXp ~ RT.ref, data = RUNS_RT,
                     span = spanvalue,
                     control=loess.control(surface="direct"))
  # direct surface allows to extrapolate outside of training data boundary while using predict.
  Loess.fit
}

#' Calculates linear fit between RT of two runs
#'
#' This function selects features from oswFiles which has m-score < maxFdrLoess. It then fit Linear model on these feature.
#' Loess mapping is established from reference to experiment run.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-14
#' @param oswFiles (list of data-frames) it is output from getOswFiles function.
#' @param ref (string) Must be a combination of "run" and an iteger e.g. "run2".
#' @param eXp (string) Must be a combination of "run" and an iteger e.g. "run2".
#' @param maxFdrGlobal (numeric) A numeric value between 0 and 1. Features should have m-score lower than this value for participation in linear fit.
#' @importFrom dplyr %>%
#' @importFrom stats lm
#' @return An object of class "lm".
#' @seealso \code{\link{getLOESSfit}, \link{getOswFiles}}
#' @keywords internal
#' @examples
#' data(oswFiles_DIAlignR, package="DIAlignR")
#' \dontrun{
#' lm.fit <- getLinearfit(oswFiles = oswFiles_DIAlignR, ref = "run1", eXp = "run2",
#'  maxFdrGlobal = 0.05)
#' }
getLinearfit <- function(oswFiles, ref, eXp, maxFdrGlobal){
  df.ref <-  oswFiles[[ref]] %>% dplyr::filter(m_score <= maxFdrGlobal & peak_group_rank == 1) %>%
    dplyr::select(transition_group_id, RT)
  df.eXp <-  oswFiles[[eXp]] %>% dplyr::filter(m_score <= maxFdrGlobal & peak_group_rank == 1) %>%
    dplyr::select(transition_group_id, RT)
  RUNS_RT <- dplyr::inner_join(df.ref, df.eXp, by = "transition_group_id", suffix = c(".ref", ".eXp"))
  # For testing we want to avoid validation peptides getting used in the fit.
  lm.fit <- lm(RT.eXp ~ RT.ref, data = RUNS_RT)
  lm.fit
}


#' Calculates global alignment between RT of two runs
#'
#' This function selects features from oswFiles which has m-score < maxFdrLoess. It then fit linear/loess regression on these feature.
#' Retention-time mapping is established from reference to experiment run.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-14
#' @param oswFiles (list of data-frames) it is output from getOswFiles function.
#' @param ref (string) Must be a combination of "run" and an iteger e.g. "run2".
#' @param eXp (string) Must be a combination of "run" and an iteger e.g. "run2".
#' @param maxFdrGlobal (numeric) A numeric value between 0 and 1. Features should have m-score lower than this value for participation in global fit.
#' @param spanvalue (numeric) Spanvalue for LOESS fit. For targeted proteomics 0.1 could be used.
#' @param fitType (char) Must be from "loess" or "linear".
#' @importFrom dplyr %>%
#' @return An object of class "loess".
#' @seealso \code{\link{getOswFiles}}
#' @examples
#' data(oswFiles_DIAlignR, package="DIAlignR")
#' Loess.fit <- getGlobalAlignment(oswFiles = oswFiles_DIAlignR, ref = "run1", eXp = "run2",
#'  maxFdrGlobal = 0.05, spanvalue = 0.1, fit = "loess")
#' @export
getGlobalAlignment <- function(oswFiles, ref, eXp, maxFdrGlobal, spanvalue = 0.1, fitType = "loess"){
  if(fitType == "loess"){
    fit <- getLOESSfit(oswFiles, ref, eXp, maxFdrGlobal, spanvalue)
  }
  else{
    fit <- getLinearfit(oswFiles, ref, eXp, maxFdrGlobal)
  }
  fit
}
