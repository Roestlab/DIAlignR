#' calculates LOESS fit between RT of two runs
#'
#' This function takes in run-pairs, names of test peptides, span paprameter for
#' LOESS function, output of OpenSWATH which include estimated retention time of
#' peptides.
#' @param pairName Names of runs joined with underscore. e.g. runA_runB, This will allow
#' alignment of runB to runA.
#' @param peptides Test peptides' names.
#' @param oswOutput list of list (OpenSWATH output).
#' @param spanvalue A numeric Spanvalue for LOESS fit. For targeted proteomics
#'   0.1 could be used.
#' @importFrom dplyr %>%
#' @return An object of class "loess".
getLOESSfit <- function(oswFiles, ref, eXp, maxFdrLoess, spanvalue = 0.1){
  df.ref <-  oswFiles[[ref]] %>% dplyr::filter(m_score <= maxFdrLoess & peak_group_rank == 1) %>%
    dplyr::select(transition_group_id, RT)
  df.eXp <-  oswFiles[[eXp]] %>% dplyr::filter(m_score <= maxFdrLoess & peak_group_rank == 1) %>%
    dplyr::select(transition_group_id, RT)
  RUNS_RT <- dplyr::inner_join(df.ref, df.eXp, by = "transition_group_id", suffix = c(".ref", ".eXp"))
  Loess.fit <- loess(RT.eXp ~ RT.ref, data = RUNS_RT,
                     span = spanvalue,
                     control=loess.control(surface="direct"))
  # direct surface allows to extrapolate outside of training data boundary while using predict.
  Loess.fit
}

#' calculates linear fit between RT of two runs
#'
#' This function takes in run-pairs, names of test peptides, output of OpenSWATH
#' which include estimated retention time of peptides.
#' @param pairName Names of runs joined with underscore. e.g. runA_runB, This will allow
#' alignment of runB to runA.
#' @param peptides Test peptides' names.
#' @param oswOutput list of list (OpenSWATH output).
#' @return  An object of class "lm".
#' @importFrom dplyr %>%
#' @export
getLinearfit <- function(oswFiles, ref, eXp, maxFdrLoess){
  df.ref <-  oswFiles[[ref]] %>% dplyr::filter(m_score <= maxFdrLoess & peak_group_rank == 1) %>%
    dplyr::select(transition_group_id, RT)
  df.eXp <-  oswFiles[[eXp]] %>% dplyr::filter(m_score <= maxFdrLoess & peak_group_rank == 1) %>%
    dplyr::select(transition_group_id, RT)
  RUNS_RT <- dplyr::inner_join(df.ref, df.eXp, by = "transition_group_id", suffix = c(".ref", ".eXp"))
  # For testing we want to avoid validation peptides getting used in the fit.
  lm.fit <- lm(RT.eXp ~ RT.ref, data = RUNS_RT)
  lm.fit
}
