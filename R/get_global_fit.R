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
#' @return An object of class "loess".
#' @export
getLOESSfit <- function(oswRef, oswExp, peptides, spanvalue = 0.1){
  run_pair <- strsplit(pairName, split = "_")[[1]]
  RUN1 <- oswOutput[[run_pair[1]]]; RUN2 <- oswOutput[[run_pair[2]]]
  cmp <- intersect(oswRef[,1], oswExp[,1]) # First column corresponds to transition_group_record
  oswRef <- oswRef[which(oswRef[,1] %in% cmp), ]
  oswExp <- oswExp[which(oswExp[,1] %in% cmp), ]
  oswRef <- oswRef[match(cmp, oswRef[,1]),]
  oswExp <- oswExp[match(cmp, oswExp[,1]),]
  RUNS_RT <- data.frame( "transition_group_record" = RUN1[,1], "RUN1" = oswRef$RT, "RUN2" = oswExp$RT)
  RUNS_RT <- RUNS_RT[order(RUNS_RT$RUN1), ]
  testPeptides <-intersect(cmp, peptides)
  # For testing we want to avoid validation peptides getting used in the fit.
  Loess.fit <- loess(RUN2 ~ RUN1, data = RUNS_RT,
                     subset = !transition_group_record %in% testPeptides,
                     span = spanvalue,
                     control=loess.control(surface="direct"))
  # direct surface allows to extrapolate outside of training data boundary while using predict.
  return(Loess.fit)
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
#' @export
getLinearfit <- function(pairName, peptides, oswOutput){
  run_pair <- strsplit(pairName, split = "_")[[1]]
  RUN1 <- oswOutput[[run_pair[1]]]; RUN2 <- oswOutput[[run_pair[2]]]
  cmp <- intersect(RUN1[,1], RUN2[,1]) # First column corresponds to transition_group_record
  RUN1 <- RUN1[which(RUN1[,1] %in% cmp), ]
  RUN2 <- RUN2[which(RUN2[,1] %in% cmp), ]
  RUN1 <- RUN1[match(cmp, RUN1[,1]),]
  RUN2 <- RUN2[match(cmp, RUN2[,1]),]
  RUNS_RT <- data.frame( "transition_group_record" = RUN1[,1], "RUN1" = RUN1$RT, "RUN2" = RUN2$RT)
  RUNS_RT <- RUNS_RT[order(RUNS_RT$RUN1), ]
  testPeptides <-intersect(cmp, peptides)
  # For testing we want to avoid validation peptides getting used in the fit.
  lm.fit <- lm(RUN2 ~ RUN1, data = RUNS_RT, subset = !transition_group_record %in% testPeptides)
  return(lm.fit)
}
