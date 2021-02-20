#' Outputs AlignObj from an alignment of two XIC-groups
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-13
#'
#' @param XICs.ref List of extracted ion chromatograms from reference run.
#' @param XICs.eXp List of extracted ion chromatograms from experiment run.
#' @param globalFit Linear or loess fit object between reference and experiment run.
#' @param alignType Available alignment methods are "global", "local" and "hybrid".
#' @param adaptiveRT (numeric) Similarity matrix is not penalized within adaptive RT.
#' @param normalization (character) Must be selected from "mean", "l2".
#' @param simType (string) Must be selected from dotProduct, cosineAngle, crossCorrelation,
#' cosine2Angle, dotProductMasked, euclideanDist, covariance and correlation.
#' @param goFactor (numeric) Penalty for introducing first gap in alignment. This value is multiplied by base gap-penalty.
#' @param geFactor (numeric) Penalty for introducing subsequent gaps in alignment. This value is multiplied by base gap-penalty.
#' @param cosAngleThresh (numeric) In simType = dotProductMasked mode, angular similarity should be higher than cosAngleThresh otherwise similarity is forced to zero.
#' @param OverlapAlignment (logical) An input for alignment with free end-gaps. False: Global alignment, True: overlap alignment.
#' @param dotProdThresh (numeric) In simType = dotProductMasked mode, values in similarity matrix higher than dotProdThresh quantile are checked for angular similarity.
#' @param gapQuantile (numeric) Must be between 0 and 1. This is used to calculate base gap-penalty from similarity distribution.
#' @param kerLen (integer) In simType = crossCorrelation, length of the kernel used to sum similarity score. Must be an odd number.
#' @param hardConstrain (logical) If FALSE; indices farther from noBeef distance are filled with distance from linear fit line.
#' @param samples4gradient (numeric) This parameter modulates penalization of masked indices.
#' @param objType (char) Must be selected from light, medium and heavy.
#' @return A S4 object. Three most-important slots are:
#' \item{indexA_aligned}{(integer) aligned indices of reference run.}
#' \item{indexB_aligned}{(integer) aligned indices of experiment run.}
#' \item{score}{(numeric) cumulative score of alignment.}
#'
#' @seealso \code{\link{alignChromatogramsCpp}}
#' @examples
#' data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
#' data(oswFiles_DIAlignR, package="DIAlignR")
#' run1 <- "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"
#' run2 <- "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"
#' XICs.ref <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[[run1]][["4618"]]
#' XICs.eXp <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[[run2]][["4618"]]
#' RUNS_RT <- getRTdf(oswFiles_DIAlignR, ref = "run1", eXp = "run2", maxFdrGlobal = 0.05)
#' globalFit <- loess(RT.eXp ~ RT.ref, data = RUNS_RT, span = 0.1, control=loess.control(surface="direct"))
#' AlignObj <- getAlignObj(XICs.ref, XICs.eXp, globalFit, alignType = "hybrid", adaptiveRT = 77.82315,
#'  normalization = "mean", simType = "dotProductMasked", goFactor = 0.125,
#'   geFactor = 40, cosAngleThresh = 0.3, OverlapAlignment = TRUE, dotProdThresh = 0.96,
#'   gapQuantile = 0.5, kerLen = 9L, hardConstrain = FALSE, samples4gradient = 100, objType = "light")
#' @export
getAlignObj <- function(XICs.ref, XICs.eXp, globalFit, alignType, adaptiveRT,
                        normalization, simType, goFactor, geFactor,
                        cosAngleThresh, OverlapAlignment,
                        dotProdThresh, gapQuantile, kerLen, hardConstrain,
                        samples4gradient, objType = "light"){
  XICs.ref1 <- xicIntersect(XICs.ref) # Fixed common time in fragment-ions
  XICs.eXp1 <- xicIntersect(XICs.eXp) # Fixed common time in fragment-ions

  tVec.ref <- XICs.ref1[[1]][, "time"] # Extracting time component
  tVec.eXp <- XICs.eXp1[[1]][, "time"] # Extracting time component
  len <- length(tVec.ref)
  B1p <- stats::predict(globalFit, data.frame("RT.ref" = tVec.ref[1]))[[1]]
  B2p <- stats::predict(globalFit, data.frame("RT.ref" = tVec.ref[len]))[[1]]

  # Set up constraints for penalizing similarity matrix
  samplingTime <- (tVec.ref[len] - tVec.ref[1])/(len-1)
  noBeef <- ceiling(adaptiveRT/samplingTime)

  # Perform dynamic programming for chromatogram alignment
  intensityList.ref <- lapply(XICs.ref1, `[`, i =, j = 2) # Extracting intensity values
  intensityList.eXp <- lapply(XICs.eXp1, `[`, i =, j = 2) # Extracting intensity values
  AlignObj <- alignChromatogramsCpp(intensityList.ref, intensityList.eXp,
                                    alignType, tVec.ref, tVec.eXp,
                                    normalization = normalization, simType = simType,
                                    B1p = B1p, B2p = B2p, noBeef = noBeef,
                                    goFactor = goFactor, geFactor = geFactor,
                                    cosAngleThresh = cosAngleThresh, OverlapAlignment = OverlapAlignment,
                                    dotProdThresh = dotProdThresh, gapQuantile = gapQuantile, kerLen = kerLen,
                                    hardConstrain = hardConstrain, samples4gradient = samples4gradient,
                                    objType = objType)
  AlignObj
}


getAlignObj2 <- function(XICs.ref, XICs.eXp, globalFit, adaptiveRT, params,
                         objType = "light"){
  XICs.ref1 <- xicIntersect(XICs.ref) # Fixed common time in fragment-ions
  XICs.eXp1 <- xicIntersect(XICs.eXp) # Fixed common time in fragment-ions

  tVec.ref <- XICs.ref1[[1]][, "time"] # Extracting time component
  tVec.eXp <- XICs.eXp1[[1]][, "time"] # Extracting time component
  len <- length(tVec.ref)
  B1p <- getPredict(globalFit, tVec.ref[1], params[["globalAlignment"]])
  B2p <- getPredict(globalFit, tVec.ref[len], params[["globalAlignment"]])

  # Set up constraints for penalizing similarity matrix
  samplingTime <- (tVec.ref[len] - tVec.ref[1])/(len-1)
  noBeef <- ceiling(adaptiveRT/samplingTime)

  # Perform dynamic programming for chromatogram alignment
  intensityList.ref <- lapply(XICs.ref1, `[`, i =, j = 2) # Extracting intensity values
  intensityList.eXp <- lapply(XICs.eXp1, `[`, i =, j = 2) # Extracting intensity values
  AlignObj <- alignChromatogramsCpp(intensityList.ref, intensityList.eXp,
                                    params[["alignType"]], tVec.ref, tVec.eXp,
                                    params[["normalization"]], params[["simMeasure"]],
                                    B1p = B1p, B2p = B2p, noBeef = noBeef,
                                    params[["goFactor"]], params[["geFactor"]],
                                    params[["cosAngleThresh"]], params[["OverlapAlignment"]],
                                    params[["dotProdThresh"]], params[["gapQuantile"]], params[["kerLen"]],
                                    params[["hardConstrain"]], params[["samples4gradient"]],
                                    objType = objType)
  AlignObj
}

#' Get mapping of reference RT on experiment run.
#'
#' This function aligns XICs of reference and experiment runs. Using alignment, it maps retention time from refernce run on experiment run.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-13
#' @param refRT Peak's retention-time in reference run.
#' @param XICs.ref List of extracted ion chromatograms from reference run.
#' @param XICs.eXp List of extracted ion chromatograms from experiment run.
#' @param globalFit Linear or loess fit object between reference and experiment run.
#' @param alignType Available alignment methods are "global", "local" and "hybrid".
#' @param adaptiveRT (numeric) Similarity matrix is not penalized within adaptive RT.
#' @param normalization (character) Must be selected from "mean", "l2".
#' @param simMeasure (string) Must be selected from dotProduct, cosineAngle, crossCorrelation,
#' cosine2Angle, dotProductMasked, euclideanDist, covariance and correlation.
#' @param goFactor (numeric) Penalty for introducing first gap in alignment. This value is multiplied by base gap-penalty.
#' @param geFactor (numeric) Penalty for introducing subsequent gaps in alignment. This value is multiplied by base gap-penalty.
#' @param cosAngleThresh (numeric) In simType = dotProductMasked mode, angular similarity should be higher than cosAngleThresh otherwise similarity is forced to zero.
#' @param OverlapAlignment (logical) An input for alignment with free end-gaps. False: Global alignment, True: overlap alignment.
#' @param dotProdThresh (numeric) In simType = dotProductMasked mode, values in similarity matrix higher than dotProdThresh quantile are checked for angular similarity.
#' @param gapQuantile (numeric) Must be between 0 and 1. This is used to calculate base gap-penalty from similarity distribution.
#' @param kerLen (integer) In simType = crossCorrelation, length of the kernel used to sum similarity score. Must be an odd number.
#' @param hardConstrain (logical) If FALSE; indices farther from noBeef distance are filled with distance from linear fit line.
#' @param samples4gradient (numeric) This parameter modulates penalization of masked indices.
#' @param objType (char) Must be selected from light, medium and heavy.
#' @return (numeric)
#' @seealso \code{\link{alignChromatogramsCpp}}
#' @keywords internal
#' @examples
#' data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
#' data(oswFiles_DIAlignR, package="DIAlignR")
#' run1 <- "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"
#' run2 <- "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"
#' XICs.ref <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[[run1]][["4618"]]
#' XICs.eXp <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[[run2]][["4618"]]
#' RUNS_RT <- getRTdf(oswFiles_DIAlignR, ref = "run1", eXp = "run2", maxFdrGlobal = 0.05)
#' globalFit <- loess(RT.eXp ~ RT.ref, data = RUNS_RT, span = 0.1, control=loess.control(surface="direct"))
#' adaptiveRT <- 77.82315 #3.5*globalFit$s
#' \dontrun{
#' getMappedRT(refRT = 5238.35, XICs.ref, XICs.eXp, globalFit, alignType = "hybrid",
#'  adaptiveRT = adaptiveRT, normalization = "mean",
#'   simMeasure = "dotProductMasked", goFactor = 0.125, geFactor = 40, cosAngleThresh = 0.3,
#'   OverlapAlignment = TRUE, dotProdThresh = 0.96, gapQuantile = 0.5, kerLen = 9L, hardConstrain = FALSE,
#'   samples4gradient = 100)
#' }
getMappedRT <- function(refRT, XICs.ref, XICs.eXp, globalFit, alignType, adaptiveRT,
                        normalization, simMeasure, goFactor, geFactor, cosAngleThresh,
                        OverlapAlignment, dotProdThresh, gapQuantile, kerLen, hardConstrain,
                        samples4gradient, objType = "light"){
  AlignObj <- getAlignObj(XICs.ref, XICs.eXp, globalFit, alignType, adaptiveRT,
                          normalization, simType = simMeasure, goFactor, geFactor,
                          cosAngleThresh, OverlapAlignment, dotProdThresh, gapQuantile,
                          kerLen, hardConstrain, samples4gradient, objType)
  tVec.ref <- XICs.ref[[1]][["time"]] # Extracting time component
  tVec.eXp <- XICs.eXp[[1]][["time"]] # Extracting time component
  eXpRT <- mappedRTfromAlignObj(refRT, tVec.ref, tVec.eXp, AlignObj)
  eXpRT
}

#' Get aligned Retention times.
#'
#' This function aligns XICs of reference and experiment runs.
#' It produces aligned retention times between refernce run and experiment run.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-13
#' @param XICs.ref List of extracted ion chromatograms from reference run.
#' @param XICs.eXp List of extracted ion chromatograms from experiment run.
#' @param globalFit Linear or loess fit object between reference and experiment run.
#' @param alignType Available alignment methods are "global", "local" and "hybrid".
#' @param adaptiveRT (numeric) Similarity matrix is not penalized within adaptive RT.
#' @param normalization (character) Must be selected from "mean", "l2".
#' @param simMeasure (string) Must be selected from dotProduct, cosineAngle, crossCorrelation,
#' cosine2Angle, dotProductMasked, euclideanDist, covariance and correlation.
#' @param goFactor (numeric) Penalty for introducing first gap in alignment. This value is multiplied by base gap-penalty.
#' @param geFactor (numeric) Penalty for introducing subsequent gaps in alignment. This value is multiplied by base gap-penalty.
#' @param cosAngleThresh (numeric) In simType = dotProductMasked mode, angular similarity should be higher than cosAngleThresh otherwise similarity is forced to zero.
#' @param OverlapAlignment (logical) An input for alignment with free end-gaps. False: Global alignment, True: overlap alignment.
#' @param dotProdThresh (numeric) In simType = dotProductMasked mode, values in similarity matrix higher than dotProdThresh quantile are checked for angular similarity.
#' @param gapQuantile (numeric) Must be between 0 and 1. This is used to calculate base gap-penalty from similarity distribution.
#' @param kerLen (integer) In simType = crossCorrelation, length of the kernel used to sum similarity score. Must be an odd number.
#' @param hardConstrain (logical) If FALSE; indices farther from noBeef distance are filled with distance from linear fit line.
#' @param samples4gradient (numeric) This parameter modulates penalization of masked indices.
#' @param objType (char) Must be selected from light, medium and heavy.
#' @return (list) the first element corresponds to the aligned reference time, the second element is the aligned experiment time.
#' @seealso \code{\link{alignChromatogramsCpp}, \link{getAlignObj}}
#' @examples
#' data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
#' data(oswFiles_DIAlignR, package="DIAlignR")
#' run1 <- "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"
#' run2 <- "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"
#' XICs.ref <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[[run1]][["4618"]]
#' XICs.eXp <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[[run2]][["4618"]]
#' RUNS_RT <- getRTdf(oswFiles_DIAlignR, ref = "run1", eXp = "run2", maxFdrGlobal = 0.05)
#' globalFit <- loess(RT.eXp ~ RT.ref, data = RUNS_RT, span = 0.1, control=loess.control(surface="direct"))
#' adaptiveRT <- 77.82315 #3.5*globalFit$s
#' getAlignedTimes(XICs.ref, XICs.eXp, globalFit, alignType = "hybrid",
#'  adaptiveRT = adaptiveRT, normalization = "mean",
#'   simMeasure = "dotProductMasked", goFactor = 0.125, geFactor = 40, cosAngleThresh = 0.3,
#'   OverlapAlignment = TRUE, dotProdThresh = 0.96, gapQuantile = 0.5, kerLen = 9L, hardConstrain = FALSE,
#'   samples4gradient = 100)
#' @export
getAlignedTimes <- function(XICs.ref, XICs.eXp, globalFit, alignType, adaptiveRT,
                        normalization, simMeasure, goFactor, geFactor, cosAngleThresh,
                        OverlapAlignment, dotProdThresh, gapQuantile, kerLen, hardConstrain,
                        samples4gradient, objType = "light"){
  alignedIndices <- getAlignedIndices(XICs.ref, XICs.eXp, globalFit, alignType, adaptiveRT,
                          normalization, simMeasure, goFactor, geFactor,
                          cosAngleThresh, OverlapAlignment, dotProdThresh, gapQuantile,
                          kerLen, hardConstrain, samples4gradient, objType)
  keep <- !is.na(alignedIndices[,"indexAligned.ref"])
  tVec.ref <- XICs.ref[[1]][,"time"] # Extracting time component
  tVec.eXp <- XICs.eXp[[1]][, "time"] # Extracting time component
  tAligned.ref <- mapIdxToTime(tVec.ref, alignedIndices[,"indexAligned.ref"])
  tAligned.eXp <- mapIdxToTime(tVec.eXp, alignedIndices[,"indexAligned.eXp"])
  tAligned.ref <- tAligned.ref[keep]
  tAligned.eXp <- tAligned.eXp[keep]
  list(tAligned.ref, tAligned.eXp)
}

#' Get aligned Retention times.
#'
#' This function aligns XICs of reference and experiment runs.
#' It produces aligned retention times between refernce run and experiment run.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2021) + GPL-3
#' Date: 2021-01-02
#' @inheritParams checkParams
#' @param XICs.ref List of extracted ion chromatograms from reference run.
#' @param XICs.eXp List of extracted ion chromatograms from experiment run.
#' @param globalFit Linear or loess fit object between reference and experiment run.
#' @param adaptiveRT (numeric) Similarity matrix is not penalized within adaptive RT.
#' @return (matrix) the first column corresponds to the aligned reference time, the second column is the aligned experiment time.
#' @seealso \code{\link{alignChromatogramsCpp}, \link{getAlignObj}}
#' @examples
#' data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
#' data(oswFiles_DIAlignR, package="DIAlignR")
#' run1 <- "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"
#' run2 <- "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"
#' XICs.ref <- lapply(XIC_QFNNTDIVLLEDFQK_3_DIAlignR[[run1]][["4618"]], as.matrix)
#' XICs.eXp <- lapply(XIC_QFNNTDIVLLEDFQK_3_DIAlignR[[run2]][["4618"]], as.matrix)
#' params <- paramsDIAlignR()
#' params[["globalAlignment"]] <- "linear"
#' globalFit <- getGlobalAlignment(oswFiles_DIAlignR, ref = "run2", eXp = "run0",
#'  fitType = params[["globalAlignment"]], maxFdrGlobal = 0.05, spanvalue = 0.1)
#' adaptiveRT <- 77.82315 #3.5*getRSE(globalFit, params[["globalAlignment"]])
#' globalFit <- coef(globalFit)
#' getAlignedTimesFast(XICs.ref, XICs.eXp, globalFit, adaptiveRT, params)
#' @export
getAlignedTimesFast <- function(XICs.ref, XICs.eXp, globalFit, adaptiveRT, params){
  B1p <- getPredict(globalFit, XICs.ref[[1]][1,1], params[["globalAlignment"]])
  len <- nrow(XICs.ref[[1]])
  B2p <- getPredict(globalFit, XICs.ref[[1]][len,1], params[["globalAlignment"]])
  if(is.na(B1p) || is.na(B2p)){
    B1p <- XICs.eXp[[1]][1,1]
    B2p <- XICs.eXp[[1]][nrow(XICs.eXp[[1]]),1]
  }
  tAligned <- getAlignedTimesCpp(XICs.ref, XICs.eXp, params[["kernelLen"]], params[["polyOrd"]], params[["alignType"]],
                     adaptiveRT, params[["normalization"]], params[["simMeasure"]],
                     B1p = B1p, B2p = B2p, params[["goFactor"]], params[["geFactor"]], params[["cosAngleThresh"]],
                     params[["OverlapAlignment"]], params[["dotProdThresh"]], params[["gapQuantile"]], 9L,
                     params[["hardConstrain"]], params[["samples4gradient"]])
  tAligned
}

#' Get aligned indices.
#'
#' This function aligns XICs of reference and experiment runs.
#' It produces aligned indices between refernce run and experiment run.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-06-07
#' @param XICs.ref List of extracted ion chromatograms from reference run.
#' @param XICs.eXp List of extracted ion chromatograms from experiment run.
#' @param globalFit Linear or loess fit object between reference and experiment run.
#' @param alignType Available alignment methods are "global", "local" and "hybrid".
#' @param adaptiveRT (numeric) Similarity matrix is not penalized within adaptive RT.
#' @param normalization (character) Must be selected from "mean", "l2".
#' @param simMeasure (string) Must be selected from dotProduct, cosineAngle, crossCorrelation,
#' cosine2Angle, dotProductMasked, euclideanDist, covariance and correlation.
#' @param goFactor (numeric) Penalty for introducing first gap in alignment. This value is multiplied by base gap-penalty.
#' @param geFactor (numeric) Penalty for introducing subsequent gaps in alignment. This value is multiplied by base gap-penalty.
#' @param cosAngleThresh (numeric) In simType = dotProductMasked mode, angular similarity should be higher than cosAngleThresh otherwise similarity is forced to zero.
#' @param OverlapAlignment (logical) An input for alignment with free end-gaps. False: Global alignment, True: overlap alignment.
#' @param dotProdThresh (numeric) In simType = dotProductMasked mode, values in similarity matrix higher than dotProdThresh quantile are checked for angular similarity.
#' @param gapQuantile (numeric) Must be between 0 and 1. This is used to calculate base gap-penalty from similarity distribution.
#' @param kerLen (integer) In simType = crossCorrelation, length of the kernel used to sum similarity score. Must be an odd number.
#' @param hardConstrain (logical) If FALSE; indices farther from noBeef distance are filled with distance from linear fit line.
#' @param samples4gradient (numeric) This parameter modulates penalization of masked indices.
#' @param objType (char) Must be selected from light, medium and heavy.
#' @return (data-frame) Aligned indices of reference and experiment runs. Gaps are introduced as NA.
#' @seealso \code{\link{alignChromatogramsCpp}, \link{getAlignObj}}
#' @keywords internal
#' @examples
#' data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
#' data(oswFiles_DIAlignR, package="DIAlignR")
#' run1 <- "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"
#' run2 <- "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"
#' XICs.ref <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[[run1]][["4618"]]
#' XICs.eXp <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[[run2]][["4618"]]
#' globalFit <- getGlobalAlignment(oswFiles_DIAlignR, ref = "run1", eXp = "run2",
#'   fitType = "loess", maxFdrGlobal = 0.05, spanvalue = 0.1)
#' adaptiveRT <- 77.82315 #3.5*globalFit$s
#' \dontrun{
#' getAlignedIndices(XICs.ref, XICs.eXp, globalFit, alignType = "hybrid",
#'  adaptiveRT = adaptiveRT, normalization = "mean",
#'   simMeasure = "dotProductMasked", goFactor = 0.125, geFactor = 40, cosAngleThresh = 0.3,
#'   OverlapAlignment = TRUE, dotProdThresh = 0.96, gapQuantile = 0.5, kerLen = 9L, hardConstrain = FALSE,
#'   samples4gradient = 100)
#'  }
getAlignedIndices <- function(XICs.ref, XICs.eXp, globalFit, alignType, adaptiveRT,
                            normalization, simMeasure, goFactor, geFactor, cosAngleThresh,
                            OverlapAlignment, dotProdThresh, gapQuantile, kerLen, hardConstrain,
                            samples4gradient, objType = "light"){
  alignObj <- getAlignObj(XICs.ref, XICs.eXp, globalFit, alignType, adaptiveRT,
                          normalization, simType = simMeasure, goFactor, geFactor,
                          cosAngleThresh, OverlapAlignment, dotProdThresh, gapQuantile,
                          kerLen, hardConstrain, samples4gradient, objType)
  alignedIndices <- cbind(alignObj@indexA_aligned,
                          alignObj@indexB_aligned,
                          alignObj@score)
  colnames(alignedIndices) <- c("indexAligned.ref", "indexAligned.eXp", "score")
  alignedIndices[, 1:2][alignedIndices[, 1:2] == 0] <- NA_integer_
  alignedIndices
}


alignedTimes2 <- function(alignObj, XICs.ref, XICs.eXp){
  alignedIndices <- cbind(alignObj@indexA_aligned,
                          alignObj@indexB_aligned,
                          alignObj@score)
  colnames(alignedIndices) <- c("indexAligned.ref", "indexAligned.eXp", "score")
  alignedIndices[, 1:2][alignedIndices[, 1:2] == 0] <- NA_integer_
  keep <- !is.na(alignedIndices[,"indexAligned.ref"])
  tVec.ref <- XICs.ref[[1]][["time"]] # Extracting time component
  tVec.eXp <- XICs.eXp[[1]][["time"]] # Extracting time component
  tAligned.ref <- mapIdxToTime(tVec.ref, alignedIndices[,"indexAligned.ref"])
  tAligned.eXp <- mapIdxToTime(tVec.eXp, alignedIndices[,"indexAligned.eXp"])
  tAligned.ref <- tAligned.ref[keep]
  tAligned.eXp <- tAligned.eXp[keep]
  data.frame(tAligned.ref, tAligned.eXp)
}
