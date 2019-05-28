#' Establishes pairwise correspondence between extracted-ion chromatograms.
#'
#' This function calculates three matrices for affine gap alignment using input
#' similarity matrix and affine gap opening and gap closing penalties. An
#' implementation of Needleman-Wunsch alignment and overlap alignment is also
#' provided. All three matrices are clubbed together in an output S4 object.
#' @param featureTable A numeric matrix
#' @param pairName A numeric value
#' @param XICs.A A numeric value
#' @param XICs.B A logical scalar
#' @param alignType A logical scalar
#' @param oswFeatureList A logical scalar
#' @param expRSE A logical scalar
#' @param samplingTime A logical scalar
#' @param samples4gradient A logical scalar
#' @param RSEdistFactor A logical scalar
#' @param hardConstrain A logical scalar
#' @param gapQuantile A logical scalar
#' @param goFactor A logical scalar
#' @param geFactor A logical scalar
#' @param normalization A logical scalar
#' @param simMeasure A logical scalar
#' @return A numeric matrix
#' @examples
#' pair <- "run1_run2"
#' MappedTime2 <- getPepPeakCorp(StrepAnnot, pair, StrepChroms[["run1"]], StrepChroms[["run2"]],
#' "hybrid", oswOutStrep, expRSE = 7.4)
#' @export
getPepPeakCorp <- function(featureTable, pairName, XICs.A, XICs.B, alignType, oswFeatureList, expRSE,
                           samplingTime = 3.4, samples4gradient = 100, RSEdistFactor = 3.5,
                           hardConstrain = FALSE, gapQuantile = 0.5, goFactor = 0.125, geFactor = 40,
                           normalization = "mean", simMeasure = "dotProductMasked"){
  peptides <- rownames(featureTable)
  run_pair <- strsplit(pairName, split = "_")[[1]]
  # In finding global fit, peptides are removed while training the fit.
  Loess.fit <- getLOESSfit(run_pair, peptides, oswFeatureList, 0.1)
  rse <- Loess.fit$s # Residual Standard Error
  MappedTime <- matrix(NA, nrow = length(peptides), ncol = 1)
  rownames(MappedTime) <- peptides
  for(peptide in peptides){
    intensityListA <- lapply(XICs.A[[peptide]], `[[`, 2) # Extracting intensity values
    intensityListB <- lapply(XICs.B[[peptide]], `[[`, 2) # Extracting intensity values
    tAVec <- XICs.A[[peptide]][[1]][["time"]] # Extracting time component
    tBVec <- XICs.B[[peptide]][[1]][["time"]] # Extracting time component
    noBeef <- ceiling(RSEdistFactor*min(rse, expRSE)/samplingTime)
    B1p <- predict(Loess.fit, tAVec[1])
    B2p <- predict(Loess.fit, tAVec[length(tAVec)])
    Alignobj <- alignChromatograms_cpp(intensityListA, intensityListB,
                                       tAVec, tBVec,
                                       normalization, simMeasure,
                                       B1p, B2p, noBeef)
    AlignedIndices <- cbind(Alignobj@indexA_aligned, Alignobj@indexB_aligned, Alignobj@score)
    colnames(AlignedIndices) <- c("indexA_aligned", "indexB_aligned", "score")
    AlignedIndices[, 1:2][AlignedIndices[, 1:2] == 0] <- NA
    tA.aligned <- mapIdxToTime(tAVec, AlignedIndices[,"indexA_aligned"])
    tB.aligned <- mapIdxToTime(tBVec, AlignedIndices[,"indexB_aligned"])
    MappedTime[peptide, 1] <- tB.aligned[which.min(abs(tA.aligned - featureTable[peptide,run_pair[1]]))]
  }
  return(MappedTime)
}
