#' Get child chromatograms from parents
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-05-23
#' @param XICs.ref (data-frame) extracted ion chromatograms from reference run.
#' @param XICs.eXp (data-frame) extracted ion chromatograms from experiment run.
#' @param alignedIndices (data-frame) must have two columns "indexAligned.ref" and "indexAligned.eXp".
#' @param method (string) must be either "spline", "sgolay" or "linear".
#' @param polyOrd (integer) must be less than kernelLen.
#' @param kernelLen (integer) must be an odd integer.
#' @param splineMethod (string) must be either "fmm" or "natural".
#' @param keepFlanks (logical) TRUE: Flanking chromatogram is not removed.
#' @param mergeStrategy (string) must be either ref, avg, refStart or refEnd.
#' @return (list) a list of chromatograms.
#' @seealso \code{\link{alignChromatogramsCpp}, \link{getAlignObj}}
#' @examples
#' data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
#' data(alignObj_DIAlignR, package="DIAlignR")
#' XICs.ref <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["run1"]][["14299_QFNNTDIVLLEDFQK/3"]]
#' XICs.eXp <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["run2"]][["14299_QFNNTDIVLLEDFQK/3"]]
#' alignedIndices <- cbind(alignObj_DIAlignR@indexA_aligned, alignObj_DIAlignR@indexB_aligned)
#' colnames(alignedIndices) <- c("indexAligned.ref", "indexAligned.eXp")
#' alignedIndices[, 1:2][alignedIndices[, 1:2] == 0] <- NA_integer_
#' newXICs <- childXICs(XICs.ref, XICs.eXp, alignedIndices)
#' plotXICgroup(XICs.ref)
#' plotXICgroup(newXICs)
#' @export
childXICs <- function(XICs.ref, XICs.eXp, alignedIndices, method = "spline", polyOrd = 4,
                     kernelLen = 9, splineMethod = "fmm", mergeStrategy = "avg", keepFlanks = FALSE){
  child_xics <- lapply(seq_along(XICs.ref), function(i){
    childXIC(XICs.ref[[i]], XICs.eXp[[i]], alignedIndices, method, polyOrd, kernelLen,
             splineMethod, mergeStrategy, keepFlanks)
  })
  for(i in seq_along(child_xics)){
    colnames(child_xics[[i]]) <- c("time", paste0("intensity", i))
  }
  child_xics
  }

#' Get child chromatogram from parents
#'
#' Three startegy to pick time to keep sampling rate equal to parents:
#'  - Use reference time
#'  - Start with the average then fixed sampling rate
#'  - Average time (Non-gap region will have parent's sampling rate)
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-05-23
#' @param XICs.ref (data-frame) extracted ion chromatograms from reference run.
#' @param XICs.eXp (data-frame) extracted ion chromatograms from experiment run.
#' @param alignedIndices (data-frame) must have two columns "indexAligned.ref" and "indexAligned.eXp".
#' @param method (string) must be either "spline", "sgolay" or "linear".
#' @param polyOrd (integer) must be less than kernelLen.
#' @param kernelLen (integer) must be an odd integer.
#' @param splineMethod (string) must be either "fmm" or "natural".
#' @param keepFlanks (logical) TRUE: Flanking chromatogram is not removed.
#' @param mergeStrategy (string) must be either ref, avg, refStart or refEnd.
#' @return (data-frame) a chromatogram.
#' @seealso \code{\link{mergeXIC}, \link{alignedXIC}}
#' @keywords internal
#' @examples
#' data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
#' data(alignObj_DIAlignR, package="DIAlignR")
#' XICs.ref <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["run1"]][["14299_QFNNTDIVLLEDFQK/3"]]
#' XICs.eXp <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["run2"]][["14299_QFNNTDIVLLEDFQK/3"]]
#' alignedIndices <- cbind(alignObj_DIAlignR@indexA_aligned, alignObj_DIAlignR@indexB_aligned)
#' colnames(alignedIndices) <- c("indexAligned.ref", "indexAligned.eXp")
#' alignedIndices[, 1:2][alignedIndices[, 1:2] == 0] <- NA_integer_
#' \dontrun{
#' childXIC(XICs.ref[[1]], XICs.eXp[[1]], alignedIndices)
#' }
childXIC <- function(XIC.ref, XIC.eXp, alignedIndices, method = "spline", polyOrd = 4,
                      kernelLen = 9, splineMethod = "fmm", mergeStrategy = "avg", keepFlanks = TRUE){
  XIC.ref.imp <- alignedXIC(XIC.ref, alignedIndices[,"indexAligned.ref"], method, polyOrd, kernelLen, splineMethod)
  XIC.eXp.imp <- alignedXIC(XIC.eXp, alignedIndices[,"indexAligned.eXp"], method, polyOrd, kernelLen, splineMethod)
  flank <- is.na(XIC.ref.imp[["time"]]) | is.na(XIC.eXp.imp[["time"]])
  skip <- flank | is.na(alignedIndices[,"indexAligned.ref"]) # Remove gaps from the reference

  newXIC <- mergeXIC(XIC.ref.imp[!skip,], XIC.eXp.imp[!skip,], mergeStrategy)

  if(keepFlanks){
    refFlank <- flank & !is.na(alignedIndices[,"indexAligned.ref"])
    eXpFlank <- flank & !is.na(alignedIndices[,"indexAligned.eXp"])

    # Add flanking sequence to the left of a chromatogram
    if(refFlank[1]){
      newXIC <- addFlankToLeft(refFlank, XIC.ref.imp, newXIC)
    } else if(eXpFlank[1]){
      newXIC <- addFlankToLeft(eXpFlank, XIC.eXp.imp, newXIC)
    }

    # Add flanking sequence to the right of a chromatogram
    np <- length(flank)
    if(refFlank[np]){
      newXIC <- addFlankToRight(refFlank, XIC.ref.imp, newXIC)
    } else if(eXpFlank[np]){
      newXIC <- addFlankToRight(eXpFlank, XIC.eXp.imp, newXIC)
    }
  }
  # check if this is a good child either by alignment or some other methods
  newXIC
}

addFlankToLeft <- function(flankSeq, XIC.imp, newXIC){
  flankStop <- min(which(!flankSeq))-1
  idx <- 1:flankStop
  flankXIC <- XIC.imp[idx,]
  t1 <- XIC.imp[flankStop+1, "time"]
  t1_new <- newXIC[1,"time"]
  flankXIC[,"time"] <- flankXIC[,"time"] - t1 + t1_new
  names(flankXIC) <- c("time", "intensity")
  rbind(flankXIC, newXIC) # Add flanking sequence to the left of a chromatogram
}

addFlankToRight <- function(flankSeq, XIC.imp, newXIC){
  np <- length(flankSeq)
  flankStart <- max(which(!flankSeq))+1
  idx <- flankStart:np
  flankXIC <- XIC.imp[idx,]
  tn <- XIC.imp[flankStart-1, "time"]
  tn_new <- newXIC[nrow(newXIC),"time"]
  flankXIC[,"time"] <- flankXIC[,"time"] - tn + tn_new
  names(flankXIC) <- c("time", "intensity")
  rbind(newXIC, flankXIC) # Add flanking sequence to the right of a chromatogram
}

#' Merge two XICs into one
#'
#' Both extracted-ion chromatograms must have same number of rows.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-05-23
#' @param XICs.ref (data-frame) extracted ion chromatograms from reference run.
#' @param XICs.eXp (data-frame) extracted ion chromatograms from experiment run.
#' @param mergeStrategy (string) must be either ref, avg, refStart or refEnd.
#' @return (data-frame) a chromatogram.
#' @seealso \code{\link{alignChromatogramsCpp}, \link{getAlignObj}}
#' @keywords internal
#' @examples
#' data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
#' XICs.ref <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["run1"]][["14299_QFNNTDIVLLEDFQK/3"]]
#' XICs.eXp <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["run2"]][["14299_QFNNTDIVLLEDFQK/3"]]
#' \dontrun{
#' mergeXIC(XICs.ref[[1]], XICs.eXp[[1]], mergeStrategy="ref")
#' }
mergeXIC <- function(XIC.ref, XIC.eXp, mergeStrategy){
  # Recalculate time and intensity for merged chromatogram
  np <- nrow(XIC.ref)
  if(mergeStrategy == "ref"){
    time <- XIC.ref[["time"]]
  } else if (mergeStrategy == "avg") {
    time <- (XIC.ref[, 1] + XIC.eXp[, 1])/2
  } else if (mergeStrategy == "refStart"){
    t1 <- (XIC.ref[["time"]][1] + XIC.eXp[["time"]][1])/2
    delta <- (XIC.ref[["time"]][np] - XIC.ref[["time"]][1])/(np -1)
    time <- seq(from = t1, to = t1+delta*(np-1), by = delta)
  } else if (mergeStrategy == "refEnd"){
    tn <- (XIC.ref[["time"]][np] + XIC.eXp[["time"]][np])/2
    delta <- (XIC.ref[["time"]][np] - XIC.ref[["time"]][1])/(np -1)
    time <- seq(from = tn-delta*(np-1), to = tn, by = delta)
  } else {
    stop("mergeStrategy is not correct.")
  }
  # Caveat: If chromatograms are not perfectly aligned, we may widen the peak.
  # TODO: Put weights on the intensity. Align child chromatogram to parent ones and see the fit.
  # If the fit is not good propagate the reference one.
  intensity <- (XIC.ref[, 2] + XIC.eXp[, 2])/2
  data.frame(time, intensity)
}

#' Create an aligned chromatogram
#'
#' Modifies chromatogram to have the same length as indices. Imputes missing values with appropriate method.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-05-23
#' @param XIC (data-frame) first column is time, second column is intensity.
#' @param indices (integer) vector of monotonically increasing integers.
#' @param method (string) must be either "spline", "sgolay" or "linear".
#' @param polyOrd (integer) must be less than kernelLen.
#' @param kernelLen (integer) must be an odd integer.
#' @param splineMethod (string) must be either "fmm" or "natural".
#' @return (data-frame) a chromatogram.
#' @seealso \code{\link{childXICs}, \link{imputeChromatogram}}
#' @keywords internal
#' @examples
#' data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
#' data(alignObj_DIAlignR, package="DIAlignR")
#' XICs.ref <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["run1"]][["14299_QFNNTDIVLLEDFQK/3"]]
#' alignedIndices <- cbind(alignObj_DIAlignR@indexA_aligned, alignObj_DIAlignR@indexB_aligned)
#' colnames(alignedIndices) <- c("indexAligned.ref", "indexAligned.eXp")
#' alignedIndices[, 1:2][alignedIndices[, 1:2] == 0] <- NA_integer_
#' \dontrun{
#' alignedXIC(XICs.ref[[1]], alignedIndices[,"indexAligned.ref"])
#' }
alignedXIC <- function(XIC, indices, method = "spline", polyOrd = 4, kernelLen = 9, splineMethod = "fmm"){
  XIC.imp <- XIC[indices, ] # Create XIC of the same length as of indices
  rownames(XIC.imp) <- NULL
  XIC.imp[["time"]] <- mapIdxToTime(XIC[["time"]], indices) # Interpolate time
  skip <- is.na(XIC.imp[["time"]]) # Find flanking missing values
  XIC.imp2 <- XIC.imp[!skip, ] # Remove flanking missing values as imputation doesn't support them
  XIC.imp2 <- imputeChromatogram(XIC.imp2, method, polyOrd, kernelLen, splineMethod)
  XIC.imp[!skip, ] <- XIC.imp2
  XIC.imp
}
