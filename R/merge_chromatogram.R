#' Get child chromatograms from parents
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-05-23
#' @inheritParams childXIC
#' @param XICs.ref (list of data-frames) extracted ion chromatograms from reference run.
#' @param XICs.eXp (list of data-frames) extracted ion chromatograms from experiment run.
#' @return  (list) the first element is a list of chromatograms. The second element is aligned parent time-vectors.
#' @seealso \code{\link{childXIC}, \link{mergeXIC}}
#' @examples
#' data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
#' data(alignObj_DIAlignR, package="DIAlignR")
#' XICs.ref <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]]
#' XICs.eXp <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]]
#' alignedIndices <- cbind(alignObj_DIAlignR@indexA_aligned, alignObj_DIAlignR@indexB_aligned)
#' colnames(alignedIndices) <- c("indexAligned.ref", "indexAligned.eXp")
#' alignedIndices[, 1:2][alignedIndices[, 1:2] == 0] <- NA_integer_
#' newXICs <- childXICs(XICs.ref, XICs.eXp, alignedIndices)[[1]]
#' plotXICgroup(XICs.ref)
#' plotXICgroup(newXICs)
#' @export
childXICs <- function(XICs.ref, XICs.eXp, alignedIndices, method = "spline", polyOrd = 4,
                     kernelLen = 9, splineMethod = "fmm", wRef = 0.5, mergeStrategy = "avg", keepFlanks = TRUE){
  child_xics <- vector(mode = "list", length = length(XICs.ref))
  alignedVec <- c()
  for(i in seq_along(child_xics)){
    output <- childXIC(XICs.ref[[i]], XICs.eXp[[i]], alignedIndices, method, polyOrd, kernelLen,
             splineMethod, wRef, mergeStrategy, keepFlanks)
    xic <- output[[1]]
    alignedVec <- output[[2]]
    colnames(xic) <- c("time", paste0("intensity", i))
    child_xics[[i]] <- xic
  }
  x <- alignedVec[, "alignedChildTime"]
  if(!all(x == cummax(v = cummax(ifelse(is.na(x), -Inf, x))), na.rm = TRUE)){
    warning("Child run does not have monotonically increasing time.")
  }
  # Interpolate middle values.
  alignedVec[,3:5] <- sapply(3:5, function(i) zoo::na.approx(alignedVec[,i], na.rm = FALSE))
  list(child_xics, alignedVec)
  }


#' Get a child chromatogram from parents
#'
#' Internal gaps in reference chromatogram and corresponding indices in experiment chromatogram
#' are discarded while merging both chromatogram to obtain a child chromatogram.
#'
#' @inheritSection mergeXIC time-merging
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-05-23
#' @inherit alignedXIC return
#' @inheritParams alignedXIC
#' @inheritParams mergeXIC
#' @param alignedIndices (data-frame) must have two columns "indexAligned.ref" and "indexAligned.eXp".
#' @param keepFlanks (logical) TRUE: Flanking chromatogram is not removed.
#' @return (list) the first element is chromatogram. The second element is aligned parent time-vectors.
#' @seealso \code{\link{mergeXIC}, \link{alignedXIC}, \link{childXICs}}
#' @keywords internal
#' @examples
#' data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
#' data(alignObj_DIAlignR, package="DIAlignR")
#' XICs.ref <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]]
#' XICs.eXp <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]]
#' alignedIndices <- cbind(alignObj_DIAlignR@indexA_aligned, alignObj_DIAlignR@indexB_aligned)
#' colnames(alignedIndices) <- c("indexAligned.ref", "indexAligned.eXp")
#' alignedIndices[, 1:2][alignedIndices[, 1:2] == 0] <- NA_integer_
#' \dontrun{
#' plot(childXIC(XICs.ref[[1]], XICs.eXp[[1]], alignedIndices)[[1]], type = 'l')
#' }
childXIC <- function(XIC.ref, XIC.eXp, alignedIndices, method = "spline", polyOrd = 4,
                      kernelLen = 9, splineMethod = "fmm", wRef = 0.5, mergeStrategy = "avg", keepFlanks = TRUE){
  # Impute missing intensities in chromatogram
  XIC.ref.imp <- alignedXIC(XIC.ref, alignedIndices[,"indexAligned.ref"], method, polyOrd, kernelLen, splineMethod)
  XIC.eXp.imp <- alignedXIC(XIC.eXp, alignedIndices[,"indexAligned.eXp"], method, polyOrd, kernelLen, splineMethod)

  flank <- is.na(XIC.ref.imp[["time"]]) | is.na(XIC.eXp.imp[["time"]])
  skip <- flank | is.na(alignedIndices[,"indexAligned.ref"]) # Remove gaps from the reference

  newXIC <- mergeXIC(XIC.ref.imp[!skip,], XIC.eXp.imp[!skip,], wRef, mergeStrategy)
  alignedChildTime <- rep(NA, nrow(alignedIndices))
  alignedChildTime[!skip] <- newXIC[["time"]]

  if(keepFlanks){
    startTime <- newXIC[["time"]][1]
    endTime <- newXIC[["time"]][nrow(newXIC)]
    refFlank <- flank & !is.na(alignedIndices[,"indexAligned.ref"]) # Flanking signal in reference
    eXpFlank <- flank & !is.na(alignedIndices[,"indexAligned.eXp"]) # Flanking signal in experiment

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

    # Copy the flanking signal to alignedChildTime
    flankStop <- match(startTime, newXIC[["time"]])
    if(flankStop != 1) alignedChildTime[1:flankStop] <- newXIC[["time"]][1:flankStop]

    flankStart <- match(endTime, newXIC[["time"]])
    len <- nrow(newXIC)
    stp <- length(alignedChildTime)
    strt <- stp-(len-flankStart)+1
    if(flankStart != len) alignedChildTime[strt:stp] <- newXIC[["time"]][(flankStart+1):len]
  }

  # Group indices, time and aligned time. This will be used to map features on child XIC.
  tAligned.ref <- XIC.ref[alignedIndices[, "indexAligned.ref"],"time"]
  tAligned.eXp <- XIC.eXp[alignedIndices[, "indexAligned.eXp"],"time"]
  alignedVec <- cbind(alignedIndices, tAligned.ref, tAligned.eXp, alignedChildTime)
  # check if this is a good child either by alignment or some other methods
  list(newXIC, alignedVec)
}


#' Add signal to the left of XIC
#'
#' This function copies flanking chromatogram from XIC and paste it to the left of newXIC.
#' @inherit addFlankToRight
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-07-16
#' @param flankSeq (logical) must be TRUE at the front of the vector.
#' @seealso \code{\link{childXIC}, \link{addFlankToRight}}
#' @keywords internal
#' @examples
#' time <- seq(from = 3003.4, to = 3048, by = 3.4)
#' y <- c(0.2050595, 0.8850070, 2.2068768, 3.7212677, 5.1652605, 5.8288915, 5.5446804,
#'        4.5671360, 3.3213154, 1.9485889, 0.9520709, 0.3294218, 0.2009581, 0.1420923)
#' chrom <- data.frame(time, y)
#' chrom2 <- data.frame(time = c(3013.4, 3016, 3020), intensity = c(1.2, 3.4, 5.6))
#' flankSeq <- as.logical(c(1,1,0,0,0,0,0,0,0,0,0,0,1,1))
#' \dontrun{
#' addFlankToLeft(flankSeq, chrom, chrom2)
#' }
addFlankToLeft <- function(flankSeq, XIC, newXIC){
  flankStop <- min(which(!flankSeq))-1
  idx <- seq(from =1, to=flankStop, by = 1L)
  flankXIC <- XIC[idx,]
  t1 <- XIC[flankStop+1, "time"]
  t1_new <- newXIC[1,"time"]
  flankXIC[,"time"] <- flankXIC[,"time"] - t1 + t1_new
  names(flankXIC) <- c("time", "intensity")
  rbind(flankXIC, newXIC) # Add flanking sequence to the left of a chromatogram
}


#' Add signal to the right of XIC
#'
#' This function copies flanking chromatogram from XIC and paste it to the right of newXIC.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-07-16
#' @param flankSeq (logical) must be TRUE at the tail of the vector.
#' @param XIC (data-frame) first column is time, second column is intensity.
#' @param newXIC (data-frame) first column is time, second column is intensity.
#' @return (dataframe) has two columns:
#' \item{time}{(numeric)}
#' \item{intensity}{(numeric)}
#' @seealso \code{\link{childXIC}, \link{addFlankToLeft}}
#' @keywords internal
#' @examples
#' time <- seq(from = 3003.4, to = 3048, by = 3.4)
#' y <- c(0.2050595, 0.8850070, 2.2068768, 3.7212677, 5.1652605, 5.8288915, 5.5446804,
#'        4.5671360, 3.3213154, 1.9485889, 0.9520709, 0.3294218, 0.2009581, 0.1420923)
#' chrom <- data.frame(time, y)
#' chrom2 <- data.frame(time = c(3013.4, 3016, 3020), intensity = c(1.2, 3.4, 5.6))
#' flankSeq <- as.logical(c(1,1,0,0,0,0,0,0,0,0,0,0,1,1))
#' \dontrun{
#' addFlankToRight(flankSeq, chrom, chrom2)
#' }
addFlankToRight <- function(flankSeq, XIC, newXIC){
  np <- length(flankSeq)
  flankStart <- max(which(!flankSeq))+1
  idx <- seq(from = flankStart, to = np, by = 1L)
  flankXIC <- XIC[idx,]
  rownames(flankXIC) <- NULL
  tn <- XIC[flankStart-1, "time"]
  tn_new <- newXIC[nrow(newXIC),"time"]
  flankXIC[,"time"] <- flankXIC[,"time"] - tn + tn_new
  names(flankXIC) <- c("time", "intensity")
  rbind(newXIC, flankXIC) # Add flanking sequence to the right of a chromatogram
}


#' Merge two XICs into one
#'
#' Both extracted-ion chromatograms must have same number of rows. Missing values are not allowed.
#' Intensities are weighted-averaged.
#' @section time-merging:
#' There are three startegies for calculating time to keep sampling rate equal to parents:\cr
#'  - Use reference time. \cr
#'  - Average time (Non-gap region will have parent's sampling rate). \cr
#'  - Start/end with the average then fixed sampling rate.
#' @inherit alignedXIC return
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-05-23
#' @param XIC.ref (data-frame) extracted ion chromatogram from reference run. Must not contain missing values.
#' @param XIC.eXp (data-frame) extracted ion chromatogram from experiment run. Must not contain missing values.
#' @param wRef (numeric) Weight of the reference XIC. Must be between 0 and 1.
#' @param mergeStrategy (string) must be either ref, avg, refStart or refEnd.
#' @seealso \code{\link{childXIC}, \link{alignedXIC}}
#' @keywords internal
#' @examples
#' data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
#' XICs.ref <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]]
#' XICs.eXp <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"]][["14299_QFNNTDIVLLEDFQK/3"]]
#' \dontrun{
#' plot(mergeXIC(XICs.ref[[1]], XICs.eXp[[1]], wRef = 0.5, mergeStrategy = "ref"), type = "l")
#' }
mergeXIC <- function(XIC.ref, XIC.eXp, wRef, mergeStrategy){
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
    stop("mergeStrategy is not correct. Must be selected from 'ref', 'avg', 'refStart' and 'refEnd'.")
  }
  # Caveat: If chromatograms are not perfectly aligned, we may widen the peak.
  # TODO: Put weights on the intensity. Align child chromatogram to parent ones and see the fit.
  # If the fit is not good propagate the reference one.
  intensity <- wRef*XIC.ref[, 2] + (1.0 - wRef)*XIC.eXp[, 2]
  data.frame(time, intensity)
}


#' Create an aligned chromatogram
#'
#' Modifies chromatogram to have the same length as indices. Imputes missing values with appropriate method.
#' Time and intensity for the flanking missing indices are set as NA.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-05-23
#' @inherit imputeChromatogram
#' @param XIC (data-frame) first column is time, second column is intensity.
#' @param indices (integer) vector of monotonically increasing integers.
#' @seealso \code{\link{childXIC}, \link{imputeChromatogram}}
#' @keywords internal
#' @examples
#' data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
#' data(alignObj_DIAlignR, package="DIAlignR")
#' XICs.ref <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]]
#' alignedIndices <- cbind(alignObj_DIAlignR@indexA_aligned, alignObj_DIAlignR@indexB_aligned)
#' colnames(alignedIndices) <- c("indexAligned.ref", "indexAligned.eXp")
#' alignedIndices[, 1:2][alignedIndices[, 1:2] == 0] <- NA_integer_
#' \dontrun{
#' plot(alignedXIC(XICs.ref[[1]], alignedIndices[,"indexAligned.ref"]), type = "l")
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
