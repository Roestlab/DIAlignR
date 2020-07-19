#' Smooth chromatogram signal
#'
#' Smoothing methods are Savitzky-Golay, Boxcar, Gaussian kernel and LOESS. Savitzky-Golay smoothing
#' is good at preserving peak-shape compared to gaussian and boxcar smoothing. However, it assumes
#' equidistant points that fortunately is the case for DIA data. This requires a quadratic memory
#' to store the fit and slower than other smoothing methods.
#'
#' Gaussian smoothing uses a gaussian function whose bandwidth is scaled by 0.3706505 to have
#' quartiles at +/- 0.25*bandwidth. The point selection cut-off is also hard at 0.3706505*4*bandwidth.
#'
#' qnorm(0.75, sd = 0.3706505)
#'
#' The definition of C_ksmooth can be found using
#' getAnywhere('C_ksmooth')
#' stats:::C_ksmooth
#'
#' @importFrom stats ksmooth predict
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL3
#' Date: 2020-02-21
#' @param chromatogram (dataframe) A dataframe of two columns. First column must always be
#' monotonically increasing.
#' @param type (char) must be either sgolay, boxcar, gaussian, loess or none.
#' @param samplingTime (numeric) Time difference between neighboring points.
#' @param kernelLen (integer) Number of data-points to consider in the kernel.
#' @param polyOrd (integer) Order of the polynomial to be fit in the kernel.
#' @return A dataframe with two columns.
#' @examples
#' data("XIC_QFNNTDIVLLEDFQK_3_DIAlignR")
#' chrom <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt"]][["4618"]][[1]]
#' \dontrun{
#' newChrom <- smoothSingleXIC(chrom, type = "sgolay", samplingTime = 3.42, kernelLen = 9,
#'  polyOrd = 3)
#' }
#' @seealso \url{https://terpconnect.umd.edu/~toh/spectrum/Smoothing.html},
#'  \url{https://rafalab.github.io/dsbook/smoothing.html},
#'  \url{https://github.com/SurajGupta/r-source/blob/master/src/library/stats/src/ksmooth.c}
smoothSingleXIC <- function(chromatogram, type, samplingTime = NULL, kernelLen = NULL, polyOrd = NULL){
  time <- chromatogram[[1]]
  if(type == "sgolay"){
    intensity <- signal::sgolayfilt(chromatogram[[2]], p = polyOrd, n = kernelLen)
  } else if (type == "boxcar"){
    intensity <- ksmooth(time, chromatogram[[2]], kernel = "box",
                         bandwidth = kernelLen*samplingTime, n.points = length(time))
    intensity <- intensity[["y"]]
  } else if (type == "gaussian"){
    intensity <- ksmooth(time, chromatogram[[2]], kernel = "normal",
                         bandwidth = kernelLen*samplingTime, n.points = length(time))
    intensity <- intensity[["y"]]
  } else if (type == "loess") {
    spanvalue <- kernelLen/length(time)
    fit <- suppressWarnings(loess(chromatogram[[2]] ~ time, span = spanvalue, degree = polyOrd))
    intensity <- predict(fit, time)
  } else {
    intensity <- chromatogram[[2]]
  }
  data.frame(time, intensity)
}

#' Smooth chromatogram signals from a list
#'
#' Smoothing methods are Savitzky-Golay, Boxcar, Gaussian kernel and LOESS. Savitzky-Golay smoothing
#' is good at preserving peak-shape compared to gaussian and boxcar smoothing. However, it assumes
#' equidistant points that fortunately is the case for DIA data. This requires a quadratic memory
#' to store the fit and slower than other smoothing methods.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL3
#' Date: 2020-02-21
#' @param XICs (A list) A list of dataframe that consists of two columns. First column must be
#' monotonically increasing.
#' @param type (char) must be either sgolay, boxcar, gaussian, loess or none.
#' @param samplingTime (numeric) Time difference between neighboring points.
#' @param kernelLen (integer) Number of data-points to consider in the kernel.
#' @param polyOrd (integer) Order of the polynomial to be fit in the kernel.
#' @return A list.
#' @examples
#' data("XIC_QFNNTDIVLLEDFQK_3_DIAlignR")
#' XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt"]][["4618"]]
#' \dontrun{
#' newXICs <- smoothXICs(XICs, type = "sgolay", samplingTime = 3.42, kernelLen = 9,
#'  polyOrd = 3)
#' }
#' @seealso \url{https://terpconnect.umd.edu/~toh/spectrum/Smoothing.html},
#'  \url{https://rafalab.github.io/dsbook/smoothing.html}
#' @export
smoothXICs <- function(XICs, type = "none", samplingTime = NULL, kernelLen = NULL, polyOrd = NULL){
  newXICs <- lapply(XICs, smoothSingleXIC, type, samplingTime, kernelLen, polyOrd)
  for(i in seq_along(newXICs)){
    colnames(newXICs[[i]]) <- c("time", paste0("intensity", i))
  }
  newXICs
}



#' Selects a part of chromatograms
#'
#' This function trims chromatograms from the end-points.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL3
#' Date: 2020-04-01
#' @param XICs (A list) A list of dataframe that consists of two columns. First column must be
#' monotonically increasing.
#' @param len (numeric) must be between 0.1 and 1.
#' @return A list.
#' @examples
#' data("XIC_QFNNTDIVLLEDFQK_3_DIAlignR")
#' XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt"]][["4618"]]
#' \dontrun{
#' newXICs <- smoothXICs(XICs, len = 0.5)
#' }
trimXICs <- function(XICs, len = 1){
  newXICs <- XICs
  if (len != 1){
    tp <- length(XICs[[1]][["time"]]) # Get the number of time-points
    leftIdx <- floor((tp/2) - (tp/2)*len)
    rightIdx <- ceiling((tp/2) + (tp/2)*len)
    indices <- seq(leftIdx, rightIdx, by =1)
    newXICs <- lapply(newXICs, function(xic){
      df <- xic[indices,]
      row.names(df) <- NULL
      df
    })
  }
  newXICs
}
