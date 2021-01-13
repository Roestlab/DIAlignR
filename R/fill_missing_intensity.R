#' Fill missing values using Savitzky–Golay
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-05-21
#' @importFrom pracma polyfit
#' @param chrom (data-frame) first column is time (must be equidistant), second column is intensity.
#' @param polyOrd (integer) must be less than kernelLen.
#' @param kernelLen (integer) must be an odd integer.
#' @return (dataframe) has two columns:
#' \item{time}{(numeric)}
#' \item{intensity}{(numeric)}
#' @examples
#' time <- seq(from = 3003.4, to = 3048, by = 3.4)
#' y <- c(0.2050595, 0.8850070, 2.2068768, 3.7212677, 5.1652605, 5.8288915, 5.5446804,
#'  4.5671360, 3.3213154, 1.9485889, 0.9520709, 0.3294218, 0.2009581, 0.1420923)
#' chrom <- data.frame(time, y)
#' chrom$y[c(1,8)] <- NA
#' \dontrun{
#' sgolayFill(chrom, polyOrd = 3, kernelLen = 9)
#' }
#' @seealso \code{\link[pracma]{polyfit}, \link[signal]{sgolayfilt}}
#' @keywords internal
sgolayFill <- function(chrom, polyOrd, kernelLen){
  intensity <- chrom[[2]]
  idx <- which(is.na(intensity))
  M <- (kernelLen -1)/2
  n <- length(intensity)
  for(i in idx){
    if(i-M <= 0){
      # Kernel can't go before the start
      if(i == 1){
        ker1 <- seq(1, 2*M-i+1, by = 1)
      } else{
        ker1 <- c(seq(-i+1, -1, by = 1), seq(1, 2*M-i+1, by = 1))
      }
    } else if(i+M > n){
      # Kernel can't go past the end
      if(i == n){
        ker1 <- seq(-2*M+n-i, -1, by = 1)
      } else{
        ker1 <- c(seq(-2*M+n-i, -1, by = 1), seq(1, n-i, by = 1))
      }
    } else {
      # Kernel is inside the vector
      ker1 <- c(seq(-M, -1, by = 1), seq(1, M, by = 1))
    }
    ker2 <- intensity[i + ker1]
    if (any(is.na(ker2))){
      # If any value is NA, polyfit output would be zero.
      intensity[i] <- NA_real_
    } else{
      coef <- polyfit(ker1, ker2, polyOrd)
      intensity[i] <- coef[polyOrd+1] # Fill the missing value
    }
  }
  data.frame(time = chrom[[1]], intensity = intensity)
}

#' Fill missing values using spline
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-05-21
#' @importFrom zoo zoo na.spline coredata
#' @param chrom (data-frames) first column is time, second column is intensity.
#' @param method (string) must be either "fmm" or "natural".
#' @return (dataframe) has two columns:
#' \item{time}{(numeric)}
#' \item{intensity}{(numeric)}
#' @examples
#' time <- seq(from = 3003.4, to = 3048, by = 3.4)
#' y <- c(0.2050595, 0.8850070, 2.2068768, 3.7212677, 5.1652605, 5.8288915, 5.5446804,
#'  4.5671360, 3.3213154, 1.9485889, 0.9520709, 0.3294218, 0.2009581, 0.1420923)
#' chrom <- data.frame(time, y)
#' chrom$y[c(1,8)] <- NA
#' \dontrun{
#' splineFill(chrom)
#' }
#' @seealso \code{\link[zoo]{na.spline}, \link[stats]{spline}}
#' @keywords internal
splineFill <- function(chrom, method = "fmm"){
  time <- chrom[[1]]
  intensity <- chrom[[2]]
  idx <- which(is.na(intensity))
  values <- na.spline(zoo(intensity, time), xout = time[idx], method = method, na.rm = F)
  intensity[idx] <- coredata(values)
  data.frame(time, intensity)
}

#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2021) + GPL-3
#' Date: 2021-01-04
akimaFill <- function(chrom){
  time <- chrom[[1]]
  inten <- chrom[[2]]
  idx <- is.na(inten)
  if(requireNamespace("akima", quietly = TRUE)){
    chrom <- data.frame(time, "intensity" = akima::aspline(time[!idx], inten[!idx], time, method = "improved")[["y"]])
  }  else {
    stop("Please install the akima package to use this functionality")
  }
  return(chrom)
}

#' Fill missing values using linear interpolation
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-05-21
#' @importFrom zoo zoo na.approx na.locf coredata
#' @param chrom (data-frames) first column is time, second column is intensity.
#' @return (dataframe) has two columns:
#' \item{time}{(numeric)}
#' \item{intensity}{(numeric)}
#' @examples
#' time <- seq(from = 3003.4, to = 3048, by = 3.4)
#' y <- c(0.2050595, 0.8850070, 2.2068768, 3.7212677, 5.1652605, 5.8288915, 5.5446804,
#'  4.5671360, 3.3213154, 1.9485889, 0.9520709, 0.3294218, 0.2009581, 0.1420923)
#' chrom <- data.frame(time, y)
#' chrom$y[c(1,8, 14)] <- NA
#' \dontrun{
#' approxFill(chrom)
#' }
#' @seealso \code{\link[zoo]{na.approx}, \link[stats]{approx}}
#' @keywords internal
approxFill <- function(chrom){
  time <- chrom[[1]]
  intensity <- chrom[[2]]
  idx <- which(is.na(intensity))
  values <- na.approx(zoo(intensity, time), xout = time[idx], na.rm = F)
  intensity[idx] <- coredata(values)
  # Extrapolation is not available in R. Therefore use last observation carried forward.
  idx <- which(is.na(intensity))
  values <- na.locf(zoo(intensity, time), xout = time[idx], na.rm = F)
  intensity[idx] <- coredata(values)
  data.frame(time, intensity)
}


#' Fill missing intensities in a chromatogram
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-05-21
#' @importFrom zoo zoo na.approx na.locf coredata
#' @param chromatogram (data-frames) first column is time, second column is intensity.
#' @param method (string) must be either "spline", "sgolay" or "linear".
#' @param polyOrd (integer) must be less than kernelLen.
#' @param kernelLen (integer) must be an odd integer.
#' @param splineMethod (string) must be either "fmm" or "natural".
#' @return (dataframe) has two columns:
#' \item{time}{(numeric)}
#' \item{intensity}{(numeric)}
#' @examples
#' time <- seq(from = 3003.4, to = 3048, by = 3.4)
#' y <- c(0.2050595, 0.8850070, 2.2068768, 3.7212677, 5.1652605, 5.8288915, 5.5446804,
#'  4.5671360, 3.3213154, 1.9485889, 0.9520709, 0.3294218, 0.2009581, 0.1420923)
#' chrom <- data.frame(time, y)
#' chrom$y[c(1,8, 14)] <- NA
#' imputeChromatogram(chrom, "sgolay")
#' imputeChromatogram(chrom, "spline")
#' imputeChromatogram(chrom, "linear")
#' @seealso \code{\link[zoo]{na.approx}, \link[zoo]{na.spline}}
#' @export
imputeChromatogram <- function(chromatogram, method = "spline", polyOrd = 4,
                               kernelLen = 9, splineMethod = "fmm"){
  if(method == "sgolay"){
    chrom <- sgolayFill(chromatogram, polyOrd, kernelLen)
    chrom <- splineFill(chrom, splineMethod)
  } else if(method == "spline"){
    chrom <- splineFill(chromatogram, splineMethod)
  } else {
    chrom <- approxFill(chromatogram)
  }
  chrom
}


# Finding the best startegy. Spline is the best strategy with "fmm" method.
# time <- seq(from = 3003.4, to = 3048, by = 3.4)
# y <- c(0.2050595, 0.8850070, 2.2068768, 3.7212677, 5.1652605, 5.8288915, 5.5446804,
#        4.5671360, 3.3213154, 1.9485889, 0.9520709, 0.3294218, 0.2009581, 0.1420923)
# num <- 1; method <- "sgolay"
# set.seed(1)
# err <- sapply(1:100, function(k){
#   chrom <- data.frame(time, y)
#   s <- seq(1,14)
#   idx <- c()
#   for(n in 1:num){
#     idx <- c(idx, sample(s, 1))
#     s <- setdiff(s, c(idx-1, idx, idx+1))
#   }
#   chrom$y[idx] <- NA_real_
#   x <- imputeChromatogram(chrom, method = method, polyOrd = 4, 9, "fmm")
#   sum(abs(y[idx] - x[idx,2]))
#   })
# sum(err)
