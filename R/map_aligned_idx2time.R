#' Establishes mapping back from index to time
#'
#' Takes a time vector and index vector of same length. This function create a
#' new time vector given indices specified in \code{idx}. For \code{NA} indices
#' it uses last index to fill time value. For \code{NA} appearing at the start
#' of \code{idx}, it uses next index to get time value.
#' @param timeVec A numeric vector
#' @param idx An integer vector
#' @return A mutated time vector
#' @importFrom zoo na.locf
#' @export
mapIdxToTime <- function(timeVec, idx){
    if (!requireNamespace("zoo", quietly = TRUE)) {
        stop("Package \"zoo\" needed for this function to work. Please install it.",
             call. = FALSE)
    }
    mutateT <- na.locf(na.locf(timeVec[idx], na.rm = FALSE), fromLast = TRUE)
    return(mutateT)
}
