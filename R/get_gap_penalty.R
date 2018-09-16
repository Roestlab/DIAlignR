# Hello, world!
#
# This function calculates gap penalty for dynamic programming.
#

#' Calculates gap penalty
#'
#' Takes a matrix, vectorize it and find the given quantile. This value is given out as gap penalty.
#' @param s A numeric matrix
#' @return The gap penalty for dynamic programming
#' @export
getGapPenalty <- function(s, gapQuantile, type = "dotProduct"){
    switch(type,
           dotProduct = {gapPenalty <- as.numeric(quantile(s, gapQuantile), na.rm = TRUE)},
           cosineAngle = {gapPenalty <- 0.95},
           cosine2Angle = {gapPenalty <- 0.95},
           cosineNangle = {gapPenalty <- 0.95},
           cosineAndDotProd = {gapPenalty <- 1.8},
           euclideanDist = {gapPenalty <- as.numeric(quantile(s, gapQuantile), na.rm = TRUE)},
           covariance = {gapPenalty <- as.numeric(quantile(s, gapQuantile), na.rm = TRUE)},
           correlation = {gapPenalty <- as.numeric(quantile(s, gapQuantile), na.rm = TRUE)},
           dotProductMasked = {gapPenalty <- as.numeric(quantile(s, gapQuantile), na.rm = TRUE)})
    return(gapPenalty)
}
