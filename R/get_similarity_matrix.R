OuterProdEuclFunc <- function(data, pep, runA, runB){
    num_of_frag <- length(data[[runA]][[pep]])
    num_of_samplesA <- length(data[[runA]][[pep]][[1]][,1])
    num_of_samplesB <- length(data[[runB]][[pep]][[1]][,1])
    MeanNormA <- sapply(data[[runA]][[pep]], function(x) sum(x[,2])/num_of_samplesA)
    MeanNormA <- mean(MeanNormA)
    MeanNormB <- sapply(data[[runB]][[pep]], function(x) sum(x[,2])/num_of_samplesB)
    MeanNormB <- mean(MeanNormB)
    outerProdList <- list()
    for (i in 1:num_of_frag){
        NormIntensityA <- data[[runA]][[pep]][[i]][,2]/MeanNormA
        NormIntensityB <- data[[runB]][[pep]][[i]][,2]/MeanNormB
        outerProdList[[i]] <- (outer(NormIntensityA, NormIntensityB, FUN = "-"))**2
    }
    return(outerProdList)}

OuterProdCovFunc <- function(data, pep, runA, runB){
    num_of_frag <- length(data[[runA]][[pep]])
    num_of_samplesA <- length(data[[runA]][[pep]][[1]][,1])
    num_of_samplesB <- length(data[[runB]][[pep]][[1]][,1])
    MeanNormA <- sapply(data[[runA]][[pep]], function(x) sum(x[,2])/num_of_samplesA)
    MeanNormA <- mean(MeanNormA)
    MeanNormB <- sapply(data[[runB]][[pep]], function(x) sum(x[,2])/num_of_samplesB)
    MeanNormB <- mean(MeanNormB)
    s <- matrix(NA, nrow = num_of_samplesA, ncol = num_of_samplesB)
    matA <- sapply(data[[runA]][[pep]], function(x) x[,2])/MeanNormA
    matB <- sapply(data[[runB]][[pep]], function(x) x[,2])/MeanNormB
    for(a in 1:num_of_samplesA){
        vecA <- matA[a,]
        for(b in 1:num_of_samplesB){
            vecB <- matB[b,]
            s[a,b] <- cov(vecA, vecB)
        }
    }
    return(s) }

OuterProdCorFunc <- function(data, pep, runA, runB){
    num_of_frag <- length(data[[runA]][[pep]])
    num_of_samplesA <- length(data[[runA]][[pep]][[1]][,1])
    num_of_samplesB <- length(data[[runB]][[pep]][[1]][,1])
    MeanNormA <- sapply(data[[runA]][[pep]], function(x) sum(x[,2])/num_of_samplesA)
    MeanNormA <- mean(MeanNormA)
    MeanNormB <- sapply(data[[runB]][[pep]], function(x) sum(x[,2])/num_of_samplesB)
    MeanNormB <- mean(MeanNormB)
    s <- matrix(NA, nrow = num_of_samplesA, ncol = num_of_samplesB)
    matA <- sapply(data[[runA]][[pep]], function(x) x[,2])/MeanNormA
    matB <- sapply(data[[runB]][[pep]], function(x) x[,2])/MeanNormB
    for(a in 1:num_of_samplesA){
        vecA <- matA[a,]
        for(b in 1:num_of_samplesB){
            vecB <- matB[b,]
            s[a,b] <- cor(vecA, vecB)
        }
    }
    return(s) }

#' @export
OuterProdMeanNormAll6Func <- function(data, pep, runA, runB){
    num_of_frag <- length(data[[runA]][[pep]])
    num_of_samplesA <- length(data[[runA]][[pep]][[1]][,1])
    num_of_samplesB <- length(data[[runB]][[pep]][[1]][,1])
    MeanNormA <- sapply(data[[runA]][[pep]], function(x) sum(x[,2])/num_of_samplesA)
    MeanNormA <- mean(MeanNormA)
    MeanNormB <- sapply(data[[runB]][[pep]], function(x) sum(x[,2])/num_of_samplesB)
    MeanNormB <- mean(MeanNormB)
    outerProdList <- list()
    for (i in 1:num_of_frag){
        NormIntensityA <- data[[runA]][[pep]][[i]][,2]/MeanNormA
        NormIntensityB <- data[[runB]][[pep]][[i]][,2]/MeanNormB
        outerProdList[[i]] <- outer(NormIntensityA, NormIntensityB)
    }
    return(outerProdList) }
# OuterProdMeanNormAll6 <- OuterProdMeanNormAll6Func(data, pep, "run1", "run2")

OuterProdL2NormAllFunc <- function(data, pep, runA, runB){
    num_of_frag <- length(data[[runA]][[pep]])
    L2NormA <- sapply(data[[runA]][[pep]], function(x) x[,2])
    L2NormA <- sqrt(rowSums(L2NormA^2))
    L2NormB <- sapply(data[[runB]][[pep]], function(x) x[,2])
    L2NormB <- sqrt(rowSums(L2NormB^2))
    outerProdList <- list()
    for (i in 1:num_of_frag){
        NormIntensityA <- data[[runA]][[pep]][[i]][,2]/L2NormA
        NormIntensityA[is.nan(NormIntensityA)] <-0
        NormIntensityB <- data[[runB]][[pep]][[i]][,2]/L2NormB
        NormIntensityB[is.nan(NormIntensityB)] <-0
        outerProdList[[i]] <- outer(NormIntensityA, NormIntensityB)
    }
    return(outerProdList) }

add <- function(x) Reduce("+", x)

#' calculates similarity matrix between two chromatogram groups
#'
#' This function takes in datafile, precursor id, names of run pair. Based on
#' similarity measure, it calculates the similarity matrix.
#' @param data Chromatographic data
#' @param pep precursor ID
#' @param runA First run of the run-pair
#' @param runB Second run of the run-pair
#' @param type A character string. must be selected from (dotProduct,
#'   cosineAngle, cosine2Angle, dotProductMasked, euclideanDist, covariance and
#'   correlation)
#' @export
getSimilarityMatrix <- function(data, pep, runA, runB, type = c("dotProductMasked", "dotProduct", "cosineAngle", "cosine2Angle", "euclideanDist", "covariance", "correlation"), dotProdThresh = 0.96, cosAngleThresh = 0.3){
    type <- match.arg(type)
    switch(type,
           dotProduct = {OuterProdNormAll6 <- OuterProdMeanNormAll6Func(data, pep, runA, runB); s <- add(OuterProdNormAll6)},
           cosineAngle = {OuterProdL2NormAll <- OuterProdL2NormAllFunc(data, pep, runA, runB); s <- add(OuterProdL2NormAll)},
           cosine2Angle = {OuterProdL2NormAll <- OuterProdL2NormAllFunc(data, pep, runA, runB); s <- cos(2*acos(pmin(add(OuterProdL2NormAll), 1)))},
           dotProductMasked = {OuterProdNormAll6 <- OuterProdMeanNormAll6Func(data, pep, runA, runB)
           s1 <- add(OuterProdNormAll6)
           OuterProdL2NormAll <- OuterProdL2NormAllFunc(data, pep, runA, runB)
           s2 <- cos(2*acos(pmin(add(OuterProdL2NormAll), 1)))
           MASK <- (s1 > quantile(s1, dotProdThresh))
           AngleGreat <- (((1*MASK)*s2) + (1-MASK)) > cosAngleThresh
           s <- s1*(1*AngleGreat)
           },
           euclideanDist = {OuterProdEucl <- OuterProdEuclFunc(data, pep, runA, runB); s <- 1/(1 + sqrt(add(OuterProdEucl)))},
           covariance = {s <- OuterProdCovFunc(data, pep, runA, runB)},
           correlation = {s <- OuterProdCorFunc(data, pep, runA, runB); s[is.na(s)] <- 0} )
    return(s)
}
