#' Output matrix with alignment path represented as binary-hot encoding.
#'
#' This function outputs matix of the same size of input similarity matrix s.
#' The output matrix has unit cell value for row and column indices given in
#' the matrix alignedIndicesObj.
#' @param alignedIndicesObj A numeric matrix
#' @param s A numeric matrix
#' @return A numeric matrix
#' @export
getAlignmentPath <- function(alignedIndicesObj, s){
    path <- matrix(0, nrow = nrow(s), ncol = ncol(s))
    removeRowIndices <- c()
    i <- 1
    while(is.na(alignedIndicesObj[i,"indexA_aligned"]) | is.na(alignedIndicesObj[i,"indexB_aligned"])){
        removeRowIndices <- c(removeRowIndices, i)
        i <- i+1
    }
    i <- nrow(alignedIndicesObj)
    while(is.na(alignedIndicesObj[i,"indexA_aligned"]) | is.na(alignedIndicesObj[i,"indexB_aligned"])){
        removeRowIndices <- c(removeRowIndices, i)
        i <- i-1
    }
    alignedIndicesObj <- alignedIndicesObj[-c(removeRowIndices), ]
    alignedIndicesObj[,"indexA_aligned"] <- na.locf(alignedIndicesObj[,"indexA_aligned"])
    alignedIndicesObj[,"indexB_aligned"] <- na.locf(alignedIndicesObj[,"indexB_aligned"])
    for(i in 1:nrow(alignedIndicesObj)) path[alignedIndicesObj[i,"indexA_aligned"], alignedIndicesObj[i,"indexB_aligned"]] <- 1
    return(path)
}
