#' Computes indices of row and column of maximum scored traceback path
#'
#' This function takes aligned object of the class NeedleObjAffine as an input.
#' With the provided traceback path, it calculates row and column index pair
#' associated with the highest scoring path through similarity matrix. Output is
#' a list of row-column index pairs of all highest scoring traceback paths.
#' @param alignObj An object of class NeedleObjAffine
#' @return A list of matrices having row and column indices of traceback path.
#' @export
getAlignment <- function(alignObj){
    NumRow <- nrow(alignObj@M); NumCol <- ncol(alignObj@M)
    if(!alignObj@FreeEndGaps){
        startRow <- NumRow
        startCol <- NumCol
        AlignmentScore <- max(alignObj@M[startRow, startCol], alignObj@A[startRow, startCol], alignObj@B[startRow, startCol])
        matrixName <- c("TrM", "TrA", "TrB")[c(alignObj@M[startRow, startCol], alignObj@A[startRow, startCol], alignObj@B[startRow, startCol]) == AlignmentScore]
    }else{
        AlignmentScore = -Inf
        for(i in 1:NumRow){
            score <- max(alignObj@M[i, NumCol], alignObj@A[i, NumCol], alignObj@B[i, NumCol])
            if(score > AlignmentScore){
                startRow <- i; startCol <- NumCol; AlignmentScore <- score
                matrixName <- c("TrM", "TrA", "TrB")[c(alignObj@M[startRow, startCol], alignObj@A[startRow, startCol], alignObj@B[startRow, startCol]) == AlignmentScore]
            }
        }
        for(j in 1:NumCol){
            score <- max(alignObj@M[NumRow, j], alignObj@A[NumRow, j], alignObj@B[NumRow, j])
            if(score > AlignmentScore){
                startRow <- NumRow; startCol <- j; AlignmentScore <- score
                matrixName <- c("TrM", "TrA", "TrB")[c(alignObj@M[startRow, startCol], alignObj@A[startRow, startCol], alignObj@B[startRow, startCol]) == AlignmentScore]
            }
        }
    }
    # startRow; startCol; AlignmentScore; matrixName

    AlignedIndices <- list()
    for(matName in matrixName){
        indexA_aligned <- vector(); indexB_aligned <- vector(); score <- vector()
        if(startCol != NumCol){
            indexA_aligned <- c(rep(NA, NumCol-startCol), indexA_aligned)
            indexB_aligned <- c(startCol:(NumCol-1), indexB_aligned)
            score <- c(rep(AlignmentScore, NumCol-startCol), score)
        } else if(startRow!= NumRow){
            indexA_aligned <- c(startRow:(NumRow-1), indexA_aligned)
            indexB_aligned <- c(rep(NA, NumRow-startRow), indexB_aligned)
            score <- c(rep(AlignmentScore, NumRow-startRow), score)
        }
        i <- startRow; j <- startCol
        pointer <- alignObj@Traceback[[matName]][i,j]
        score <- c(slot(alignObj, substr(matName, start = 3, stop = 3))[i,j], score)
        while(!is.na(pointer)){
            MAB <- substr(pointer, start = 2, stop = 2)
            # D: Diagonal, T: Top, L: Left
            if(grepl("D", pointer)) {
                i <- i-1; j <-j-1
                indexA_aligned <- c(i, indexA_aligned)
                indexB_aligned <- c(j, indexB_aligned)
                score <- c(slot(alignObj, MAB)[i,j], score)
            } else if(grepl("T", pointer)) {
                i <- i-1
                indexA_aligned <- c(i, indexA_aligned)
                indexB_aligned <- c(NA, indexB_aligned)
                score <- c(slot(alignObj, MAB)[i,j], score)
            } else {
                j <- j-1
                indexA_aligned <- c(NA, indexA_aligned)
                indexB_aligned <- c(j, indexB_aligned)
                score <- c(slot(alignObj, MAB)[i,j], score)
            }
            pointer <- alignObj@Traceback[[paste0("Tr",MAB)]][i,j]
        }
        score <- score[-c(1)]
        AlignedIndices[[matName]] <- cbind(indexA_aligned, indexB_aligned, score)
    }
    return(AlignedIndices)
}

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
