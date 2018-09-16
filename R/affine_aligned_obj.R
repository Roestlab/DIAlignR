#' An S4 object to group matrices produced by the alignmnet
#'
#' This class combines M, A, B matrices and Traceback list into one object.
#' @param Traceback A list of
#' @param M A numeric matrix
#' @param A A numeric matrix
#' @param B A numeric matrix
#' @param GapOpen A numeric scalar value
#' @param GapExten A numeric scalar value
#' @param FreeEndGaps A logical value
#' @importFrom methods setClass
#' @export
NeedleObjAffine <- setClass(Class="NeedleObjAffine", representation(Traceback = "list", M = "matrix", A = "matrix", B = "matrix", GapOpen = "numeric", GapExten = "numeric", FreeEndGaps = "logical") )

#' Creates affine alignment object
#'
#' This function calculates three matrices for affine gap alignment using input
#' similarity matrix and affine gap opening and gap closing penalties. An
#' implementation of Needleman-Wunsch alignment and overlap alignment is also
#' provided. All three matrices are clubbed together in an output S4 object.
#' @param s A numeric matrix
#' @param go A numeric value
#' @param ge A numeric value
#' @param FreeEndGaps A logical scalar
#' @return A S4 object
#' @examples
#' M = 10; MM = -2; go = 22; ge = 7
#' seq1 = "GCAT"
#' seq2 = "CAGTG"
#' s = matrix(NA, nrow = nchar(seq1), ncol= nchar(seq2))
#' for(i in 1:nrow(s)){
#'    for(j in 1: ncol(s)){
#'        if(substr(seq1, start = i, stop = i) == substr(seq2, start = j, stop = j)) s[i,j] <- M
#'        else s[i,j] <- MM
#'    }
#' }
#' alignObj <- getAffineAlignObj(s, go, ge, FreeEndGaps=FALSE)
#' @export
getAffineAlignObj <- function(s, go, ge, FreeEndGaps=TRUE){
    ChromA_Len <- dim(s)[1]
    ChromB_Len <- dim(s)[2]

    M <- matrix(NA, nrow = ChromA_Len+1, ncol = ChromB_Len+1)
    A <- matrix(NA, nrow = ChromA_Len+1, ncol = ChromB_Len+1)
    B <- matrix(NA, nrow = ChromA_Len+1, ncol = ChromB_Len+1)
    Traceback <- vector("list", 3)
    Traceback[[1]] <- matrix(NA, nrow = ChromA_Len+1, ncol = ChromB_Len+1)
    Traceback[[2]] <- matrix(NA, nrow = ChromA_Len+1, ncol = ChromB_Len+1)
    Traceback[[3]] <- matrix(NA, nrow = ChromA_Len+1, ncol = ChromB_Len+1)
    names(Traceback) <- c("TrM", "TrA", "TrB")

    M[1,] <- -Inf; M[,1] <- -Inf; M[1,1] <- 0
    A[1,] <- -Inf; B[,1] <- -Inf
    Traceback[["TrA"]][2:(ChromA_Len+1), 1] <- "TA"
    Traceback[["TrB"]][1, 2:(ChromB_Len+1)] <- "LB"
    if(FreeEndGaps){
        A[2:(ChromA_Len+1), 1] <- 0
        B[1, 2:(ChromB_Len+1)] <- 0
    } else{
        A[2:(ChromA_Len+1), 1] <- -(0:(ChromA_Len-1))*ge - go
        B[1, 2:(ChromB_Len+1)] <- -(0:(ChromB_Len-1))*ge - go
    }

    for(j in 2:(ChromB_Len+1)){
        for(i in 2:(ChromA_Len+1)){
            Diago <- M[i-1,j-1] + s[i-1, j-1]
            gapInA <- A[i-1,j-1] + s[i-1, j-1]
            gapInB <- B[i-1,j-1] + s[i-1, j-1]
            M[i,j] <- max(Diago, gapInA, gapInB)
            if(M[i,j] == Diago) Traceback[["TrM"]][i, j] <- "DM" # D: Diagonal
            else if(M[i,j] == gapInA) Traceback[["TrM"]][i, j] <- "DA"
            else Traceback[["TrM"]][i, j] <- "DB"

            A[i,j] <- max(M[i-1,j] - go, A[i-1,j] - ge, B[i-1,j] - go)
            if(A[i,j] == M[i-1,j] - go) Traceback[["TrA"]][i,j] <- "TM" # T: Top
            else if(A[i,j] == A[i-1,j] - ge) Traceback[["TrA"]][i,j] <- "TA"
            else Traceback[["TrA"]][i,j] <- "TB"

            B[i,j] <- max(M[i,j-1] - go, A[i,j-1] - go, B[i,j-1] - ge)
            if(B[i,j] == M[i,j-1] - go) Traceback[["TrB"]][i,j] <- "LM" # L: Left
            else if(B[i,j] == A[i,j-1] - go) Traceback[["TrB"]][i,j] <- "LA"
            else Traceback[["TrB"]][i,j] <- "LB"
        }}

    return(new("NeedleObjAffine", Traceback = Traceback, M = M, A = A, B = B,
               GapOpen = go,  GapExten = ge, FreeEndGaps = FreeEndGaps))
}
