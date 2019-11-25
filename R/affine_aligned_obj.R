#' An S4 object for class AffineAlignObj
#'
#' @importFrom methods setClass
#' @export
AffineAlignObj <- setClass(Class="AffineAlignObj",
                   representation(s = "matrix", M = "matrix", A = "matrix", B = "matrix",
                                  Traceback = "matrix",
                                  path = "matrix",
                                  signalA_len = "numeric", signalB_len = "numeric",
                                  GapOpen = "numeric", GapExten = "numeric",
                                  FreeEndGaps = "logical",
                                  indexA_aligned = "numeric", indexB_aligned = "numeric",
                                  score = "numeric", simScore_forw = "numeric",
                                  nGaps = "numeric")
                   )

#' An S4 object for class AffineAlignObj. It only contains aligned indices.
#'
#' @importFrom methods setClass
#' @export
AffineAlignObjLight <- setClass(Class="AffineAlignObj",
                           representation(indexA_aligned = "numeric", indexB_aligned = "numeric",
                                          score = "numeric")
)

#' An S4 object for class AffineAlignObj. It only contains similarity matrix and aligned indices.
#'
#' @importFrom methods setClass
#' @export
AffineAlignObjMedium <- setClass(Class="AffineAlignObj",
                        representation(s = "matrix", path = "matrix",
                                       indexA_aligned = "numeric", indexB_aligned = "numeric",
                                               score = "numeric")
)


#' An S4 object for class AlignObj
#'
#' @importFrom methods setClass
#' @export
AlignObj <- setClass(Class="AlignObj",
                     representation(s = "matrix", M = "matrix", Traceback = "matrix",
                                    path = "matrix", optionalPaths = "matrix",
                                    M_forw = "matrix",
                                    signalA_len = "numeric", signalB_len = "numeric",
                                    GapOpen = "numeric", GapExten = "numeric",
                                    FreeEndGaps = "logical",
                                    indexA_aligned = "numeric", indexB_aligned = "numeric",
                                    score = "numeric", score_forw = "numeric", simScore_forw = "numeric",
                                    nGaps = "numeric")
)
