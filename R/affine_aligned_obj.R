#' An S4 object for class AffineAlignObj
#'
#' @importFrom methods setClass
#' @export
AffineAlignObj <- setClass(Class="AffineAlignObj",
                   representation(M = "matrix", A = "matrix", B = "matrix", Traceback = "matrix",
                                  signalA_len = "numeric", signalB_len = "numeric",
                                  GapOpen = "numeric", GapExten = "numeric",
                                  FreeEndGaps = "logical",
                                  indexA_aligned = "numeric", indexB_aligned = "numeric",
                                  score = "numeric")
                   )

#' An S4 object for class AlignObj
#'
#' @importFrom methods setClass
#' @export
AlignObj <- setClass(Class="AlignObj",
                            representation(M = "matrix", Traceback = "matrix",
                                           signalA_len = "numeric", signalB_len = "numeric",
                                           GapOpen = "numeric", GapExten = "numeric",
                                           FreeEndGaps = "logical",
                                           indexA_aligned = "numeric", indexB_aligned = "numeric",
                                           score = "numeric")
)
