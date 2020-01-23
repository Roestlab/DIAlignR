#' An S4 object for class AffineAlignObj
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-14
#' @importFrom methods setClass new
#' @seealso \code{\link{doAffineAlignmentCpp}}
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
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-14
#' @importFrom methods setClass new
#' @seealso \code{\link{doAffineAlignmentCpp}}
#' @export
AffineAlignObjLight <- setClass(Class="AffineAlignObjLight",
                           representation(indexA_aligned = "numeric", indexB_aligned = "numeric",
                                          score = "numeric")
)

#' An S4 object for class AffineAlignObj. It only contains similarity matrix and aligned indices.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-14
#' @importFrom methods setClass new
#' @seealso \code{\link{doAffineAlignmentCpp}}
#' @export
AffineAlignObjMedium <- setClass(Class="AffineAlignObjMedium",
                        representation(s = "matrix", path = "matrix",
                                       indexA_aligned = "numeric", indexB_aligned = "numeric",
                                               score = "numeric")
)


#' An S4 object for class AlignObj
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-14
#' @importFrom methods setClass new
#' @seealso \code{\link{doAlignmentCpp}}
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
