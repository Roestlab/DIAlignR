#' An S4 object for class AffineAlignObj
#'
#' s is a point-wise similarity matrix between signalA and signalB.
#' Intermediate matrices M,A,B are calculated from s for affine-alignment. Each cell of the Traceback
#' matrix has coordinate of its parent cell. path matrix is a binary matrix with ones indicating path
#' of maximum cumulative score.
#' GapOpen and GapExten are gap-opening and gap-extension penalties used by affine alignment algorithm.
#' indexA_aligned and indexB_aligned are aligned indices of signalA and SignalB. The cumulative
#' alignment score is in score vector.
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
#' indexA_aligned and indexB_aligned are aligned indices of signalA and SignalB. The cumulative
#' alignment score is in score vector.
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
#' s is a point-wise similarity matrix between signalA and signalB.
#' path matrix is a binary matrix with ones indicating path of maximum cumulative score.
#' GapOpen and GapExten are gap-opening and gap-extension penalties used by affine alignment algorithm.
#' indexA_aligned and indexB_aligned are aligned indices of signalA and SignalB. The cumulative
#' alignment score is in score vector.
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
#' s is a point-wise similarity matrix between signalA and signalB.
#' Intermediate matrices M is calculated from s for alignment. Each cell of the Traceback
#' matrix has coordinate of its parent cell. path matrix is a binary matrix with ones indicating path
#' of maximum cumulative score.
#' GapOpen and GapExten are gap-opening and gap-extension penalties used by alignment algorithm. They
#' must be the same. indexA_aligned and indexB_aligned are aligned indices of signalA and SignalB.
#' The cumulative alignment score is in score vector.
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
