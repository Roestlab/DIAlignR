#' Extracted-ion chromatograms (XICs) of a peptide
#'
#' @description
#'  XICs of peptide QFNNTDIVLLEDFQK/3 from three SWATH runs:
#'
#' run0 : hroest_K120808_Strep10\%PlasmaBiolRepl1_R03_SW_filt.chrom.mzML\cr
#' run1 : hroest_K120809_Strep0\%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML\cr
#' run2 : hroest_K120809_Strep10\%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML
#'
#' @format A list of three elements where each element consists of a list of six data frames. Each
#'  data frame has two columns:
#' \describe{
#'   \item{time}{Retention time of ananlyte in the run, in sec}
#'   \item{intensity}{Intensity of signal for the transition}
#' }
#' @source Raw files are downloaded from \href{http://www.peptideatlas.org/PASS/PASS01508}{Peptide Atlas}.
#' File test_GenerateData.R has \href{https://github.com/shubham1637/DIAlignR/tree/master/data-raw}{source code}
#' to generate the example data.
"XIC_QFNNTDIVLLEDFQK_3_DIAlignR"

#' Analytes information from osw files
#'
#' @description
#'  analytes info from three SWATH runs:
#'
#' run0 : hroest_K120808_Strep10\%PlasmaBiolRepl1_R03_SW_filt.chrom.mzML\cr
#' run1 : hroest_K120809_Strep0\%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML\cr
#' run2 : hroest_K120809_Strep10\%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML
#'
#' @format A list of three elements where each element consists of a dataframe:
#' \describe{
#'   \item{transition_group_id}{ID of each peptide}
#'   \item{RT}{Retention time, in sec}
#'   \item{intensity}{Inensity of associated feature}
#'   \item{leftWidth}{Left width of the peak, in sec}
#'   \item{rightWidth}{Right width of the peak, in sec}
#'   \item{peak_group_rank}{Ranking of associated feature}
#'   \item{m_score}{qvalue of associated feature}
#' }
#' @source Raw files are downloaded from \href{http://www.peptideatlas.org/PASS/PASS01508}{Peptide Atlas}.
#' File test_GenerateData.R has \href{https://github.com/shubham1637/DIAlignR/tree/master/data-raw}{source code}
#' to generate the example data.
"oswFiles_DIAlignR"


#' Analytes information from multipeptide.
#'
#' @description
#'  analytes info from three SWATH runs:
#'
#' run0 : hroest_K120808_Strep10\%PlasmaBiolRepl1_R03_SW_filt.chrom.mzML\cr
#' run1 : hroest_K120809_Strep0\%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML\cr
#' run2 : hroest_K120809_Strep10\%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML
#'
#' @format A list of 199 elements where each element represents a precursor and consists of a dataframe:
#' \describe{
#'   \item{transition_group_id}{ID of each precursor. Same as the name of the list}
#'   \item{RT}{Retention time, in sec}
#'   \item{intensity}{Inensity of associated feature}
#'   \item{leftWidth}{Left width of the peak, in sec}
#'   \item{rightWidth}{Right width of the peak, in sec}
#'   \item{peak_group_rank}{Ranking of associated feature}
#'   \item{m_score}{qvalue of associated feature}
#'   \item{run}{Name of the run, feature is from}
#'   \item{alignment_rank}{Rank of the feature after alignment}
#' }
#' @source Raw files are downloaded from \href{http://www.peptideatlas.org/PASS/PASS01508}{Peptide Atlas}.
#' File test_GenerateData.R has \href{https://github.com/shubham1637/DIAlignR/tree/master/data-raw}{source code}
#' to generate the example data.
"multipeptide_DIAlignR"


#' Alignment object of a peptide.
#'
#' @description
#'  Aligned XICs of peptide 14299_QFNNTDIVLLEDFQK/3 across two SWATH runs:
#'
#' run1 : hroest_K120809_Strep0\%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML\cr
#' run2 : hroest_K120809_Strep10\%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML
#'
#' @format A S4 object of 16 slots:
#' \describe{
#'   \item{s}{similarity score matrix.}
#'   \item{M}{Match or Mismatch matrix, residues of A and B are aligned without a gap. M(i,j) = Best score upto (i,j) given Ai is aligned to Bj.}
#'   \item{A}{Insert in sequence A, residue in A is aligned to gap in B. A(i,j) is the best score given that Ai is aligned to a gap in B.}
#'   \item{B}{Insert in sequence B, residue in B is aligned to gap in A. B(i,j) is the best score given that Bj is aligned to a gap in A.}
#'   \item{Traceback}{Traceback matrices store source matrix name and direction as matrices are filled with dynamic programming.}
#'   \item{path}{Path matrix would represent alignment path through similarity matrix as binary-hot encoding.}
#'   \item{signalA_len}{Number of data-points in signal A.}
#'   \item{signalB_len}{Number of data-points in signal B.}
#'   \item{GapOpen}{Penalty for Gap opening. For n consecutive gaps: Penalty = GapOpen + (n-1)*GapExten.}
#'   \item{GapExten}{Penalty for Gap extension. For n consecutive gaps: Penalty = GapOpen + (n-1)*GapExten.}
#'   \item{FreeEndGaps}{True for Overlap alignment.}
#'   \item{indexA_aligned}{Aligned signalA indices after affine alignment.}
#'   \item{indexB_aligned}{Aligned signalB indices after affine alignment.}
#'   \item{score}{Cumulative score along the aligned path.}
#'   \item{simScore_forw}{Not needed, will be removed.}
#'   \item{nGaps}{Total number of gaps in the alignment path.}
#' }
#' @source C++ code is exaplained at \href{https://shubham1637.github.io/DIAlignR/src/doc/html/structDIAlign_1_1AffineAlignObj.html}{DIAlign namespace}.
#' File test_GenerateData.R has \href{https://github.com/shubham1637/DIAlignR/tree/master/data-raw}{source code}
#' to generate the example data.
"alignObj_DIAlignR"
