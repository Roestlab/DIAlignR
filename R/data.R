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
#' @source data-raw/test_GenerateData.R has source code.
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
#'   \item{filename}{Name of corresponding mzML file}
#'   \item{RT}{Retention time, in sec}
#'   \item{delta_rt}{Retention time difference to library, in sec}
#'   \item{assay_RT}{Library retention time, in sec}
#'   \item{Intensity}{Inensity of associated feature}
#'   \item{leftWidth}{Left width of the peak, in sec}
#'   \item{rightWidth}{Right width of the peak, in sec}
#'   \item{peak_group_rank}{Ranking of associated feature}
#'   \item{m_score}{qvalue of associated feature}
#'   \item{chromatogramIndex}{Indices of XICs of the feature in corresponding mzML}
#'   \item{transition_ids}{Transition IDs (chromatogram IDs) of the feature}
#' }
#' @source data-raw/test_GenerateData.R has source code.
"oswFiles_DIAlignR"
