#' Fetch the reference run for each precursor
#'
#' Provides the reference run based on lowest m-score.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-04-08
#' @param multipeptide (list of data-frames) Each element of the list is collection of features associated with a precursor.
#' @return (dataframe) has two columns:
#' \item{transition_group_id}{(integer) a unique id for each precursor.}
#' \item{run}{(string) run identifier.}
#' @examples
#' data("multipeptide_DIAlignR", package = "DIAlignR")
#' \dontrun{
#' getRefRun(multipeptide_DIAlignR)
#' }
#' @seealso \code{\link{getMultipeptide}}
#' @keywords internal
getRefRun <- function(multipeptide){
  refRun <- data.frame()
  for(pep in multipeptide){
    idx <- which.min(pep$m_score)
    if(length(idx)==0) {
      df <- pep[1, c("transition_group_id", "run")]
      df$run <- NA_character_
    } else{
      df <- pep[idx, c("transition_group_id", "run")]
    }
    refRun <- rbind(refRun, df, make.row.names = FALSE)
  }
  refRun
}

#' Get multipeptides
#'
#' Each element of the multipeptide is a collection of features associated with a precursor.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-04-08
#' @importFrom bit64 NA_integer64_
#' @param precursors (data-frames) Contains precursors and associated transition IDs.
#' @param features (list of data-frames) Contains features and their properties identified in each run.
#' @return (list) of dataframes having following columns:
#' \item{transition_group_id}{(integer) a unique id for each precursor.}
#' \item{run}{(string) run identifier.}
#' \item{RT}{(numeric) retention time as in FEATURE.EXP_RT of osw files.}
#' \item{Intensity}{(numeric) peak intensity as in FEATURE_MS2.AREA_INTENSITY of osw files.}
#' \item{leftWidth}{(numeric) as in FEATURE.LEFT_WIDTH of osw files.}
#' \item{rightWidth}{(numeric) as in FEATURE.RIGHT_WIDTH of osw files.}
#' \item{peak_group_rank}{(integer) rank of each feature associated with transition_group_id.}
#' \item{m_score}{(numeric) q-value of each feature associated with transition_group_id.}
#' \item{alignment_rank}{(integer) rank of each feature post-alignment.}
#'
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' fileInfo <- getRunNames(dataPath, oswMerged = TRUE)
#' precursors <- getPrecursors(fileInfo, oswMerged = TRUE)
#' features <- getFeatures(fileInfo, maxFdrQuery = 0.05)
#' multipeptide <- getMultipeptide(precursors, features)
#' @seealso \code{\link{getPrecursors}, \link{getFeatures}}
#' @export
getMultipeptide <- function(precursors, features){
  multipeptide <- vector(mode = "list", length = length(precursors[["transition_group_id"]]))
  for(i in seq_along(multipeptide)){
    analyte <- precursors[["transition_group_id"]][i]
    multipeptide[[i]] <- data.frame()
    for(run in names(features)){
      # Match precursor in each run, if not found add NA
      index <- which(features[[run]][["transition_group_id"]] == analyte)
      if(length(index) != 0){
        df <- features[[run]][index, ]
        df["run"] <- run
      } else {
        df <- data.frame("transition_group_id" = analyte, "feature_id" = NA_integer64_, "run" = run,
                         "RT" = NA_real_, "intensity" = NA_real_,
                         "leftWidth" = NA_real_, "rightWidth" = NA_real_,
                         "peak_group_rank" = NA_integer_, "m_score" = NA_real_,
                         stringsAsFactors = FALSE)
      }
      multipeptide[[i]] <- rbind(multipeptide[[i]], df, make.row.names = FALSE)
    }
    multipeptide[[i]][["alignment_rank"]] <- NA_integer_
  }

  # Convert precursors as character. Add names to the multipeptide list.
  names(multipeptide) <- as.character(precursors[["transition_group_id"]])
  multipeptide
}

#' Writes the output table post-alignment
#'
#' Selects all features from multipeptide with alignment rank = 1, and write them onto the disk.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-04-14
#' @importFrom rlang .data
#' @importFrom tidyr pivot_longer
#' @param filename (string) Name of the output file.
#' @param fileInfo (data-frame) Output of getRunNames function.
#' @param multipeptide (list of data-frames) Each element of the list is collection of features associated with a precursor.
#' @param precursors (data-frame) for each transition_group_id, contains peptide sequence and charge.
#' @return An output table with following columns: precursor, run, intensity, RT, leftWidth, rightWidth,
#'  peak_group_rank, m_score, alignment_rank, peptide_id, sequence, charge, group_label.
#' @seealso \code{\link{getRunNames}, \link{getMultipeptide}, \link{getPrecursors}}
#' @keywords internal
#'
#' @examples
#' data(oswFiles_DIAlignR, package="DIAlignR")
#' \dontrun{
#' writeTables(fileInfo, multipeptide, precursors)
#' }
writeTables <- function(fileInfo, multipeptide, precursors){
  analytes <- precursors$transition_group_id
  runs <- rownames(fileInfo)
  idx <- grep("^master[0-9]+$", runs, invert = TRUE)
  runs <- runs[idx]
  runName <- fileInfo$runName[idx]

  ######### Initilize output tables. #######
  rtTbl <- lapply(1:length(analytes), function(i) rep(NA, length(runs)))
  intesityTbl <- lapply(1:length(analytes), function(i) rep(NA, length(runs)))
  lwTbl <- lapply(1:length(analytes), function(i) rep(NA, length(runs)))
  rwTbl <- lapply(1:length(analytes), function(i) rep(NA, length(runs)))
  prTbl <- lapply(1:length(analytes), function(i) rep(NA, length(runs)))
  mTbl <- lapply(1:length(analytes), function(i) rep(NA, length(runs)))
  arTbl <- lapply(1:length(analytes), function(i) rep(NA, length(runs))) # alignment rank table

  #### Fill tables. ###########
  for(i in seq_along(analytes)){
    analyte_chr <- as.character(analytes[i])
    df <- multipeptide[[analyte_chr]]
    for(run in seq_along(runs)){
      idx <- which(df[["run"]] == runs[run] & df[["alignment_rank"]] == 1L)
      if(length(idx) == 0) next
      rtTbl[[i]][run] <- df[["RT"]][idx]
      intesityTbl[[i]][run] <- df[["intensity"]][idx]
      lwTbl[[i]][run] <- df[["leftWidth"]][idx]
      rwTbl[[i]][run] <- df[["rightWidth"]][idx]
      prTbl[[i]][run] <- df[["peak_group_rank"]][idx]
      mTbl[[i]][run] <- df[["m_score"]][idx]
      arTbl[[i]][run] <- df[["alignment_rank"]][idx]
    }
  }

  ##### Convert lists into tables. ############
  rtTbl <- as.data.frame(do.call(rbind, rtTbl))
  intesityTbl <- as.data.frame(do.call(rbind, intesityTbl))
  lwTbl <- as.data.frame(do.call(rbind, lwTbl))
  rwTbl <- as.data.frame(do.call(rbind, rwTbl))
  prTbl <- as.data.frame(do.call(rbind, prTbl))
  mTbl <- as.data.frame(do.call(rbind, mTbl))
  arTbl <- as.data.frame(do.call(rbind, arTbl))

  ######## Assign column names. ###############
  colnames(rtTbl) <- runName
  colnames(intesityTbl) <- runName
  colnames(lwTbl) <- runName
  colnames(rwTbl) <- runName
  colnames(prTbl) <- runName
  colnames(mTbl) <- runName
  colnames(arTbl) <- runName

  ######### Add precursor column. ##############
  rtTbl[["precursor"]] <- analytes
  intesityTbl[["precursor"]] <- analytes
  lwTbl[["precursor"]] <- analytes
  rwTbl[["precursor"]] <- analytes
  prTbl[["precursor"]] <- analytes
  mTbl[["precursor"]] <- analytes
  arTbl[["precursor"]] <- analytes

  ##### Changing format from wide to long. ############
  rtTbl <- pivot_longer(rtTbl, -.data$precursor, names_to = "run", values_to = "RT")
  intesityTbl <- pivot_longer(intesityTbl, -.data$precursor, names_to = "run", values_to = "intensity")
  lwTbl <- pivot_longer(lwTbl, -.data$precursor, names_to = "run", values_to = "leftWidth")
  rwTbl <- pivot_longer(rwTbl, -.data$precursor, names_to = "run", values_to = "rightWidth")
  prTbl <- pivot_longer(prTbl, -.data$precursor, names_to = "run", values_to = "peak_group_rank")
  mTbl <- pivot_longer(mTbl, -.data$precursor, names_to = "run", values_to = "m_score")
  arTbl <- pivot_longer(arTbl, -.data$precursor, names_to = "run", values_to = "alignment_rank")

  ##### Merging all into one table. ###################
  finalTbl <- merge(x = intesityTbl, y = rtTbl, by = c("precursor", "run"), all = TRUE)
  finalTbl <- merge(x = finalTbl, y = lwTbl, by = c("precursor", "run"), all = TRUE)
  finalTbl <- merge(x = finalTbl, y = rwTbl, by = c("precursor", "run"), all = TRUE)
  finalTbl <- merge(x = finalTbl, y = prTbl, by = c("precursor", "run"), all = TRUE)
  finalTbl <- merge(x = finalTbl, y = mTbl, by = c("precursor", "run"), all = TRUE)
  finalTbl <- merge(x = finalTbl, y = arTbl, by = c("precursor", "run"), all = TRUE)

  ##### Merging precursor information. ###################
  finalTbl <- merge(x = finalTbl, y = precursors[,-which(names(precursors) %in% c("transition_ids"))],
                    by.x = "precursor", by.y = "transition_group_id", all.x = TRUE)
  finalTbl
}

# Alignment of precursor 4618 , sequence = QFNNTDIVLLEDFQK/3 across runs
# ref = hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt
# eXp = hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"
# Example: test_getAlignObj
testAlignObj <- function(){
  AlignObj <- new("AffineAlignObjLight",
                    indexA_aligned = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,0,20,0,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,0,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,0,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,0,162,0,163,164,165,166,0,167,0,168,169,170,171,172,173,174,175,176),
                    indexB_aligned = c(0,0,0,1,0,2,0,3,0,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,0,92,0,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176),
                    score = c(0,0,0,2.675751,2.385165,4.745081,4.454496,6.67706,6.386474,9.641135,15.48268,26.87968,46.67034,77.78939,121.8607,170.8698,214.2358,244.4554,262.8537,262.5631,270.8538,270.5632,288.2477,319.3287,364.6418,413.4873,453.1931,479.6588,496.5207,506.6607,512.7442,515.267,516.8242,518.2747,519.8424,521.2872,522.6472,524.2912,525.6285,526.5892,526.9713,527.1531,527.4022,527.1116,530.9457,547.7525,588.2834,658.8819,748.3079,833.8337,898.3289,935.5809,948.8015,952.0709,952.8035,953.4267,954.0863,954.8143,955.4842,956.0834,956.802,957.535,958.2853,959.0355,959.7972,960.7983,961.8922,963.0142,964.2597,965.5837,966.878,968.0037,968.4412,968.1507,968.1958,968.9242,985.3144,1085.833,1364.976,1846.928,2409.31,2869.416,3132.509,3231.061,3257.015,3264.422,3269.377,3275.003,3282.515,3290.524,3297.864,3304.43,3310.324,3314.403,3316.806,3317.992,3318.933,3318.642,3319.328,3319.038,3320.17,3321.781,3323.71,3325.64,3327.855,3330.382,3332.989,3335.319,3337.555,3339.96,3342.381,3344.48,3346.456,3348.605,3350.446,3352.092,3353.829,3355.911,3358.256,3360.576,3363.292,3367.099,3372.687,3380.124,3389.957,3401.498,3414.81,3428.762,3441.046,3451.052,3459.235,3466.392,3473.212,3480.14,3490.173,3506.584,3530.062,3561.003,3595.718,3624.828,3642.574,3650.352,3653.893,3656.295,3658.798,3661.361,3663.704,3665.936,3667.714,3669.478,3670.721,3671.991,3673.278,3674.689,3676.068,3677.317,3678.688,3680.062,3681.513,3683.097,3684.786,3686.565,3688.24,3689.741,3690.859,3690.568,3691.496,3691.205,3692.31,3693.465,3694.458,3695.352,3695.061,3695.892,3695.602,3696.512,3697.468,3698.18,3698.799,3699.363,3699.94,3700.634,3701.585,3702.988))

  AlignObj
}

#' Checks all the alignment parameters
#'
#' Parameters are defined in the \code{\link{paramsDIAlignR}} function. If parameters are found
#' incongruous, the program will terminate.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-07-01
#' @param params (list) parameters are entered as list. Output of the \code{\link{paramsDIAlignR}} function.
#' @examples
#' params <- paramsDIAlignR()
#' \dontrun{
#' checkParams(params)
#' }
#' @seealso \code{\link{paramsDIAlignR}}
#' @keywords internal
checkParams <- function(params){
  if(params[["context"]] != "experiment-wide" & params[["context"]] != "global"){
    stop("context must either be experiment-wide or global for the alignment.")
  }

  if(params[["maxFdrQuery"]]>1 | params[["maxFdrQuery"]] < 0){
    stop("maxFdrQuery must be between 0 and 1. Recommended value is 0.05")
  }

  if(params[["maxPeptideFdr"]] < 0 | params[["maxPeptideFdr"]] > 1){
    stop("maxPeptideFdr must be between 0 and 1. Recommended value is 0.05")
  }

  if(!any(params[["XICfilter"]] %in% c("sgolay",  "boxcar", "gaussian", "loess", "none"))){
    stop("XICfilter must either be sgolay, boxcar, gaussian, loess or none.")
  }

  if(params[["XICfilter"]] == "sgolay"){
    if( (params[["kernelLen"]] %% 2) != 1){
      stop("For XICfilter = sgolay, kernelLen can only be odd number.")
    }
    if(params[["polyOrd"]] > params[["kernelLen"]]) {
      stop("For XICfilter = sgolay, polyOrd must be less than kernelLen.")
    }
  }

  if(!any(params[["globalAlignment"]] %in% c("loess",  "linear"))){
    stop("globalAlignment must either be loess or linear.")
  }

  if(params[["globalAlignment"]] == "loess" & (params[["globalAlignmentSpan"]] >1 | params[["globalAlignmentSpan"]] < 0)){
    stop("For globalAlignment = loess, globalAlignmentSpan must be between 0 and 1.")
  }

  if(params[["globalAlignmentFdr"]] < 0 | params[["globalAlignmentFdr"]] > params[["maxFdrQuery"]]){
    stop("globalAlignmentFdr must be between 0 and maxFdrQuery.")
  }

  if(!any(params[["normalization"]] %in% c("mean", "l2", "none"))){
    stop("normalization must be either mean, l2 or none.")
  }

  if(!any(params[["simMeasure"]] %in% c("dotProduct", "cosineAngle", "crossCorrelation", "cosine2Angle",
                            "dotProductMasked", "euclideanDist", "covariance", "correlation"))){
    stop("Choose appropriate simMeasure.")
  }

  if(!any(params[["alignType"]] %in% c("global", "local", "hybrid"))){
    stop("alignType must be either global, local or hybrid. Recommended value is hybrid.")
  }

  if(params[["unalignedFDR"]] < 0 | params[["unalignedFDR"]] > params[["maxFdrQuery"]]){
    stop("unalignedFDR must be between 0 and maxFdrQuery. Recommended value is 0.01")
  }

  if(params[["alignedFDR"]] < params[["unalignedFDR"]] | params[["alignedFDR"]] > params[["maxFdrQuery"]]){
    stop("alignedFDR must be between unalignedFDR and maxFdrQuery. Recommended value is 0.05")
  }

  if(!any(params[["integrationType"]] %in% c("intensity_sum", "trapezoid", "simpson"))){
    stop("integrationType must be either intensity_sum, trapezoid or simpson.")
  }

  if(!any(params[["baselineType"]] %in% c("none","base_to_base", "vertical_division_min", "vertical_division_max"))){
    stop("baselineType must be either base_to_base, vertical_division_min or vertical_division_max.")
  }

  # DOI: 10.1002/cem.1343
  if(params[["fitEMG"]] == TRUE){
    warning("Exponentially Modified Gaussian is not implemented yet.")
  }

  if(params[["analyteFDR"]] < 0 | params[["analyteFDR"]] > 1){
    # Find some logic for analyteFDR.
  }
  invisible(NULL)
}

#' Parameters for the alignment functions
#'
#' Retention alignment requires OpenSWATH/pyProphet extracted features and chromatograms. This function provides
#' a suite of parameters used for selecting features and manipulating chromatograms. Chromatogram
#' alignment can be performed via reference based or progressively via rooted or unrooted tree. This
#' function provides sensible parameters for these tasks.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-07-11
#' @return A list of parameters:
#' \item{runType}{(string) must be one of the strings "DIA_proteomics", "DIA_Metabolomics".}
#' \item{maxFdrQuery}{(numeric) a numeric value between 0 and 1. It is used to filter peptides from osw file which have SCORE_MS2.QVALUE less than itself.}
#' \item{maxPeptideFdr}{(numeric) a numeric value between 0 and 1. It is used to filter peptides from osw file which have SCORE_PEPTIDE.QVALUE less than itself.}
#' \item{analyteFDR}{(numeric) defines the upper limit of FDR on a precursor to be considered for multipeptide.}
#' \item{context}{(string) used in pyprophet peptide. Must be either "run-specific", "experiment-wide", or "global".}
#' \item{unalignedFDR}{(numeric) must be between 0 and maxFdrQuery. Features below unalignedFDR are
#'  considered for quantification even without the RT alignment.}
#' \item{alignedFDR}{(numeric) must be between unalignedFDR and 1. Features below alignedFDR are
#'  considered for quantification after the alignment.}
#' \item{integrationType}{(string) method to ompute the area of a peak contained in XICs. Must be
#'  from "intensity_sum", "trapezoid", "simpson".}
#' \item{baseSubtraction}{{logical} TRUE: remove background from peak signal using estimated noise levels.}
#' \item{baselineType}{(string) method to estimate the background of a peak contained in XICs. Must be
#'  from "none", "base_to_base", "vertical_division_min", "vertical_division_max".}
#' \item{fitEMG}{(logical) enable/disable exponentially modified gaussian peak model fitting.}
#' \item{recalIntensity}{(logical) recalculate intensity for all analytes.}
#' \item{fillMissing}{(logical) calculate intensity for ananlytes for which features are not found.}
#' \item{XICfilter}{(string) must be either sgolay, boxcar, gaussian, loess or none.}
#' \item{polyOrd}{(integer) order of the polynomial to be fit in the kernel.}
#' \item{kernelLen}{(integer) number of data-points to consider in the kernel.}
#' \item{globalAlignment}{(string) must be either "loess" or "linear".}
#' \item{globalAlignmentFdr}{(numeric) a numeric value between 0 and 1. Features should have m-score lower than this value for participation in LOESS fit.}
#' \item{globalAlignmentSpan}{(numeric) spanvalue for LOESS fit. For targeted proteomics 0.1 could be used.}
#' \item{RSEdistFactor}{(numeric) defines how much distance in the unit of rse remains a noBeef zone.}
#' \item{normalization}{(string) must be selected from "mean", "l2".}
#' \item{simMeasure}{(string) must be selected from dotProduct, cosineAngle, crossCorrelation,
#'   cosine2Angle, dotProductMasked, euclideanDist, covariance and correlation.}
#' \item{alignType}{(numeric) available alignment methods are "global", "local" and "hybrid".}
#' \item{goFactor}{(numeric) penalty for introducing first gap in alignment. This value is multiplied by base gap-penalty.}
#' \item{geFactor}{(numeric) penalty for introducing subsequent gaps in alignment. This value is multiplied by base gap-penalty.}
#' \item{cosAngleThresh}{(numeric) in simType = dotProductMasked mode, angular similarity should be higher than cosAngleThresh otherwise similarity is forced to zero.}
#' \item{OverlapAlignment}{(logical) an input for alignment with free end-gaps. False: Global alignment, True: overlap alignment.}
#' \item{dotProdThresh}{(numeric) in simType = dotProductMasked mode, values in similarity matrix higher than dotProdThresh quantile are checked for angular similarity.}
#' \item{gapQuantile}{(numeric) must be between 0 and 1. This is used to calculate base gap-penalty from similarity distribution.}
#' \item{kerLen}{(integer) In simType = crossCorrelation, length of the kernel used to sum similarity score. Must be an odd number.}
#' \item{hardConstrain}{(logical) if FALSE; indices farther from noBeef distance are filled with distance from linear fit line.}
#' \item{samples4gradient}{(numeric) modulates penalization of masked indices.}
#' \item{fillMethod}{(string) must be either "spline", "sgolay" or "linear".}
#' \item{splineMethod}{(string) must be either "fmm" or "natural".}
#' \item{mergeTime}{(string) must be either "ref", "avg", "refStart" or "refEnd".}
#' \item{keepFlanks}{(logical) TRUE: Flanking chromatogram is not removed.}
#' @seealso \code{\link{checkParams}, \link{alignTargetedRuns}}
#' @examples
#' params <- paramsDIAlignR()
#' @export
paramsDIAlignR <- function(){
  params <- list( runType = "DIA_proteomics", maxFdrQuery = 0.05, maxPeptideFdr = 0.05, analyteFDR = 0.05,
                  context = "experiment-wide", unalignedFDR = 0.01, alignedFDR = 0.05,
                  integrationType = "intensity_sum", baselineType = "base_to_base", fitEMG = FALSE,
                  recalIntensity = FALSE, fillMissing = TRUE, baseSubtraction = TRUE,
                  XICfilter = "sgolay", polyOrd = 4, kernelLen = 9,
                  globalAlignment = "loess", globalAlignmentFdr = 0.01, globalAlignmentSpan = 0.1,
                  RSEdistFactor = 3.5, normalization = "mean", simMeasure = "dotProductMasked",
                  alignType = "hybrid", goFactor = 0.125, geFactor = 40,
                  cosAngleThresh = 0.3, OverlapAlignment = TRUE,
                  dotProdThresh = 0.96, gapQuantile = 0.5, kerLen = 9,
                  hardConstrain = FALSE, samples4gradient = 100,
                  fillMethod = "spline", splineMethod = "fmm", mergeTime = "avg", smoothPeakArea = FALSE,
                  keepFlanks = FALSE)
  params
}

#' Prints alignment summary
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-07-15
#' @param finalTbl (dataframe) an output of  \code{\link{writeTables}} function.
#' @param params (list) must have following elements: alignedFDR and unalignedFDR.
#' @seealso \code{\link{paramsDIAlignR}, \link{writeTables}}
#' @keywords internal
alignmentStats <- function(finalTbl, params){
  # Without alignment at unaligned FDR:
  num1 <- sum(finalTbl$m_score <= params[["unalignedFDR"]] & finalTbl$alignment_rank == 1L &
                !is.na(finalTbl$intensity), na.rm = TRUE)
  message("The number of quantified precursors at ", params[["unalignedFDR"]], " FDR: ", num1)

  # Without alignment at aligned FDR (Gain):
  num2 <- sum(finalTbl$m_score > params[["unalignedFDR"]] & finalTbl$alignment_rank == 1L &
                !is.na(finalTbl$intensity), na.rm = TRUE)
  message("The increment in the number of quantified precursors at ", params[["alignedFDR"]],
          " FDR: ", num2)

  # Corrected peptides by Alignmet
  idx <- finalTbl$peak_group_rank != 1L & finalTbl$m_score > params[["unalignedFDR"]] &
    finalTbl$alignment_rank == 1L & !is.na(finalTbl$intensity)
  num3 <- sum(idx, na.rm = TRUE)
  message("Out of ", num2, " DIAlignR corrects the peaks for ", num3, " precursors.")
  message("Hence, it has corrected quantification of ", round(num3*100/(num2 + num1), 3), "% precursors.")

  # Gain by calculating area of missing features:
  num4 <- sum(!is.na(finalTbl$intensity) & is.na(finalTbl$m_score) & finalTbl$alignment_rank == 1L, na.rm = TRUE)
  message("DIAlignR has calculated quantification for ", num4, " precursors, for which peaks were not identified.")
  message("Thus, it provides a gain of ", round(num4*100/(num2 + num1 + num4), 3), "%.")

  if(params[["fillMissing"]]) message(sum(is.na(finalTbl$intensity)),
                          " precursors had part of the aligned peak out of the chromatograms or missing chromatograms, hence could not be quantified.")
}

#' Prints messages if a certain number of analytes are aligned
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-07-26
#' @keywords internal
updateOnalignTargetedRuns <-  function(i){
  if(i < 5){
    message(i, " precursors have been aligned.")
  } else if(i < 1000){
    if(i %% 100 == 0) message(i, " precursors have been aligned.")
  } else {
    if(i %% 1000 == 0) message(i, " precursors have been aligned.")
  }
}

envName <- "TricEnvr"
# helper function to skip tests if we don't have the 'foo' module
skip_if_no_pyopenms <- function() {
  ropenms <- try(get_ropenms(condaEnv = envName, useConda=TRUE))
  no_ropenms <- is(ropenms, "try-error")
  if(no_ropenms)
    skip("ropenms not available for testing. A conda environment with name TricEnvr is MUST for testing.")
}
