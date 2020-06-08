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
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' \dontrun{
#' data("multipeptide_DIAlignR", package = "DIAlignR")
#' getRefRun(multipeptide_DIAlignR)
#' }
#' @seealso \code{\link{getMultipeptide}}
#' @keywords internal
getRefRun <- function(multipeptide){
  refRun <- data.frame()
  for(pep in multipeptide){
    df <- pep[which.min(pep$m_score), c("transition_group_id", "run")]
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
#' writeTables("DIAlignR.csv", fileInfo, multipeptide, precursors)
#' }
writeTables <- function(filename, fileInfo, multipeptide, precursors){
  allIDs <- as.integer(names(multipeptide))
  runs <- rownames(fileInfo)

  ######### Initilize output tables. #######
  rtTbl <- lapply(1:length(allIDs), function(i) rep(NA, length(runs)))
  intesityTbl <- lapply(1:length(allIDs), function(i) rep(NA, length(runs)))
  lwTbl <- lapply(1:length(allIDs), function(i) rep(NA, length(runs)))
  rwTbl <- lapply(1:length(allIDs), function(i) rep(NA, length(runs)))
  prTbl <- lapply(1:length(allIDs), function(i) rep(NA, length(runs)))
  mTbl <- lapply(1:length(allIDs), function(i) rep(NA, length(runs)))
  arTbl <- lapply(1:length(allIDs), function(i) rep(NA, length(runs)))

  #### Fill tables. ###########
  for(i in seq_along(multipeptide)){
    for(run in seq_along(runs)){
      idx <- which(multipeptide[[i]][["run"]] == runs[run] & multipeptide[[i]][["alignment_rank"]] == 1L)
      if(length(idx) == 0) next
      rtTbl[[i]][run] <- multipeptide[[i]][["RT"]][idx]
      intesityTbl[[i]][run] <- multipeptide[[i]][["intensity"]][idx]
      lwTbl[[i]][run] <- multipeptide[[i]][["leftWidth"]][idx]
      rwTbl[[i]][run] <- multipeptide[[i]][["rightWidth"]][idx]
      prTbl[[i]][run] <- multipeptide[[i]][["peak_group_rank"]][idx]
      mTbl[[i]][run] <- multipeptide[[i]][["m_score"]][idx]
      arTbl[[i]][run] <- multipeptide[[i]][["alignment_rank"]][idx]
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
  colnames(rtTbl) <- fileInfo[["runName"]]
  colnames(intesityTbl) <- fileInfo[["runName"]]
  colnames(lwTbl) <- fileInfo[["runName"]]
  colnames(rwTbl) <- fileInfo[["runName"]]
  colnames(prTbl) <- fileInfo[["runName"]]
  colnames(mTbl) <- fileInfo[["runName"]]
  colnames(arTbl) <- fileInfo[["runName"]]

  ######### Add precursor column. ##############
  rtTbl[["precursor"]] <- allIDs
  intesityTbl[["precursor"]] <- allIDs
  lwTbl[["precursor"]] <- allIDs
  rwTbl[["precursor"]] <- allIDs
  prTbl[["precursor"]] <- allIDs
  mTbl[["precursor"]] <- allIDs
  arTbl[["precursor"]] <- allIDs

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

  utils::write.table(finalTbl, file = filename, sep = ",", row.names = FALSE)
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
