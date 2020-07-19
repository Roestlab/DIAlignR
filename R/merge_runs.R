#' Create a child run from two parent runs
#'
#' Get merged features and merged chromatograms from parent runs. Chromatograms are written on the disk
#' at dataPath/mzml. For each precursor aligned parent time-vectors and corresponding child time-vector
#' are also calculated and written as *_av.rda at dataPath.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-06-06
#' @inherit getChildXICs params
#' @param mergeName (string) name of the node that is generated with merging of runA and runB.
#' @param adaptiveRTs (environment) an empty environment used to store data for downstream analysis.
#' @param refRuns (environment) an empty environment used to store data for downstream analysis.
#' @param ropenms (pyopenms module) get this python module through get_ropenms().
#' @return (None)
#' @seealso \code{\link{childXICs}, \link{getChildXICs}, \link{traverseUp}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' params <- paramsDIAlignR()
#' fileInfo <- DIAlignR::getRunNames(dataPath = dataPath)
#' mzPntrs <- list2env(getMZMLpointers(fileInfo))
#' features <- list2env(getFeatures(fileInfo, maxFdrQuery = 1.00, runType = "DIA_proteomics"))
#' precursors <- getPrecursors(fileInfo, oswMerged = TRUE, runType = params[["runType"]],
#'  context = "experiment-wide", maxPeptideFdr = params[["maxPeptideFdr"]])
#' prec2chromIndex <- list2env(getChromatogramIndices(fileInfo, precursors, mzPntrs))
#' mergeName <- "master2"
#' adaptiveRTs <- new.env()
#' refRuns <- new.env()
#' \dontrun{
#' getNodeRun(runA = "run2", runB = "run0", mergeName = mergeName, dataPath = ".", fileInfo, features,
#'  mzPntrs, prec2chromIndex, precursors, params, adaptiveRTs)
#' file.remove(file.path(".", "mzml", paste0(mergeName, ".chrom.mzML")))
#' }
getNodeRun <- function(runA, runB, mergeName, dataPath, fileInfo, features, mzPntrs, prec2chromIndex,
                       precursors, params, adaptiveRTs, refRuns, ropenms){
  analytes <- precursors$transition_group_id
  refRun <- rep(NA_integer_, length(analytes))

  ##### Select reference for each feature. #####
  for(i in seq_along(analytes)){
    mscoreA <- suppressWarnings(min(features[[runA]]$m_score[features[[runA]]$transition_group_id == analytes[i]], na.rm = T))
    mscoreB <- suppressWarnings(min(features[[runB]]$m_score[features[[runB]]$transition_group_id == analytes[i]], na.rm = T))
    # If features are missing, min() will output Inf.
    if(mscoreA == Inf & mscoreB != Inf) {
      refRun[i] <- 2L
    } else if (mscoreB == Inf & mscoreA != Inf){
      refRun[i] <- 1L
    } else {
      refRun[i] <- ifelse(mscoreA > mscoreB, 2L, 1L)
    }
  }

  ##### Get childXICs #####
  mergedXICs_alignedVec <- getChildXICs(runA, runB, fileInfo, features, mzPntrs, precursors, prec2chromIndex,
                                        refRun, params)
  mergedXICs <- mergedXICs_alignedVec[[1]]
  alignedVecs <- mergedXICs_alignedVec[[2]]
  adaptiveRT <- mergedXICs_alignedVec[[3]]

  ##### Get merged features, calculate intensities, left width, right width, m-score. #####
  # we can also run pyopenms feature finder on new chromatogram.
  childFeature <- data.frame()
  for(i in seq_along(analytes)){
    if(is.null(mergedXICs[[i]])) next
    analyte <- analytes[i]
    alignedVec <- alignedVecs[[i]]
    XICs <- mergedXICs[[i]]
    df.A <-dplyr::filter(features[[runA]], .data$transition_group_id == analyte)
    df.B <-dplyr::filter(features[[runB]], .data$transition_group_id == analyte)
    if(refRun[i] == 1L) {
      df.ref <- df.A
      df.eXp <- df.B
    } else{
      df.ref <- df.B
      df.eXp <- df.A
    }
    if((nrow(df.ref) + nrow(df.eXp)) == 0) next
    alignedVec[, "alignedChildTime"] <- zoo::na.approx(alignedVec[, "alignedChildTime"], na.rm = FALSE)
    rows <- getChildFeature(XICs, alignedVec, df.ref, df.eXp, params)
    childFeature <- dplyr::bind_rows(childFeature, rows)
  }

  ##### Add node features #####
  assign("temp", childFeature, envir = parent.frame(n = 1))
  with(parent.frame(n = 1), features[[mergeName]] <- temp)

  ##### Write node mzML file #####
  mzmlName <- file.path(dataPath, "mzml", paste0(mergeName, ".chrom.mzML"))
  createMZML(ropenms, mzmlName, mergedXICs, precursors$transition_ids)

  ##### Add node run to fileInfo #####
  row <- data.frame(runName = mergeName, spectraFile = NA_character_, spectraFileID = NA_character_,
                    featureFile = NA_character_, chromatogramFile = mzmlName, row.names = mergeName)
  df <- dplyr::bind_rows(fileInfo, row)
  assign("temp", df, envir = parent.frame(n = 1))
  with(parent.frame(n = 1), fileInfo <- temp)

  ##### Add node run to mzPntrs #####
  mzPntr <- mzR::openMSfile(mzmlName, backend = "pwiz")
  assign("temp", mzPntr, envir = parent.frame(n = 1))
  with(parent.frame(n = 1), mzPntrs[[mergeName]] <- temp)

  ##### Add node run to prec2chromIndex #####
  prec2transition <- dplyr::select(precursors, .data$transition_group_id, .data$transition_ids) %>%
    tidyr::unnest(.data$transition_ids) %>% as.data.frame()
  chromHead <- mzR::chromatogramHeader(mzPntr) #TODO: Make sure that chromatogramIndex is read as integer64
  chromatogramIdAsInteger(chromHead) # Select only chromatogramId, chromatogramIndex
  df <- mapPrecursorToChromIndices(prec2transition, chromHead) # Get chromatogram Index for each precursor.
  assign("temp", df, envir = parent.frame(n = 1))
  with(parent.frame(n = 1), prec2chromIndex[[mergeName]] <- temp)

  ##### Write AlignedVecs #####
  filename <- file.path(dataPath, paste0(mergeName, "_av.rds"))
  saveRDS(alignedVecs, file = filename)

  ##### Add adaptiveRT to adaptiveRTs #####
  ab <- paste(runA, runB, sep = "_")
  ba <- paste(runB, runA, sep = "_")
  adaptiveRTs[[ab]] <- adaptiveRT[["ab"]]
  adaptiveRTs[[ba]] <- adaptiveRT[["ba"]]

  ##### Add refRun to refRuns #####
  refRuns[[mergeName]] <- refRun

  ##### Done #####
  message("Created a child run: ", mergeName)
}


#' Transform features to child time-domain
#'
#' This function transforms the peaks' times to child run's time-domain. The feature intensity is
#' calculated with appropriate method stated in params.
#' Internal missing values are not allowed in timeParent.
#' @inheritParams checkParams
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-07-17
#' @param XICs (list of data-frames) extracted ion chromatograms from the child run.
#' @param alignedVec (data-frame) aligned parent time-vectors. Must have five columns
#' @param df.ref (data-frame) contains reference-run features to be transformed. It has a format of \code{\link{getFeatures}} output.
#' @param df.eXp (data-frame) contains experiment-run features to be transformed. It has a format of \code{\link{getFeatures}} output.
#' @return (data-frame) this has a format of \code{\link{getFeatures}} output.
#' @keywords internal
#' @seealso \code{\link{trfrParentFeature}, \link{getNodeRun}}
#' @examples
#' data(masterXICs_DIAlignR, package="DIAlignR")
#' newXICs <- masterXICs_DIAlignR
#' params <- paramsDIAlignR()
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' fileInfo <- DIAlignR::getRunNames(dataPath = dataPath)
#' features <- getFeatures(fileInfo, maxFdrQuery = 1.00, runType = "DIA_proteomics")
#' df.ref <- features$run1[features$run1$transition_group_id == 4618L, ]
#' df.eXp <- features$run2[features$run2$transition_group_id == 4618L, ]
#' \dontrun{
#' getChildFeature(newXICs[[1]], newXICs[[2]], df.ref, df.eXp, params)
#' }
getChildFeature <- function(XICs, alignedVec, df.ref, df.eXp, params){
  # Convert Ref features to childXIC
  timeParent <- alignedVec[, c("tAligned.ref", "alignedChildTime")]
  colnames(timeParent) <- c("tAligned", "alignedChildTime")
  if(nrow(df.ref)!=0) df.ref <- trfrParentFeature(XICs, timeParent, df.ref, params)

  # Convert eXp features to childXIC
  timeParent <-alignedVec[, c("tAligned.eXp", "alignedChildTime")]
  colnames(timeParent) <- c("tAligned", "alignedChildTime")
  if(nrow(df.eXp)!=0) df.eXp <- trfrParentFeature(XICs, timeParent, df.eXp, params)

  # Convert Ref features to childXIC
  df <- dplyr::bind_rows(df.ref, df.eXp)
  df$child_pg_rank <- NA_integer_
  df$child_m_score <- df$m_score

  r <- 0L
  while (any(is.na(df$child_pg_rank))) {
    idx <- which.min(df$child_m_score)
    r <- r+1L
    df$child_pg_rank[idx] <- r
    df$child_m_score[idx] <- NA_real_
    pk <- c(df$leftWidth[idx], df$rightWidth[idx])
    # Remove other peaks that have conflicting boundaries
    rmv <- sapply(1:nrow(df), function(i) checkOverlap(pk, c(df$leftWidth[i], df$rightWidth[i])))
    rmv[idx] <- FALSE
    df <- df[!rmv,]
  }
  # First keep the peak based on m-score
  df$peak_group_rank <- df$child_pg_rank
  df <- df[order(df$peak_group_rank), -c(9,10)]
  rownames(df) <- NULL
  df
}


#' Transform features to child time-domain
#'
#' This function transforms the peaks' times to child run's time-domain. The feature intensity is
#' calculated with appropriate method stated in params.
#' Internal missing values are not allowed in timeParent.
#' @inheritParams checkParams
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-07-17
#' @param XICs (list of data-frames) extracted ion chromatograms from the child run.
#' @param timeParent (data-frame) has two columns: tAligned and alignedChildTime. tAligned is time vector from one of the parent run.
#' @param df (data-frame) contains features to be transformed. It has a format of \code{\link{getFeatures}} output.
#' @return (data-frame) this has a format of \code{\link{getFeatures}} output.
#' @keywords internal
#' @seealso \code{\link{getChildFeature}}
#' @examples
#' data(masterXICs_DIAlignR, package="DIAlignR")
#' newXICs <- masterXICs_DIAlignR
#' timeParent <- newXICs[[2]][, c("tAligned.ref", "alignedChildTime")]
#' colnames(timeParent) <- c("tAligned", "alignedChildTime")
#' params <- paramsDIAlignR()
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' fileInfo <- DIAlignR::getRunNames(dataPath = dataPath)
#' features <- getFeatures(fileInfo, maxFdrQuery = 1.00, runType = "DIA_proteomics")
#' df <- features$run1[features$run1$transition_group_id == 4618L, ]
#' \dontrun{
#' trfrParentFeature(newXICs[[1]], timeParent, df, params)
#' }
trfrParentFeature <- function(XICs, timeParent, df, params){
  # Transform the left boundary
  idx <- sapply(df$leftWidth, function(t) which.min(abs(timeParent[, "tAligned"] - t)))
  peak.left <- timeParent[idx, "alignedChildTime"]
  df$leftWidth <- peak.left

  # Transform the right boundary
  idx <- sapply(df$rightWidth, function(t) which.min(abs(timeParent[, "tAligned"] - t)))
  peak.right <- timeParent[idx, "alignedChildTime"]
  df$rightWidth <- peak.right

  # Transform the retention time
  idx <- sapply(df$RT, function(t) which.min(abs(timeParent[, "tAligned"] - t)))
  peak.RT <- timeParent[idx, "alignedChildTime"]
  df$RT <- peak.RT

  # Calculate peak area
  area <- sapply(1:nrow(df), function(i) calculateIntensity(XICs, df$leftWidth[i], df$rightWidth[i],
                  params[["integrationType"]], params[["baselineType"]],
                  params[["fitEMG"]], params[["baseSubtraction"]]))
  df$intensity <- area
  # If feature is outside of XICs we will get missing values
  df[!(is.na(df$leftWidth) | is.na(df$rightWidth) | is.na(df$intensity)),]
}

#' Overlap of two time ranges
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-07-17
#' @param x (numeric) must have only two values and x[2] > x[1].
#' @param y (numeric) must have only two values and y[2] > y[1].
#' @return (logical) TRUE: both time ranges overlap. FALSE: both time ranges do not overlap.
#' @keywords internal
#' @seealso \code{\link{getChildFeature}}
#' @examples
#' \dontrun{
#' checkOverlap(c(9.1, 13.1), c(2.1, 3.1))
#' checkOverlap(c(1.1, 3.1), c(3.2, 7.1))
#' }
checkOverlap <- function(x, y){
  # y has left boundary between x
  leftOverlap <- (y[1] - x[1]) >= 0 & (y[1] - x[2]) <= 0
  # y has right boundary between x
  rightOverlap <- (y[2] - x[1]) >= 0 & (y[2] - x[2]) <= 0
  # y is over-arching x
  overArch <- (y[2] - x[2]) >= 0 & (y[1] - x[1]) <= 0
  olap <- (leftOverlap | rightOverlap | overArch)
  olap
}

#' Develop child XICs for precursors
#'
#' This function performs the chromatogram alignment of all precursors across runA and runB. Aligned
#' chromatograms are merged into a child chromatogram. Aligned time vector and resulting child time
#' vector for each precursor is also returned.
#' @inherit traverseUp params
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-06-06
#' @param runA (string) name of a run to be merged with runB. Must be in the rownames of fileInfo.
#' @param runB (string) name of a run to be merged with runA. Must be in the rownames of fileInfo.
#' @param refRun (integer) must be of the same length as of precursors. 1: reference is runA, 2: reference is runB.
#' @return (list) has three elements. The first element has child XICs for all the precursors.
#' The second element has corresponding aligned time vectors. Third element contains Residual Standard
#' Errors (RSE) of global fits amongst runA and runB.
#' @seealso \code{\link{childXICs}, \link{getNodeRun}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' params <- paramsDIAlignR()
#' fileInfo <- DIAlignR::getRunNames(dataPath = dataPath)
#' mzPntrs <- getMZMLpointers(fileInfo)
#' features <- getFeatures(fileInfo, maxFdrQuery = 1.00, runType = "DIA_proteomics")
#' precursors <- getPrecursors(fileInfo, oswMerged = TRUE, runType = "DIA_proteomics",
#'  context = "experiment-wide", maxPeptideFdr = 0.05)
#' prec2chromIndex <- getChromatogramIndices(fileInfo, precursors, mzPntrs)
#' refRun <- rep(1L, nrow(precursors))
#' mergedXICs <- getChildXICs(runA = "run1", runB = "run2", fileInfo, features, mzPntrs,
#'   precursors, prec2chromIndex,refRun, params)
#' rm(mzPntrs)
#' @export
getChildXICs <- function(runA, runB, fileInfo, features, mzPntrs, precursors, prec2chromIndex, refRun,
                         params){
  analytes <- precursors$transition_group_id

  #### Get global alignment between runs ####
  pair <- paste(runA, runB, sep = "_")
  globalFit1 <- getGlobalAlignment(features, runA, runB,
                                   params[["globalAlignment"]], params[["globalAlignmentFdr"]], params[["globalAlignmentSpan"]])
  adaptiveRT1 <-  params[["RSEdistFactor"]]*getRSE(globalFit1)
  globalFit2 <- getGlobalAlignment(features, runB, runA,
                                   params[["globalAlignment"]], params[["globalAlignmentFdr"]], params[["globalAlignmentSpan"]])
  adaptiveRT2 <-  params[["RSEdistFactor"]]*getRSE(globalFit2)

  #### Get merged XICs ####
  mergedXICs <- vector(mode = "list", length = length(analytes))
  alignedVecs <- vector(mode = "list", length = length(analytes))
  for(i in seq_along(analytes)){
    analyte <- analytes[i]
    # TODO: ref and eXp should be analyte specific not global
    chromIndices.A <- prec2chromIndex[[runA]][["chromatogramIndex"]][[i]]
    chromIndices.B <- prec2chromIndex[[runB]][["chromatogramIndex"]][[i]]
    if(any(is.na(c(chromIndices.A, chromIndices.B)))){
      warning("Chromatogram indices for ", analyte, " are missing.")
      message("Skipping ", analyte, ".")
      mergedXICs[[i]] <- NULL
      next
    } else {
      XICs.A <- extractXIC_group(mz = mzPntrs[[runA]], chromIndices = chromIndices.A)
      XICs.A.s <- smoothXICs(XICs.A, type = params[["XICfilter"]], kernelLen = params[["kernelLen"]], polyOrd = params[["polyOrd"]])
      XICs.B <- extractXIC_group(mz = mzPntrs[[runB]], chromIndices = chromIndices.B)
      XICs.B.s <- smoothXICs(XICs.B, type = params[["XICfilter"]], kernelLen = params[["kernelLen"]], polyOrd = params[["polyOrd"]])
    }

    if(refRun[i] == 1L){
      XICs.ref <- XICs.A.s
      XICs.eXp <- XICs.B.s
      globalFit <- globalFit1
      adaptiveRT <- adaptiveRT1
    } else{
      XICs.ref <- XICs.B.s
      XICs.eXp <- XICs.A.s
      globalFit <- globalFit2
      adaptiveRT <- adaptiveRT2
    }

    #### Align chromatograms  ####
    alignedIndices <- getAlignedIndices(XICs.ref, XICs.eXp, globalFit, params[["alignType"]], adaptiveRT,
                                        params[["normalization"]], params[["simMeasure"]], params[["goFactor"]], params[["geFactor"]],
                                        params[["cosAngleThresh"]], params[["OverlapAlignment"]], params[["dotProdThresh"]], params[["gapQuantile"]],
                                        params[["kerLen"]], params[["hardConstrain"]], params[["samples4gradient"]], objType = "light")

    alignedIndices <- alignedIndices[, 1:2]
    #### Merge chromatogramsm ####
    merged_xics <- childXICs(XICs.ref, XICs.eXp, alignedIndices, params[["fillMethod"]], params[["polyOrd"]],
                             params[["kernelLen"]], params[["splineMethod"]], params[["mergeTime"]], params[["keepFlanks"]])
    mergedXICs[[i]] <- merged_xics[[1]]
    alignedVecs[[i]] <- merged_xics[[2]]
  }
  list(mergedXICs, alignedVecs, adaptiveRT = list(ab= adaptiveRT1, ba = adaptiveRT2))
}
