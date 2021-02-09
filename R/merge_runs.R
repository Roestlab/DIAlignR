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
#' @importFrom data.table data.table set
#' @import RSQLite DBI
#' @inherit getChildXICs params
#' @inheritParams traverseUp
#' @param mergeName (string) name of the node that is generated with merging of runA and runB.
#' @param adaptiveRTs (environment) an empty environment used to store data for downstream analysis.
#' @param refRuns (environment) an empty environment used to store data for downstream analysis.
#' @return (None)
#' @seealso \code{\link{childXICs}, \link{getChildXICs}, \link{traverseUp}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' params <- paramsDIAlignR()
#' fileInfo <- getRunNames(dataPath = dataPath)
#' mzPntrs <- list2env(getMZMLpointers(fileInfo))
#' features <- list2env(getFeatures(fileInfo, maxFdrQuery = 1.00, runType = "DIA_proteomics"))
#' precursors <- getPrecursors(fileInfo, oswMerged = TRUE, runType = params[["runType"]],
#'  context = "experiment-wide", maxPeptideFdr = params[["maxPeptideFdr"]])
#' precursors <- dplyr::arrange(precursors, .data$peptide_id, .data$transition_group_id)
#' peptideIDs <- unique(precursors$peptide_id)
#' peptideScores <- getPeptideScores(fileInfo, peptideIDs, oswMerged = TRUE, params[["runType"]], params[["context"]])
#'  masters <- paste("master", 1:(nrow(fileInfo)-1), sep = "")
#' peptideScores <- lapply(peptideIDs, function(pep) {x <- peptideScores[.(pep)][,-c(1L)]
#'   x <- rbindlist(list(x, data.table("run" = masters, "score" = NA_real_, "pvalue" = NA_real_,
#'     "qvalue" = NA_real_)), use.names=TRUE)
#'   setkeyv(x, "run"); x})
#' names(peptideScores) <- as.character(peptideIDs)
#' multipeptide <- getMultipeptide(precursors, features)
#' prec2chromIndex <- list2env(getChromatogramIndices(fileInfo, precursors, mzPntrs))
#' mergeName <- "master1"
#' adaptiveRTs <- new.env()
#' refRuns <- new.env()
#' \dontrun{
#' multipeptide <- getNodeRun(runA = "run2", runB = "run0", mergeName = mergeName, dataPath = ".", fileInfo, features,
#'  mzPntrs, prec2chromIndex, precursors, params, adaptiveRTs, refRuns, multipeptide, peptideScores, ropenms = NULL)
#' rm(mzPntrs)
#' file.remove(file.path(".", "mzml", paste0(mergeName, ".chrom.sqMass")))
#' file.remove(list.files(".", pattern = "*_av.rds", full.names = TRUE))
#' }
getNodeRun <- function(runA, runB, mergeName, dataPath, fileInfo, features, mzPntrs, prec2chromIndex,
                       precursors, params, adaptiveRTs, refRuns, multipeptide, peptideScores, ropenms, applyFun = lapply){
  peptides <- unique(precursors$peptide_id)

  ##### Select reference for each peptide and update peptideScores. #####
  var1 <- rep(NA_integer_, length(peptides))
  var2 <- rep(NA_integer_, length(peptides))
  for(i in seq_along(peptides)){
    temp <- peptideScores[[i]][.(c(runA, runB)),]
    pvalA <- .subset2(temp, "pvalue")[[1]]
    pvalB <- .subset2(temp, "pvalue")[[1]]
    if(is.na(pvalA)) pvalA <- 1
    if(is.na(pvalB)) pvalB <- 1
    var1[i] <- ifelse(pvalA > pvalB, 2L, 1L)
    if(nrow(temp) > 0){
      idx <- which(peptideScores[[i]][["run"]] == mergeName)
      set(peptideScores[[i]], idx, c(2L, 3L, 4L),
          list(max(temp$score), min(temp$pvalue), min(temp$qvalue)))
    }
    temp <- multipeptide[[i]][.(c(runA, runB)),]
    idx <- which.min(temp$m_score)
    idx <- ifelse(length(idx)==0, 1, idx)
    var2[i] <- .subset2(temp, "transition_group_id")[[idx]]
  }
  refRun <- data.table("var1" = var1, "var2" = as.character(var2))

  ##### Get childXICs #####
  message("Getting merged chromatograms for run ", mergeName)
  mergedXICs_alignedVec <- getChildXICs(runA, runB, fileInfo, features, mzPntrs, precursors, prec2chromIndex,
                                        refRun, peptideScores, params, applyFun)
  mergedXICs <- mergedXICs_alignedVec[[1]]
  alignedVecs <- mergedXICs_alignedVec[[2]]
  adaptiveRT <- mergedXICs_alignedVec[[3]]

  ##### Get merged features, calculate intensities, left width, right width, m-score. #####
  # we can also run pyopenms feature finder on new chromatogram.
  message("Getting merged features for run ", mergeName)
  childFeatures <- applyFun(seq_along(peptides), function(i){
    analytes <- as.integer(names(mergedXICs[[i]]))
    alignedVec <- alignedVecs[[i]]
    childFeature <- data.table()
    for(j in seq_along(mergedXICs[[i]])){
      if(is.null(mergedXICs[[i]][[j]])) next
      XICs <- mergedXICs[[i]][[j]]
      analyte <- analytes[j]
      df.A <- features[[runA]][.(analyte),]
      df.B <- features[[runB]][.(analyte),]
      if(refRun[i, 1] == 1L) {
        df.ref <- df.A
        df.eXp <- df.B
      } else{
        df.ref <- df.B
        df.eXp <- df.A
      }
      if((nrow(df.ref) + nrow(df.eXp)) == 0) next
      rows <- getChildFeature(XICs, alignedVec, df.ref, df.eXp, params)
      childFeature <- dplyr::bind_rows(childFeature, rows)
    }
    childFeature
  })

  ##### Add features to multipeptide #####
  newMP <- applyFun(seq_along(peptides), function(i){
    peptide_chr <- as.character(peptides[i])
    temp <- multipeptide[[i]]
    df <- childFeatures[[i]]
    if(nrow(df) == 0) df <- data.frame("transition_group_id" = temp$transition_group_id[1],
                                       "feature_id" = bit64::NA_integer64_,
                                       "RT" = NA_real_, "intensity" = NA_real_,
                                       "leftWidth" = NA_real_, "rightWidth" = NA_real_,
                                       "peak_group_rank" = NA_integer_, "m_score" = NA_real_, stringsAsFactors = FALSE)
    df$run <- mergeName
    df$alignment_rank <- NA_integer_
    dplyr::bind_rows(temp, df)
  })
  names(newMP) <- names(multipeptide)

  rm(temp)
  ##### Add node features #####
  childFeatures <- dplyr::bind_rows(childFeatures)
  assign("temp", childFeatures, envir = parent.frame(n = 1))
  with(parent.frame(n = 1), features[[mergeName]] <- temp)

  ##### Write node mzML file #####
  mergedXICs <- unlist(mergedXICs, recursive = FALSE, use.names = FALSE)
  if(params[["chromFile"]] =="mzML"){
    fileName <- file.path(dataPath, "mzml", paste0(mergeName, ".chrom.mzML"))
    createMZML(ropenms, fileName, mergedXICs, precursors$transition_ids)
  } else if(params[["chromFile"]] =="sqMass"){
    fileName <- file.path(dataPath, "mzml", paste0(mergeName, ".chrom.sqMass"))
    createSqMass(fileName, mergedXICs, precursors$transition_ids, params[["lossy"]])
  }

  ##### Add node run to fileInfo #####
  row <- data.frame(runName = mergeName, spectraFile = NA_character_, spectraFileID = NA_character_,
                    featureFile = NA_character_, chromatogramFile = fileName, row.names = mergeName)
  df <- dplyr::bind_rows(fileInfo, row)
  assign("temp", df, envir = parent.frame(n = 1))
  with(parent.frame(n = 1), fileInfo <- temp)

  ##### Add node run to mzPntrs #####
  if(params[["chromFile"]] =="mzML"){
    mzPntr <- mzR::openMSfile(fileName, backend = "pwiz")
    chromHead <- mzR::chromatogramHeader(mzPntr) #TODO: Make sure that chromatogramIndex is read as integer64
  } else if(params[["chromFile"]] =="sqMass"){
    mzPntr <- DBI::dbConnect(RSQLite::SQLite(), dbname = fileName)
    chromHead <- readSqMassHeader(mzPntr)
  }
  assign("temp", mzPntr, envir = parent.frame(n = 1))
  with(parent.frame(n = 1), mzPntrs[[mergeName]] <- temp)

  ##### Add node run to prec2chromIndex #####
  prec2transition <- dplyr::select(precursors, .data$transition_group_id, .data$transition_ids) %>%
    tidyr::unnest(.data$transition_ids) %>% as.data.frame()
  chromatogramIdAsInteger(chromHead) # Select only chromatogramId, chromatogramIndex
  df <- mapPrecursorToChromIndices(prec2transition, chromHead) # Get chromatogram Index for each precursor.
  df <- df[match(precursors$transition_group_id, df$transition_group_id),]
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
  newMP
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
  timeParent <- alignedVec[, c(1L, 3L)]
  colnames(timeParent) <- c("tAligned", "alignedChildTime")
  if(nrow(df.ref)!=0) df.ref <- trfrParentFeature(XICs, timeParent, df.ref, params)

  # Convert eXp features to childXIC
  timeParent <-alignedVec[, c(2L, 3L)]
  colnames(timeParent) <- c("tAligned", "alignedChildTime")
  if(nrow(df.eXp)!=0) df.eXp <- trfrParentFeature(XICs, timeParent, df.eXp, params)

  # Assemble converted features. Pick non-overlapping features based on minimum m-score.
  df <- dplyr::bind_rows(df.ref, df.eXp)
  if(nrow(df) == 0) return(NULL)
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
#' precursors <- dplyr::arrange(precursors, .data$peptide_id, .data$transition_group_id)
#' peptideIDs <- unique(precursors$peptide_id)
#' peptideScores <- getPeptideScores(fileInfo, peptideIDs, oswMerged = TRUE, params[["runType"]], params[["context"]])
#' peptideScores <- lapply(peptideIDs, function(pep) dplyr::filter(peptideScores, .data$peptide_id == pep))
#' names(peptideScores) <- as.character(peptideIDs)
#' prec2chromIndex <- getChromatogramIndices(fileInfo, precursors, mzPntrs)
#' var2 <- as.character(sapply(peptideIDs, function(p) precursors$transition_group_id[which(precursors$peptide_id == p)[1]]))
#' refRun <- data.frame(rep(1L, length(peptideIDs)), var2)
#' mergedXICs <- getChildXICs(runA = "run1", runB = "run2", fileInfo, features, mzPntrs,
#'   precursors, prec2chromIndex, refRun, peptideScores, params)
#' rm(mzPntrs)
#' @export
getChildXICs <- function(runA, runB, fileInfo, features, mzPntrs, precursors, prec2chromIndex, refRun,
                         peptideScores, params, applyFun = lapply){
  peptides <- unique(precursors$peptide_id)

  #### Get global alignment between runs ####
  pair <- paste(runA, runB, sep = "_")
  globalFit1 <- getGlobalAlignment(features, runA, runB,
                                   params[["globalAlignment"]], params[["globalAlignmentFdr"]], params[["globalAlignmentSpan"]])
  adaptiveRT1 <-  params[["RSEdistFactor"]]*getRSE(globalFit1, params[["globalAlignment"]])
  globalFit1 <- extractFit(globalFit1, params[["globalAlignment"]])
  globalFit2 <- getGlobalAlignment(features, runB, runA,
                                   params[["globalAlignment"]], params[["globalAlignmentFdr"]], params[["globalAlignmentSpan"]])
  adaptiveRT2 <-  params[["RSEdistFactor"]]*getRSE(globalFit2, params[["globalAlignment"]])
  globalFit2 <- extractFit(globalFit2, params[["globalAlignment"]])

  #### Get merged XICs ####
  num_of_batch <- ceiling(length(peptides)/params[["batchSize"]])
  temp <- lapply(1:num_of_batch, parFUN1, runA, runB, peptides, precursors, prec2chromIndex, mzPntrs,
                 params, peptideScores, refRun, globalFit1, globalFit2, adaptiveRT1, adaptiveRT2, applyFun)
  temp <- unlist(temp, recursive = FALSE)
  mergedXICs <- lapply(temp, `[[`, 1)
  alignedVecs <- lapply(temp, `[[`, 2)
  list(mergedXICs, alignedVecs, adaptiveRT = list(ab= adaptiveRT1, ba = adaptiveRT2))
}

parFUN1 <- function(iBatch, runA, runB, peptides, precursors, prec2chromIndex, mzPntrs, params,
                    peptideScores, refRun, globalFit1, globalFit2, adaptiveRT1, adaptiveRT2, applyFun){
  batchSize <- params[["batchSize"]]
  if(params[["chromFile"]] =="mzML") fetchXIC = extractXIC_group
  if(params[["chromFile"]] =="sqMass") fetchXIC = extractXIC_group2
  strt <- ((iBatch-1)*batchSize+1)
  stp <- min((iBatch*batchSize), length(peptides))
  ##### Get XICs for the batch from both runs #####
  XICs <- lapply(strt:stp, function(rownum){
    ##### Get transition_group_id for that peptideID #####
    idx <- which(precursors$peptide_id == peptides[rownum])
    analytes <- .subset2(precursors, "transition_group_id")[idx]
    ##### Get XIC_group from runA and runB. If missing, add NULL #####
    chromIndices.A <- prec2chromIndex[[runA]][["chromatogramIndex"]][idx]
    chromIndices.B <- prec2chromIndex[[runB]][["chromatogramIndex"]][idx]
    nope <- any(is.na(c(unlist(chromIndices.A), unlist(chromIndices.B))))
    nope <- nope | is.null(unlist(chromIndices.A)) | is.null(unlist(chromIndices.B))
    if(nope) return(NULL)
    XICs.A <- lapply(chromIndices.A, function(i1) fetchXIC(mzPntrs[[runA]], i1))
    XICs.B <- lapply(chromIndices.B, function(i1) fetchXIC(mzPntrs[[runB]], i1))
    names(XICs.A) <- names(XICs.B) <- as.character(analytes)
    list(XICs.A, XICs.B)
  })

  ##### Get child XICs for the batch from both runs #####
  cluster <- applyFun(strt:stp, function(rownum){
    peptide <- peptides[rownum]
    idx <- (rownum - (iBatch-1)*batchSize)
    if(is.null(XICs[[idx]])){
      warning("Chromatogram indices for ", peptide, " are missing.")
      message("Skipping peptide ", peptide, ".")
      analytes <- precursors[precursors$peptide_id == peptide, "transition_group_id"]
      return(list(vector(mode = "list", length = length(analytes)), NULL))
    } else {
      XICs.A <- XICs[[idx]][[1]]
      XICs.B <- XICs[[idx]][[2]]
    }

    ##### Calculate the weights of XICs from runA and runB #####
    temp <- peptideScores[[as.character(peptide)]]
    pA <- temp$pvalue[temp$run == runA]
    pB <- temp$pvalue[temp$run == runB]
    wA <- ifelse(length(pA)==0, 0.3, -log10(pA))
    wB <- ifelse(length(pB)==0, 0.3, -log10(pB))
    wA <- wA/(wA + wB)

    ##### Decide the reference run from runA and runB. #####
    if(refRun[rownum, 1] == 1L){
      XICs.ref <- XICs.A
      XICs.eXp <- XICs.B
      globalFit <- globalFit1
      adaptiveRT <- adaptiveRT1
      wRef <- wA
    } else{
      XICs.ref <- XICs.B
      XICs.eXp <- XICs.A
      globalFit <- globalFit2
      adaptiveRT <- adaptiveRT2
      wRef <- (1-wA)
    }

    ##### Select 1) all precursors OR 2) high quality precursor. #####
    analytes_chr <- names(XICs.A)
    analyte_chr <- .subset2(refRun, 2L)[[rownum]]
    if(FALSE){
      # Turned off as precursor XICs have different time ranges.
      XICs.ref.pep <- unlist(XICs.ref, recursive = FALSE, use.names = FALSE)
      XICs.eXp.pep <- unlist(XICs.eXp, recursive = FALSE, use.names = FALSE)
    } else {
      XICs.ref.pep <- XICs.ref[[analyte_chr]]
      XICs.eXp.pep <- XICs.eXp[[analyte_chr]]
    }

    B1p <- getPredict(globalFit, XICs.ref.pep[[1]][1,1], params[["globalAlignment"]])
    len <- nrow(XICs.ref.pep[[1]])
    B2p <- getPredict(globalFit, XICs.ref.pep[[1]][len,1], params[["globalAlignment"]])
    #### Merge chromatograms  ####
    merged_xics <- getChildXICpp(XICs.ref.pep, XICs.eXp.pep, params[["kernelLen"]], params[["polyOrd"]],
                  params[["alignType"]], adaptiveRT, params[["normalization"]],
                  params[["simMeasure"]], B1p, B2p, params[["goFactor"]], params[["geFactor"]],
                  params[["cosAngleThresh"]], params[["OverlapAlignment"]],
                  params[["dotProdThresh"]], params[["gapQuantile"]], params[["kerLen"]],
                  params[["hardConstrain"]], params[["samples4gradient"]], wRef,
                  params[["splineMethod"]], params[["mergeTime"]], params[["keepFlanks"]])
    merged_xics[[1]] <- list(merged_xics[[1]])
    names(merged_xics[[1]]) <- analyte_chr
    otherPrecs <- setdiff(analytes_chr, analyte_chr)
    if(length(otherPrecs) !=0){
      for(name in otherPrecs){
        merged_xics[[1]][[name]] <- otherChildXICpp(XICs.ref[[name]], XICs.eXp[[name]], 0L,
                  params[["polyOrd"]], merged_xics[[2]], merged_xics[[1]][[1]][[1]][,1],
                  wRef, params[["splineMethod"]])
      }
    }
    merged_xics[[1]] <- merged_xics[[1]][match(analytes_chr, names(merged_xics[[1]]))]
    merged_xics # 1st element has list of precursors. 2nd element has aligned time vectors.
  })
  cluster
}
