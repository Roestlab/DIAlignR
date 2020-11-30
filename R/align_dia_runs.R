#' Outputs intensities for each analyte from aligned Targeted-MS runs
#'
#' This function expects osw and mzml directories at dataPath. It first reads osw files and fetches chromatogram indices for each analyte.
#' It then align XICs of its reference XICs. Best peak, which has lowest m-score, about the aligned retention time is picked for quantification.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-14
#' @importFrom dplyr %>%
#' @inheritParams checkParams
#' @param dataPath (string) path to mzml and osw directory.
#' @param outFile (string) name of the output file.
#' @param oswMerged (logical) TRUE for experiment-wide FDR and FALSE for run-specific FDR by pyprophet.
#' @param runs (a vector of string) names of mzml file without extension.
#' @param refRun (string) reference for alignment. If no run is provided, m-score is used to select reference run.
#' @param applyFun (function) value must be either lapply or BiocParallel::bplapply.
#' @return An output table with following columns: precursor, run, intensity, RT, leftWidth, rightWidth,
#'  peak_group_rank, m_score, alignment_rank, peptide_id, sequence, charge, group_label.
#'
#' @seealso \code{\link{getRunNames}, \link{getFeatures}, \link{setAlignmentRank}, \link{getMultipeptide}}
#' @examples
#' params <- paramsDIAlignR()
#' params[["context"]] <- "experiment-wide"
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' BiocParallel::register(BiocParallel::MulticoreParam(workers = 4, progressbar = TRUE))
#' alignTargetedRuns(dataPath, outFile = "testDIAlignR.tsv", params = params, applyFun = BiocParallel::bplapply)
#' @references Gupta S, Ahadi S, Zhou W, Röst H. "DIAlignR Provides Precise Retention Time Alignment Across Distant Runs in DIA and Targeted Proteomics." Mol Cell Proteomics. 2019 Apr;18(4):806-817. doi: https://doi.org/10.1074/mcp.TIR118.001132 Epub 2019 Jan 31.
#'
#' @export
alignTargetedRuns <- function(dataPath, outFile = "DIAlignR.tsv", params = paramsDIAlignR(), oswMerged = TRUE, runs = NULL,
                              refRun = NULL, applyFun = lapply){
  #### Check if all parameters make sense.  #########
  checkParams(params)

  #### Get filenames from .osw file and check consistency between osw and mzML files. #################
  fileInfo <- getRunNames(dataPath, oswMerged)
  fileInfo <- updateFileInfo(fileInfo, runs)
  runs <- rownames(fileInfo)
  message("Following runs will be aligned:")
  print(fileInfo[, "runName"], sep = "\n")

  #### Get Precursors from the query and respectve chromatogram indices. ######
  # Get all the precursor IDs, transition IDs, Peptide IDs, Peptide Sequence Modified, Charge.
  precursors <- getPrecursors(fileInfo, oswMerged, params[["runType"]], params[["context"]], params[["maxPeptideFdr"]], params[["level"]])
  precursors <- dplyr::arrange(precursors, .data$peptide_id, .data$transition_group_id)

  #### Get Peptide scores, pvalue and qvalues. ######
  # Some peptides may not be found due to using a subset of runs. Appends NA for them.
  # This translates as "Chromatogram indices for peptide ID are missing in NA"
  peptideIDs <- unique(precursors$peptide_id)
  peptideScores <- getPeptideScores(fileInfo, peptideIDs, oswMerged, params[["runType"]], params[["context"]])
  peptideScores <- lapply(peptideIDs, function(pep) dplyr::filter(peptideScores, .data$peptide_id == pep))
  names(peptideScores) <- as.character(peptideIDs)

  #### Get reference run for each precursor ########
  idx <- which(fileInfo$runName == refRun)
  if(length(idx) == 0){
    message("Calculating reference run for each peptide.")
    refRuns <- getRefRun(peptideScores, applyFun)
  } else{
    run <- rownames(fileInfo)[idx]
    refRuns <- data.frame("peptide_id" = peptideIDs, "run" = run)
  }

  #### Get OpenSWATH peak-groups and their retention times. ##########
  if(params[["transitionIntensity"]]){
    features <- getTransitions(fileInfo, params[["maxFdrQuery"]], params[["runType"]], applyFun)
  } else{
    features <- getFeatures(fileInfo, params[["maxFdrQuery"]], params[["runType"]], applyFun)
  }

  #### Collect pointers for each mzML file. #######
  message("Collecting metadata from mzML files.")
  mzPntrs <- getMZMLpointers(fileInfo)
  message("Metadata is collected from mzML files.")

  #### Get chromatogram Indices of precursors across all runs. ############
  message("Collecting chromatogram indices for all precursors.")
  prec2chromIndex <- getChromatogramIndices(fileInfo, precursors, mzPntrs, applyFun)

  #### Convert features into multi-peptide #####
  message("Building multipeptide.")
  if(params[["transitionIntensity"]]){
    multipeptide <- getMultipeptide2(precursors, features, applyFun)
  } else{
    multipeptide <- getMultipeptide(precursors, features, applyFun)
  }
  message(length(multipeptide), " peptides are in the multipeptide.")

  #### Container to save Global alignments.  #######
  message("Calculating global alignments.")
  globalFits <- getGlobalFits(refRuns, features, fileInfo, params[["globalAlignment"]],
                              params[["globalAlignmentFdr"]], params[["globalAlignmentSpan"]], applyFun)
  RSE <- applyFun(globalFits, getRSE)

  # TODO: Check dimensions of multipeptide, PeptideIDs, precursors etc makes sense.
  #### Perform pairwise alignment ###########
  message("Performing reference-based alignment.")
  start_time <- Sys.time()
  num_of_batch <- ceiling(length(multipeptide)/params[["batchSize"]])
  multipeptide <- lapply(1:num_of_batch, perBatch, peptideIDs, multipeptide, refRuns, precursors,
                           prec2chromIndex, fileInfo, mzPntrs, params, globalFits, RSE, applyFun)
  multipeptide <- unlist(multipeptide, recursive = FALSE)
  names(multipeptide) <- as.character(peptideIDs)

  #### Cleanup.  #######
  rm(mzPntrs)
  end_time <- Sys.time() # Report the execution time for hybrid alignment step.
  message("The execution time for alignment:")
  print(end_time - start_time)

  #### Write tables to the disk  #######
  finalTbl <- writeTables(fileInfo, multipeptide, precursors)
  if(params[["transitionIntensity"]]){
    finalTbl$intensity <- lapply(finalTbl$intensity,round, 3)
    finalTbl$intensity <- sapply(finalTbl$intensity, function(x) paste(unlist(x), collapse=", "))
  }
  utils::write.table(finalTbl, file = outFile, sep = "\t", row.names = FALSE, quote = FALSE)
  message("Retention time alignment across runs is done.")
  message(paste0(outFile, " file has been written."))

  #### Write alignment summary  #######
  alignmentStats(finalTbl, params)
}

#' AlignObj for analytes between a pair of runs
#'
#' This function expects osw and mzml directories at dataPath. It first reads osw files and fetches chromatogram indices for each requested analyte.
#' It then align XICs of each analyte to its reference XICs. AlignObj is returned which contains aligned indices and cumulative score along the alignment path.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-14
#' @importFrom rlang .data
#' @inheritParams alignTargetedRuns
#' @param analytes (vector of integers) transition_group_ids for which features are to be extracted.
#' @param objType (char) Must be selected from light, medium and heavy.
#' @return A list of fileInfo and AlignObjs. Each AlignObj is an S4 object. Three most-important slots are:
#' \item{indexA_aligned}{(integer) aligned indices of reference run.}
#' \item{indexB_aligned}{(integer) aligned indices of experiment run.}
#' \item{score}{(numeric) cumulative score of alignment.}
#' @seealso \code{\link{plotAlignedAnalytes}, \link{getRunNames}, \link{getFeatures}, \link{getXICs4AlignObj}, \link{getAlignObj}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' params <- paramsDIAlignR()
#' params[["context"]] <- "experiment-wide"
#' runs <- c("hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt",
#'  "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
#'  "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
#' analytes <- c(32L, 898L, 2474L)
#' AlignObjOutput <- getAlignObjs(analytes, runs, dataPath = dataPath)
#' plotAlignedAnalytes(AlignObjOutput)
#'
#' @references Gupta S, Ahadi S, Zhou W, Röst H. "DIAlignR Provides Precise Retention Time Alignment Across Distant Runs in DIA and Targeted Proteomics." Mol Cell Proteomics. 2019 Apr;18(4):806-817. doi: https://doi.org/10.1074/mcp.TIR118.001132 Epub 2019 Jan 31.
#'
#' @export
getAlignObjs <- function(analytes, runs, dataPath = ".", refRun = NULL, oswMerged = TRUE,
                         params = paramsDIAlignR(), objType = "light"){
  #### Check if all parameters make sense.  #########
  checkParams(params)

  ##### Get filenames from osw files and check if names are consistent between osw and mzML files. ######
  filenames <- getRunNames(dataPath, oswMerged)
  filenames <- updateFileInfo(filenames, runs)
  missingRun <- setdiff(runs, filenames$runName)
  if(length(missingRun) != 0){
    return(stop(missingRun, " runs are not found."))
  }
  message("Following runs will be aligned:")
  print(filenames[, "runName"], sep = "\n")

  ######### Collect pointers for each mzML file. #######
  message("Collecting metadata from mzML files.")
  mzPntrs <- getMZMLpointers(filenames)
  message("Metadata is collected from mzML files.")

  ######### Get Precursors from the query and respectve chromatogram indices. ######
  precursors <- getPrecursorByID(analytes, filenames)

  #### Precursors for which features are identified. ##############
  features <- getFeatures(filenames, params[["maxFdrQuery"]], params[["runType"]])

  ###### Report analytes that are not found ########
  refAnalytes <- analytesFromFeatures(features, analyteFDR = params[["analyteFDR"]], commonAnalytes = FALSE)
  analytesFound <- intersect(analytes, refAnalytes)
  analytesNotFound <- setdiff(analytes, analytesFound)
  if(length(analytesNotFound)>0){
    message(paste(analytesNotFound, "not found with FDR cut-off."))
  }
  analytes <- analytesFound
  precursors <- precursors[precursors[["transition_group_id"]] %in% analytes, ]
  if(nrow(precursors) == 0){
    stop("No precursors are found below ", params[["analyteFDR"]])
  }

  ############# Get chromatogram Indices of precursors across all runs. ############
  prec2chromIndex <- getChromatogramIndices(filenames, precursors, mzPntrs)

  #### Get Peptide scores, pvalue and qvalues. ######
  peptideIDs <- unique(precursors$peptide_id)
  peptideScores <- getPeptideScores(filenames, peptideIDs, oswMerged, params[["runType"]], params[["context"]])
  peptideScores <- lapply(peptideIDs, function(pep) dplyr::filter(peptideScores, .data$peptide_id == pep))
  names(peptideScores) <- as.character(peptideIDs)

  ############## Get reference run for each precursor ########
  idx <- which(filenames$runName == refRun)
  if(length(idx) == 0){
    print("Finding reference run using SCORE_PEPTIDE table")
    refRun <- data.frame("transition_group_id" = precursors$transition_group_id,
                         "run" = NA_character_)
    temp <- getRefRun(peptideScores)
    refRun$run <- temp$run[match(precursors$peptide_id, temp$peptide_id)]
  } else{
    run <- rownames(filenames)[idx]
    refRun <- data.frame("transition_group_id" = precursors$transition_group_id,
                            "run" = run)
  }

  ####################### Get XICs ##########################################
  # Get Chromatogram for each peptide in each run.
  message("Fetching Extracted-ion chromatograms from runs")
  XICs <- getXICs4AlignObj(mzPntrs, filenames, filenames[, "runName"], prec2chromIndex, analytes)
  rm(mzPntrs)

  ####################### Perfrom alignment ##########################################
  AlignObjs <- vector("list", length(analytes))
  globalFits <- list()
  RSE <- list()
  runs <- rownames(filenames)
  message("Perfroming alignment")
  for(analyteIdx in seq_along(analytes)){
    analyte <- as.character(analytes[analyteIdx])
    ref <- refRun[["run"]][analyteIdx]
    AlignObjs[[analyteIdx]] <- list()

    # Get XIC_group from reference run
    XICs.ref <- XICs[[filenames[ref,"runName"]]][[analyte]]
    if(is.null(XICs.ref)){
      warning("Chromatogram indices for ", analyte, " are missing in ", filenames[ref, "runName"])
      message("Skipping ", analyte)
      AlignObjs[[analyteIdx]] <- NULL
      next
    }
    XICs.ref.s <- smoothXICs(XICs.ref, type = params[["XICfilter"]], kernelLen = params[["kernelLen"]],
                             polyOrd = params[["polyOrd"]])
    exps <- setdiff(runs, ref)

    # Align experiment run to reference run
    for(eXp in exps){
      pair <- paste(ref, eXp, sep = "_")
      AlignObjs[[analyteIdx]][[pair]] <- list()
      # Get XIC_group from experiment run
      XICs.eXp <- XICs[[filenames[eXp,"runName"]]][[analyte]]
      if(is.null(XICs.eXp)){
        warning("Chromatogram indices for ", analyte, " are missing in ", filenames[eXp, "runName"])
        message("Skipping ", analyte)
        AlignObjs[[analyteIdx]][[pair]] <- NULL
        next
      }
      XICs.eXp.s <- smoothXICs(XICs.eXp, type = params[["XICfilter"]], kernelLen = params[["kernelLen"]],
                               polyOrd = params[["polyOrd"]])
      # Get the loess fit for hybrid alignment
      if(any(pair %in% names(globalFits))){
        globalFit <- globalFits[[pair]]
      } else{
        globalFit <- getGlobalAlignment(features, ref, eXp, params[["globalAlignment"]],
                                        params[["globalAlignmentFdr"]], params[["globalAlignmentSpan"]])
        globalFits[[pair]] <- globalFit
        RSE[[pair]] <- getRSE(globalFit)
      }
      adaptiveRT <- params[["RSEdistFactor"]]*RSE[[pair]]

      # Fetch alignment object between XICs.ref and XICs.eXp
      AlignObj <- getAlignObj(XICs.ref.s, XICs.eXp.s, globalFit, params[["alignType"]], adaptiveRT,
                              params[["normalization"]], params[["simMeasure"]], params[["goFactor"]],
                              params[["geFactor"]], params[["cosAngleThresh"]], params[["OverlapAlignment"]],
                              params[["dotProdThresh"]], params[["gapQuantile"]], params[["kerLen"]],
                              params[["hardConstrain"]], params[["samples4gradient"]],
                              objType)
      # Attach AlignObj for the analyte.
      AlignObjs[[analyteIdx]][[pair]][["AlignObj"]] <- AlignObj
      # Attach intensities of reference XICs.
      AlignObjs[[analyteIdx]][[pair]][["ref"]] <- XICs.ref
      # Attach intensities of experiment XICs.
      AlignObjs[[analyteIdx]][[pair]][["eXp"]] <- XICs.eXp
      # Attach peak boundaries to the object.
      AlignObjs[[analyteIdx]][[pair]][["peak"]] <- features[[ref]] %>%
        dplyr::filter(.data$transition_group_id == as.integer(analyte) & .data$peak_group_rank == 1) %>%
        dplyr::select(.data$leftWidth, .data$RT, .data$rightWidth) %>%
        as.vector()
    }
  }
  names(AlignObjs) <- as.character(analytes)

  ####################### Return AlignedObjs ##########################################
  message("Alignment done. Returning AlignObjs")
  list(filenames, AlignObjs)
}


#' Aligns an analyte from an experiment to the reference run
#'
#' df contains unaligned features for an analyte across multiple runs. This function aligns eXp run to
#' ref run and updates corresponding features.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-07-26
#' @keywords internal
#' @inheritParams alignIthAnalyte
#' @param eXp (string) name of the run to be aligned to reference run. Must be in the rownames of fileInfo.
#' @param ref (string) name of the reference run. Must be in the rownames of fileInfo.
#' @param analyte_chr (string) Precursor ID of the requested analyte.
#' @param XICs.ref.s (list of dataframes) Smoothed fragment-ion chromatograms of the analyte_chr from the reference run.
#' @param df (dataframe) a collection of features related to analyte_chr.
#' @return (dataframe) aligned features of analyte_chr in eXp run.
#' @seealso \code{\link{alignTargetedRuns}, \link{alignIthAnalyte}, \link{setAlignmentRank}, \link{getMultipeptide}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
alignToRef <- function(eXp, ref, preIdx, analytes, fileInfo, XICs.ref.s, params, prec2chromIndex,
                       mzPntrs, df, globalFits, RSE){
  # Get XIC_group from experiment run. if missing, go to next run.
  chromIndices <- prec2chromIndex[[eXp]][["chromatogramIndex"]][preIdx]
  if(any(is.na(unlist(chromIndices))) | is.null(unlist(chromIndices))){
    message("Chromatogram indices for precursor ", analytes, " are missing in ", fileInfo[eXp, "runName"])
    message("Skipping precursor ", analytes, " in ", fileInfo[eXp, "runName"], ".")
    df.eXp <- df[df[["run"]] == eXp, ]
    return(df.eXp)
  } else {
    XICs.eXp <- lapply(chromIndices, function(i) extractXIC_group(mz = mzPntrs[[eXp]], chromIndices = i))
    XICs.eXp.s <- lapply(XICs.eXp, smoothXICs, type = params[["XICfilter"]], kernelLen = params[["kernelLen"]],
                         polyOrd = params[["polyOrd"]])
    names(XICs.eXp.s) <- names(XICs.eXp) <- as.character(analytes)
  }
  if(params[["smoothPeakArea"]]) XICs.eXp <- XICs.eXp.s

  # Select 1) all precursors OR 2) high quality precursor
  if(FALSE){
    # Turned off as precursor XICs have different time ranges.
    XICs.ref.pep <- unlist(XICs.ref.s, recursive = FALSE, use.names = FALSE)
    XICs.eXp.pep <- unlist(XICs.eXp.s, recursive = FALSE, use.names = FALSE)
  } else {
    temp <- df[df$run == ref, c("transition_group_id", "m_score")]
    analyte_chr <- as.character(temp[which.min(temp$m_score), "transition_group_id"])
    XICs.ref.pep <- XICs.ref.s[[analyte_chr]]
    XICs.eXp.pep <- XICs.eXp.s[[analyte_chr]]
  }

  ##### Get the aligned Indices #####
  pair <- paste(ref, eXp, sep = "_")
  globalFit <- globalFits[[pair]]
  adaptiveRT <- params[["RSEdistFactor"]]*RSE[[pair]]
  tAligned <- getAlignedTimes(XICs.ref.pep, XICs.eXp.pep, globalFit, params[["alignType"]], adaptiveRT,
                               params[["normalization"]], params[["simMeasure"]], params[["goFactor"]],
                               params[["geFactor"]], params[["cosAngleThresh"]], params[["OverlapAlignment"]],
                               params[["dotProdThresh"]], params[["gapQuantile"]], params[["kerLen"]],
                               params[["hardConstrain"]], params[["samples4gradient"]], objType = "light")
  df.eXp <- setAlignmentRank(df, ref, eXp, tAligned, XICs.eXp, params, adaptiveRT)
  df.eXp <- setOtherPrecursors(df.eXp, XICs.eXp, analytes, params)
  if(params[["recalIntensity"]]) df.eXp <- reIntensity(df.eXp, XICs.eXp, params)
  df.eXp
}

#' Aligns an analyte across runs
#'
#' For the ith analyte in multipeptide, this function aligns all runs to the reference run. The result is
#' a dataframe that contains aligned features corresponding to the analyte across all runs.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-07-26
#' @keywords internal
#' @import dplyr
#' @inheritParams checkParams
#' @param rownum (integer) represnts the index of the multipepetide to be aligned.
#' @param peptideIDs (integer) vector of peptideIDs.
#' @param multipeptide (list) contains multiple data-frames that are collection of features
#'  associated with analytes. This is an output of \code{\link{getMultipeptide}}.
#' @param refRuns (data-frame) output of \code{\link{getRefRun}}. Must have two columsn : transition_group_id and run.
#' @param precursors (data-frame) atleast two columns transition_group_id and transition_ids are required.
#' @param prec2chromIndex (list) a list of dataframes having following columns: \cr
#' transition_group_id: it is PRECURSOR.ID from osw file. \cr
#' chromatogramIndex: index of chromatogram in mzML file.
#' @param fileInfo (data-frame) output of \code{\link{getRunNames}}.
#' @param mzPntrs (list) a list of mzRpwiz.
#' @param globalFits (list) each element is either of class lm or loess. This is an output of \code{\link{getGlobalFits}}.
#' @param RSE (list) Each element represents Residual Standard Error of corresponding fit in globalFits.
#' @return (dataframe) a collection of aligned features for a peptide( = peptideIDs[rownum]).
#' @seealso \code{\link{alignTargetedRuns}, \link{alignToRef}, \link{getAlignedTimes}, \link{getMultipeptide}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
alignIthAnalyte <- function(rownum, peptideIDs, multipeptide, refRuns, precursors, prec2chromIndex,
                            fileInfo, mzPntrs, params, globalFits, RSE){
  ##### Get peptideID, its reference run and multipeptide #####
  # print(rownum)
  peptide <- peptideIDs[rownum]
  df <- multipeptide[[rownum]]
  ref <- refRuns[["run"]][rownum]

  ##### Get transition_group_id for that peptideID #####
  idx <- which(precursors$peptide_id == peptide)
  analytes <- precursors[idx, "transition_group_id"]

  ##### Get XIC_group from reference run. if missing, return unaligned features #####
  chromIndices <- prec2chromIndex[[ref]][["chromatogramIndex"]][idx]
  if(any(is.na(unlist(chromIndices))) | is.null(unlist(chromIndices))){
    message("Chromatogram indices for peptide ", peptide, " are missing in ", fileInfo[ref, "runName"])
    message("Skipping peptide ", peptide, " across all runs.")
    return(df)
  } else {
    XICs.ref <- lapply(chromIndices, function(i) extractXIC_group(mz = mzPntrs[[ref]], chromIndices = i))
    XICs.ref.s <- lapply(XICs.ref, smoothXICs, type = params[["XICfilter"]], kernelLen = params[["kernelLen"]],
                             polyOrd = params[["polyOrd"]])
    names(XICs.ref.s) <- names(XICs.ref) <- as.character(analytes)
  }
  if(params[["smoothPeakArea"]]) XICs.ref <- XICs.ref.s

  ##### Set alignment rank for all precrusors of the peptide in the reference run #####
  refIdx <- which(df[["run"]] == ref & df[["peak_group_rank"]] == 1)
  refIdx <- refIdx[which.min(df$m_score[refIdx])]
  if(length(refIdx)==0) {
    message("Features for peptide ", peptide, " is missing in ", fileInfo[ref, "runName"])
    message("Skipping peptide ", peptide, " across all runs.")
    return(df)
  }
  df[["alignment_rank"]][refIdx] <- 1L
  df.ref <- df[df$run == ref,]
  df.ref <- setOtherPrecursors(df.ref, XICs.ref, analytes, params)
  # Update multipeptide reference intensity if recal is true
  if(params[["recalIntensity"]]) df.ref <- reIntensity(df.ref, XICs.ref, params)

  ##### Align all runs to reference run and set their alignment rank #####
  exps <- setdiff(rownames(fileInfo), ref)
  df.exps <- lapply(exps, alignToRef, ref, idx, analytes, fileInfo, XICs.ref.s, params, prec2chromIndex, mzPntrs,
                    df, globalFits, RSE)

  ##### Return the dataframe with alignment rank set to TRUE #####
  updateOnalignTargetedRuns(rownum)
  newDF <- dplyr::bind_rows(df.ref, df.exps)
  newDF
}



perBatch <- function(iBatch, peptideIDs, multipeptide, refRuns, precursors, prec2chromIndex,
                     fileInfo, mzPntrs, params, globalFits, RSE, applyFun = lapply){
  message("Processing Batch ", iBatch)
  batchSize <- params[["batchSize"]]
  strt <- ((iBatch-1)*batchSize+1)
  stp <- min((iBatch*batchSize), length(multipeptide))
  ##### Get XICs for the batch across all runs #####
  XICs <- lapply(strt:stp, function(rownum){
    ##### Get transition_group_id for that peptideID #####
    idx <- which(precursors$peptide_id == peptideIDs[rownum])
    analytes <- precursors[idx, "transition_group_id"]
    xics <- lapply(names(mzPntrs), function(run){
      chromIndices <- prec2chromIndex[[run]][["chromatogramIndex"]][idx]
      if(any(is.na(unlist(chromIndices))) | is.null(unlist(chromIndices))) return(NULL)
      temp <- lapply(chromIndices, function(i1) extractXIC_group(mzPntrs[[run]], i1))
      names(temp) <- as.character(analytes)
      temp
    })
    names(xics) <- names(mzPntrs)
    xics
  })

  mp <- applyFun(strt:stp, function(rownum){
    peptide <- peptideIDs[rownum]
    df <- multipeptide[[rownum]]
    ref <- refRuns[["run"]][rownum]

    idx <- (rownum - (iBatch-1)*batchSize)
    XICs.ref <- XICs[[idx]][[ref]]
    if(is.null(XICs.ref)){
      message("Chromatogram indices for peptide ", peptide, " are missing in ", fileInfo[ref, "runName"])
      message("Skipping peptide ", peptide, " across all runs.")
      return(df)
    } else {
      XICs.ref.s <- lapply(XICs.ref, smoothXICs, type = params[["XICfilter"]], kernelLen = params[["kernelLen"]],
                           polyOrd = params[["polyOrd"]])
      names(XICs.ref.s) <- names(XICs.ref)
    }
    if(params[["smoothPeakArea"]]) XICs.ref <- XICs.ref.s

    ##### Set alignment rank for all precrusors of the peptide in the reference run #####
    analytes <- as.integer(names(XICs.ref))
    refIdx <- which(df[["run"]] == ref & df[["peak_group_rank"]] == 1)
    refIdx <- refIdx[which.min(df$m_score[refIdx])]
    if(length(refIdx)==0) {
      message("Features for peptide ", peptide, " is missing in ", fileInfo[ref, "runName"])
      message("Skipping peptide ", peptide, " across all runs.")
      return(df)
    }
    df[["alignment_rank"]][refIdx] <- 1L
    df.ref <- df[df$run == ref,]
    df.ref <- setOtherPrecursors(df.ref, XICs.ref, analytes, params)
    # Update multipeptide reference intensity if recal is true
    if(params[["recalIntensity"]]) df.ref <- reIntensity(df.ref, XICs.ref, params)

    ##### Align all runs to reference run and set their alignment rank #####
    exps <- setdiff(rownames(fileInfo), ref)
    df.exps <- lapply(exps, alignToRef2, ref, idx, analytes, fileInfo, XICs, XICs.ref.s, params,
                      df, globalFits, RSE)

    ##### Return the dataframe with alignment rank set to TRUE #####
    updateOnalignTargetedRuns(rownum)
    newDF <- dplyr::bind_rows(df.ref, df.exps)
    newDF
  })
  mp
}


alignToRef2 <- function(eXp, ref, idx, analytes, fileInfo, XICs, XICs.ref.s, params,
                       df, globalFits, RSE){
  # Get XIC_group from experiment run. if missing, go to next run.
  XICs.eXp <- XICs[[idx]][[eXp]]
  df.eXp <- df[df[["run"]] == eXp, ]
  if(is.null(XICs.eXp)){
    message("Chromatogram indices for precursor ", analytes, " are missing in ", fileInfo[eXp, "runName"])
    message("Skipping precursor ", analytes, " in ", fileInfo[eXp, "runName"], ".")
    return(df.eXp)
  } else {
    XICs.eXp.s <- lapply(XICs.eXp, smoothXICs, type = params[["XICfilter"]], kernelLen = params[["kernelLen"]],
                         polyOrd = params[["polyOrd"]])
    names(XICs.eXp.s) <- names(XICs.eXp)
  }
  if(params[["smoothPeakArea"]]) XICs.eXp <- XICs.eXp.s

  # Select 1) all precursors OR 2) high quality precursor
  if(FALSE){
    # Turned off as precursor XICs have different time ranges.
    XICs.ref.pep <- unlist(XICs.ref.s, recursive = FALSE, use.names = FALSE)
    XICs.eXp.pep <- unlist(XICs.eXp.s, recursive = FALSE, use.names = FALSE)
  } else {
    temp <- df[df$run == ref, c("transition_group_id", "m_score")]
    analyte_chr <- as.character(temp[which.min(temp$m_score), "transition_group_id"])
    XICs.ref.pep <- XICs.ref.s[[analyte_chr]]
    XICs.eXp.pep <- XICs.eXp.s[[analyte_chr]]
  }

  ##### Get the aligned Indices #####
  pair <- paste(ref, eXp, sep = "_")
  globalFit <- globalFits[[pair]]
  adaptiveRT <- params[["RSEdistFactor"]]*RSE[[pair]]

  tAligned <- tryCatch(expr = getAlignedTimes(XICs.ref.pep, XICs.eXp.pep, globalFit, params[["alignType"]], adaptiveRT,
           params[["normalization"]], params[["simMeasure"]], params[["goFactor"]],
           params[["geFactor"]], params[["cosAngleThresh"]], params[["OverlapAlignment"]],
           params[["dotProdThresh"]], params[["gapQuantile"]], params[["kerLen"]],
           params[["hardConstrain"]], params[["samples4gradient"]], objType = "light"),
             error = function(e){
             message("\nError in the alignment of ", paste0(analytes, sep = " "), "in runs ",
                     fileInfo[ref, "runName"], " and ", fileInfo[eXp, "runName"])
             warning(e)
             return(df.eXp)
           })
  df.eXp <- tryCatch(expr = setAlignmentRank(df, ref, eXp, tAligned, XICs.eXp, params, adaptiveRT),
             error = function(e){
             message("\nError in setting alignment rank of ", paste0(analytes, sep = " "), "in runs ",
                     fileInfo[eXp, "runName"], " and ", fileInfo[eXp, "runName"])
             warning(e)
             return(df.eXp)
           })
  df.eXp <- setOtherPrecursors(df.eXp, XICs.eXp, analytes, params)
  if(params[["recalIntensity"]]) df.eXp <- reIntensity(df.eXp, XICs.eXp, params)
  df.eXp
}
