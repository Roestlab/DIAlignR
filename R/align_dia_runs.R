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
#' @param dataPath (string) path to mzml and osw directory.
#' @param outFile (string) name of the output file.
#' @param oswMerged (logical) TRUE for experiment-wide FDR and FALSE for run-specific FDR by pyprophet.
#' @param runs (A vector of string) names of mzml file without extension.
#' @param runType (string) must be one of the strings "DIA_proteomics", "DIA_Metabolomics".
#' @param context (string) Context used in pyprophet peptide. Must be either "run-specific", "experiment-wide", or "global".
#' @param maxPeptideFdr (numeric) A numeric value between 0 and 1. It is used to filter peptides from osw file which have SCORE_PEPTIDE.QVALUE less than itself.
#' @param maxFdrQuery (numeric) a numeric value between 0 and 1. It is used to filter features from osw file which have SCORE_MS2.QVALUE less than itself.
#' @param XICfilter (string) must be either sgolay, boxcar, gaussian, loess or none.
#' @param polyOrd (integer) order of the polynomial to be fit in the kernel.
#' @param kernelLen (integer) number of data-points to consider in the kernel.
#' @param globalAlignment (string) must be from "loess" or "linear".
#' @param globalAlignmentFdr (numeric) a numeric value between 0 and 1. Features should have m-score lower than this value for participation in LOESS fit.
#' @param globalAlignmentSpan (numeric) spanvalue for LOESS fit. For targeted proteomics 0.1 could be used.
#' @param RSEdistFactor (numeric) defines how much distance in the unit of rse remains a noBeef zone.
#' @param normalization (character) must be selected from "mean", "l2".
#' @param simMeasure (string) must be selected from dotProduct, cosineAngle, crossCorrelation,
#'   cosine2Angle, dotProductMasked, euclideanDist, covariance and correlation.
#' @param alignType available alignment methods are "global", "local" and "hybrid".
#' @param goFactor (numeric) penalty for introducing first gap in alignment. This value is multiplied by base gap-penalty.
#' @param geFactor (numeric) penalty for introducing subsequent gaps in alignment. This value is multiplied by base gap-penalty.
#' @param cosAngleThresh (numeric) in simType = dotProductMasked mode, angular similarity should be higher than cosAngleThresh otherwise similarity is forced to zero.
#' @param OverlapAlignment (logical) an input for alignment with free end-gaps. False: Global alignment, True: overlap alignment.
#' @param dotProdThresh (numeric) in simType = dotProductMasked mode, values in similarity matrix higher than dotProdThresh quantile are checked for angular similarity.
#' @param gapQuantile (numeric) must be between 0 and 1. This is used to calculate base gap-penalty from similarity distribution.
#' @param kerLen (integer) In simType = crossCorrelation, length of the kernel used to sum similarity score. Must be an odd number.
#' @param hardConstrain (logical) if FALSE; indices farther from noBeef distance are filled with distance from linear fit line.
#' @param samples4gradient (numeric) modulates penalization of masked indices.
#' @param analyteFDR (numeric) defines the upper limit of FDR on a precursor to be considered for multipeptide.
#' @param unalignedFDR (numeric) must be between 0 and maxFdrQuery. Features below unalignedFDR are
#'  considered for quantification even without the RT alignment.
#' @param alignedFDR (numeric) must be between unalignedFDR and 1. Features below alignedFDR are
#'  considered for quantification after the alignment.
#' @param baselineType (string) method to estimate the background of a peak contained in XICs. Must be
#'  from "base_to_base", "vertical_division_min", "vertical_division_max".
#' @param integrationType (string) method to ompute the area of a peak contained in XICs. Must be
#'  from "intensity_sum", "trapezoid", "simpson".
#' @param fitEMG (logical) enable/disable exponentially modified gaussian peak model fitting.
#' @param recalIntensity (logical) recalculate intensity for all analytes.
#' @param fillMissing (logical) calculate intensity for ananlytes for which features are not found.
#' @param smoothPeakArea (logical) FALSE: raw chromatograms will be used for quantification. TRUE: smoothed chromatograms will be used for quantification.
#' @param applyFun (function) value must be either lapply or BiocParallel::bplapply.
#' @return An output table with following columns: precursor, run, intensity, RT, leftWidth, rightWidth,
#'  peak_group_rank, m_score, alignment_rank, peptide_id, sequence, charge, group_label.
#'
#' @seealso \code{\link{getRunNames}, \link{getFeatures}, \link{setAlignmentRank}, \link{getMultipeptide}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' alignTargetedRuns(dataPath, outFile = "testDIAlignR.tsv", oswMerged = TRUE,
#'  context = "experiment-wide", maxPeptideFdr = 0.05)
#' @references Gupta S, Ahadi S, Zhou W, Röst H. "DIAlignR Provides Precise Retention Time Alignment Across Distant Runs in DIA and Targeted Proteomics." Mol Cell Proteomics. 2019 Apr;18(4):806-817. doi: https://doi.org/10.1074/mcp.TIR118.001132 Epub 2019 Jan 31.
#'
#' @export
alignTargetedRuns <- function(dataPath, outFile = "DIAlignR.tsv", oswMerged = TRUE, runs = NULL,
                              runType = "DIA_Proteomics", context = "global", maxPeptideFdr = 0.05,
                              maxFdrQuery = 0.05, XICfilter = "sgolay", polyOrd = 4, kernelLen = 9,
                              globalAlignment = "loess", globalAlignmentFdr = 0.01, globalAlignmentSpan = 0.1,
                              RSEdistFactor = 3.5, normalization = "mean", simMeasure = "dotProductMasked",
                              alignType = "hybrid", goFactor = 0.125, geFactor = 40,
                              cosAngleThresh = 0.3, OverlapAlignment = TRUE,
                              dotProdThresh = 0.96, gapQuantile = 0.5, kerLen = 9,
                              hardConstrain = FALSE, samples4gradient = 100,
                              analyteFDR = 0.01,
                              unalignedFDR = 0.01, alignedFDR = 0.05,
                              baselineType = "base_to_base", integrationType = "intensity_sum",
                              fitEMG = FALSE, recalIntensity = FALSE, fillMissing = TRUE, smoothPeakArea = FALSE,
                              applyFun = lapply){
  #### Check if filter length is odd for Savitzky-Golay filter.  #########
  params <- list(runType = runType, context = context, maxPeptideFdr = maxPeptideFdr,
                 maxFdrQuery = maxFdrQuery, XICfilter = XICfilter, polyOrd = polyOrd, kernelLen = kernelLen,
                 globalAlignment = globalAlignment, globalAlignmentFdr = globalAlignmentFdr, globalAlignmentSpan = globalAlignmentSpan,
                 RSEdistFactor = RSEdistFactor, normalization = normalization, simMeasure = simMeasure,
                 alignType = alignType, goFactor = goFactor, geFactor = geFactor,
                 cosAngleThresh = cosAngleThresh, OverlapAlignment = OverlapAlignment,
                 dotProdThresh = dotProdThresh, gapQuantile = gapQuantile, kerLen = kerLen,
                 hardConstrain = hardConstrain, samples4gradient = samples4gradient,
                 analyteFDR = analyteFDR, unalignedFDR = unalignedFDR, alignedFDR = alignedFDR,
                 baselineType = baselineType, integrationType = integrationType,
                 fitEMG = fitEMG, recalIntensity = recalIntensity, fillMissing = fillMissing,
                 smoothPeakArea = smoothPeakArea)
  checkParams(params)

  #### Get filenames from .osw file and check consistency between osw and mzML files. #################
  fileInfo <- getRunNames(dataPath, oswMerged)
  fileInfo <- updateFileInfo(fileInfo, runs)
  runs <- rownames(fileInfo)
  message("Following runs will be aligned:")
  print(fileInfo[, "runName"], sep = "\n")

  #### Get Precursors from the query and respectve chromatogram indices. ######
  # Get all the precursor IDs, transition IDs, Peptide IDs, Peptide Sequence Modified, Charge.
  precursors <- getPrecursors(fileInfo, oswMerged, runType, context, maxPeptideFdr)

  #### Get OpenSWATH peak-groups and their retention times. ##########
  features <- getFeatures(fileInfo, maxFdrQuery, runType)

  #### Precursors for which features are identified. ##############
  allIDs <- unique(unlist(lapply(features, function(df) df[df[["m_score"]] <= analyteFDR,
                                                           "transition_group_id"]),
                          recursive = FALSE, use.names = FALSE))
  precursors <- precursors[precursors[["transition_group_id"]] %in% allIDs, ]

  #### Collect pointers for each mzML file. #######
  message("Collecting metadata from mzML files.")
  mzPntrs <- getMZMLpointers(fileInfo)
  message("Metadata is collected from mzML files.")

  #### Get chromatogram Indices of precursors across all runs. ############
  message("Collecting chromatogram indices for all precursors.")
  prec2chromIndex <- getChromatogramIndices(fileInfo, precursors, mzPntrs)

  #### Convert features into multi-precursor #####
  message("Building multipeptide.")
  multipeptide <- getMultipeptide(precursors, features)
  message(length(multipeptide), " precursors are in the multipeptide")

  #### Get reference run for each precursor ########
  message("Calculating reference run for each precursor.")
  refRuns <- getRefRun(multipeptide)

  #### Container to save Global alignments.  #######
  globalFits <- getGlobalFits(refRuns, features, fileInfo, globalAlignment,
                              globalAlignmentFdr, globalAlignmentSpan)
  RSE <- lapply(globalFits, getRSE)

  #### Perform pairwise alignment ###########
  message("Performing reference-based alignment.")
  num_of_prec <- length(multipeptide)
  start_time <- Sys.time()
  multipeptide <- lapply(seq_along(multipeptide), alignIthAnalyte, precursors, refRuns, fileInfo, prec2chromIndex,
               multipeptide, params, mzPntrs, globalFits, RSE, applyFun)
  names(multipeptide) <- as.character(precursors[["transition_group_id"]])

  #### Cleanup.  #######
  rm(mzPntrs)
  end_time <- Sys.time() # Report the execution time for hybrid alignment step.
  message("The execution time for alignment:")
  print(end_time - start_time)

  #### Write tables to the disk  #######
  finalTbl <- writeTables(fileInfo, multipeptide, precursors)
  utils::write.table(finalTbl, file = outFile, sep = "\t", row.names = FALSE)
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
#' @param analytes (vector of integers) transition_group_ids for which features are to be extracted.
#' @param runs (A vector of string) Names of mzml file without extension.
#' @param dataPath (char) Path to mzml and osw directory.
#' @param refRun (string) reference for alignment. If no run is provided, m-score is used to select reference run.
#' @param oswMerged (logical) TRUE for experiment-wide FDR and FALSE for run-specific FDR by pyprophet.
#' @param runType (char) This must be one of the strings "DIA_proteomics", "DIA_Metabolomics".
#' @param maxFdrQuery (numeric) A numeric value between 0 and 1. It is used to filter features from osw file which have SCORE_MS2.QVALUE less than itself.
#' @param analyteFDR (numeric) only analytes that have m-score less than this, will be included in the output.
#' @param XICfilter (string) must be either sgolay, boxcar, gaussian, loess or none.
#' @param polyOrd (integer) order of the polynomial to be fit in the kernel.
#' @param kernelLen (integer) number of data-points to consider in the kernel.
#' @param globalAlignment (string) must be from "loess" or "linear".
#' @param globalAlignmentFdr (numeric) a numeric value between 0 and 1. Features should have m-score lower than this value for participation in LOESS fit.
#' @param globalAlignmentSpan (numeric) spanvalue for LOESS fit. For targeted proteomics 0.1 could be used.
#' @param RSEdistFactor (numeric) defines how much distance in the unit of rse remains a noBeef zone.
#' @param normalization (character) must be selected from "mean", "l2".
#' @param simMeasure (string) must be selected from dotProduct, cosineAngle, crossCorrelation,
#'   cosine2Angle, dotProductMasked, euclideanDist, covariance and correlation.
#' @param alignType available alignment methods are "global", "local" and "hybrid".
#' @param goFactor (numeric) penalty for introducing first gap in alignment. This value is multiplied by base gap-penalty.
#' @param geFactor (numeric) penalty for introducing subsequent gaps in alignment. This value is multiplied by base gap-penalty.
#' @param cosAngleThresh (numeric) in simType = dotProductMasked mode, angular similarity should be higher than cosAngleThresh otherwise similarity is forced to zero.
#' @param OverlapAlignment (logical) an input for alignment with free end-gaps. False: Global alignment, True: overlap alignment.
#' @param dotProdThresh (numeric) in simType = dotProductMasked mode, values in similarity matrix higher than dotProdThresh quantile are checked for angular similarity.
#' @param gapQuantile (numeric) must be between 0 and 1. This is used to calculate base gap-penalty from similarity distribution.
#' @param kerLen (integer) In simType = crossCorrelation, length of the kernel used to sum similarity score. Must be an odd number.
#' @param hardConstrain (logical) if FALSE; indices farther from noBeef distance are filled with distance from linear fit line.
#' @param samples4gradient (numeric) modulates penalization of masked indices.
#' @param objType (char) Must be selected from light, medium and heavy.
#' @return A list of fileInfo and AlignObjs. Each AlignObj is an S4 object. Three most-important slots are:
#' \item{indexA_aligned}{(integer) aligned indices of reference run.}
#' \item{indexB_aligned}{(integer) aligned indices of experiment run.}
#' \item{score}{(numeric) cumulative score of alignment.}
#' @seealso \code{\link{plotAlignedAnalytes}, \link{getRunNames}, \link{getFeatures}, \link{getXICs4AlignObj}, \link{getAlignObj}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
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
                         runType = "DIA_Proteomics", maxFdrQuery = 0.05, analyteFDR = 0.01,
                         XICfilter = "sgolay", polyOrd = 4, kernelLen = 9,
                         globalAlignment = "loess", globalAlignmentFdr = 0.01, globalAlignmentSpan = 0.1,
                         RSEdistFactor = 3.5, normalization = "mean", simMeasure = "dotProductMasked",
                         alignType = "hybrid", goFactor = 0.125, geFactor = 40,
                         cosAngleThresh = 0.3, OverlapAlignment = TRUE,
                         dotProdThresh = 0.96, gapQuantile = 0.5, kerLen = 9,
                         hardConstrain = FALSE, samples4gradient = 100,
                         objType = "light"){
  if( (kernelLen %% 2) != 1){
    print("kernelLen can only be odd number")
    return(NULL)
  }
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
  features <- getFeatures(filenames, maxFdrQuery, runType)

  ###### Report analytes that are not found ########
  refAnalytes <- analytesFromFeatures(features, analyteFDR = analyteFDR, commonAnalytes = FALSE)
  analytesFound <- intersect(analytes, refAnalytes)
  analytesNotFound <- setdiff(analytes, analytesFound)
  if(length(analytesNotFound)>0){
    message(paste(analytesNotFound, "not found with FDR cut-off."))
  }
  analytes <- analytesFound
  precursors <- precursors[precursors[["transition_group_id"]] %in% analytes, ]
  if(nrow(precursors) == 0){
    stop("No precursors are found below ", analyteFDR)
  }

  ############# Get chromatogram Indices of precursors across all runs. ############
  prec2chromIndex <- getChromatogramIndices(filenames, precursors, mzPntrs)

  ############ Convert features into multi-precursor #####
  multipeptide <- getMultipeptide(precursors, features)

  ############## Get reference run for each precursor ########
  idx <- which(filenames$runName == refRun)
    if(length(idx) == 0){
      print("Finding reference run using m-score.")
      refRun <- getRefRun(multipeptide)
    } else{
      run <- rownames(filenames)[idx]
      refRun <- data.frame("transition_group_id" = as.integer(names(multipeptide)),
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
    XICs.ref.s <- smoothXICs(XICs.ref, type = XICfilter, kernelLen = kernelLen, polyOrd = polyOrd)
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
      XICs.eXp.s <- smoothXICs(XICs.eXp, type = XICfilter, kernelLen = kernelLen, polyOrd = polyOrd)
      # Get the loess fit for hybrid alignment
      if(any(pair %in% names(globalFits))){
        globalFit <- globalFits[[pair]]
      } else{
        globalFit <- getGlobalAlignment(features, ref, eXp,
                                        globalAlignment, globalAlignmentFdr, globalAlignmentSpan)
        globalFits[[pair]] <- globalFit
        RSE[[pair]] <- getRSE(globalFit)
      }
      adaptiveRT <- RSEdistFactor*RSE[[pair]]

      # Fetch alignment object between XICs.ref and XICs.eXp
      AlignObj <- getAlignObj(XICs.ref.s, XICs.eXp.s, globalFit, alignType, adaptiveRT,
                              normalization, simType = simMeasure, goFactor, geFactor,
                              cosAngleThresh, OverlapAlignment, dotProdThresh, gapQuantile,
                              kerLen, hardConstrain, samples4gradient,
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
alignToRef <- function(eXp, ref, analyte_chr, rownum, fileInfo, XICs.ref.s, params, prec2chromIndex, mzPntrs, df, globalFits, RSE){
  # Get XIC_group from experiment run. if missing, go to next run.
  chromIndices <- prec2chromIndex[[eXp]][["chromatogramIndex"]][[rownum]]
  if(any(is.na(chromIndices))){
    warning("Chromatogram indices for ", analyte_chr, " are missing in ", fileInfo[eXp, "runName"])
    message("Skipping ", analyte_chr, " in ", fileInfo[eXp, "runName"], ".")
    df <- df[df[["run"]] == eXp, ]
    return(df)
  } else {
    XICs.eXp <- extractXIC_group(mzPntrs[[eXp]], chromIndices)
    XICs.eXp.s <- smoothXICs(XICs.eXp, type = params[["XICfilter"]],
                                     kernelLen = params[["kernelLen"]], polyOrd = params[["polyOrd"]])
  }

  pair <- paste(ref, eXp, sep = "_")
  globalFit <- globalFits[[pair]]
  adaptiveRT <- params[["RSEdistFactor"]]*RSE[[pair]]
  # Get the aligned Indices
  tAligned <- getAlignedTimes(XICs.ref.s, XICs.eXp.s, globalFit, params[["alignType"]], adaptiveRT,
                               params[["normalization"]], params[["simMeasure"]], params[["goFactor"]],
                               params[["geFactor"]], params[["cosAngleThresh"]], params[["OverlapAlignment"]],
                               params[["dotProdThresh"]], params[["gapQuantile"]], params[["kerLen"]],
                               params[["hardConstrain"]], params[["samples4gradient"]], objType = "light")

  if(params[["smoothPeakArea"]]){
    newdf <- setAlignmentRank(df, ref, eXp, tAligned, XICs.eXp.s, params, adaptiveRT)
  } else {
    newdf <- setAlignmentRank(df, ref, eXp, tAligned, XICs.eXp, params, adaptiveRT)
  }
  newdf <- newdf[newdf[["run"]] == eXp, ]
  newdf
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
#' @inheritParams traverseUp
#' @param rownum (integer) represnts the index of the multipepetide to be aligned.
#' @param refRuns (data-frame) output of \code{\link{getRefRun}}. Must have two columsn : transition_group_id and run.
#' @param multipeptide (list) contains multiple data-frames that are collection of features
#'  associated with analytes. This is an output of \code{\link{getMultipeptide}}.
#' @param globalFits (list) each element is either of class lm or loess. This is an output of \code{\link{getGlobalFits}}.
#' @param RSE (list) Each element represents Residual Standard Error of corresponding fit in globalFits.
#' @param applyFun (function) value must be either lapply or BiocParallel::bplapply.
#' @return (dataframe) a collection of aligned features related to analyte_chr.
#' @seealso \code{\link{alignTargetedRuns}, \link{alignToRef}, \link{getAlignedTimes}, \link{getMultipeptide}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
alignIthAnalyte <- function(rownum, precursors, refRuns, fileInfo, prec2chromIndex, multipeptide, params,
                    mzPntrs, globalFits, RSE, applyFun){
  analyte <- precursors[["transition_group_id"]][rownum]
  analyte_chr <- as.character(analyte)
  ref <- refRuns[["run"]][rownum]
  exps <- setdiff(rownames(fileInfo), ref)
  chromIndices <- prec2chromIndex[[ref]][["chromatogramIndex"]][[rownum]]

  # Get XIC_group from reference run. if missing, return unaligned features.
  if(any(is.na(chromIndices)) | is.null(chromIndices)){
    warning("Chromatogram indices for ", analyte, " are missing in ", fileInfo[ref, "runName"])
    message("Skipping ", analyte, " across all runs.")
    return(multipeptide[[analyte_chr]])
  } else {
    XICs.ref <- extractXIC_group(mz = mzPntrs[[ref]], chromIndices = chromIndices)
    XICs.ref.s <- smoothXICs(XICs.ref, type = params[["XICfilter"]], kernelLen = params[["kernelLen"]],
                             polyOrd = params[["polyOrd"]])
  }

  df <- multipeptide[[analyte_chr]]
  refIdx <- which(df[["run"]] == ref & df[["peak_group_rank"]] == 1)
  df[["alignment_rank"]][refIdx] <- 1L

  df.ref <- df[df[["run"]] == ref,]
  refIdx <- which(df.ref[["alignment_rank"]] == 1)
  # Update multipeptide reference intensity if recal is true
  if(recalIntensity){
    if(params[["smoothPeakArea"]]) XICs.ref <- XICs.ref.s
    area <- calculateIntensity(XICs.ref, df.ref[refIdx, "leftWidth"], df.ref[refIdx, "rightWidth"],
                               params[["integrationType"]], params[["baselineType"]], params[["fitEMG"]])
    df.ref$intensity[refIdx] <- area
  }

  # Align all runs to reference run
  df.exps <- applyFun(exps, alignToRef, ref, analyte_chr, rownum, fileInfo, XICs.ref.s, params, prec2chromIndex, mzPntrs,
                    df, globalFits, RSE)
  updateOnalignTargetedRuns(rownum)
  dplyr::bind_rows(df.ref, df.exps)
}
