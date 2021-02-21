#' Extract features and generate pairwise alignments.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2021) + GPL-3
#' Date: 2021-02-20
#' @inheritParams alignTargetedRuns
#' @return NULL
#'
#' @seealso \code{\link{alignTargetedRuns}}
#' @examples
#' params <- paramsDIAlignR()
#' params[["context"]] <- "experiment-wide"
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' BiocParallel::register(BiocParallel::MulticoreParam(workers = 4, progressbar = TRUE))
#' script1(dataPath, outFile = "testDIAlignR", params = params, applyFun = BiocParallel::bplapply)
#' file.remove(file.path(dataPath, "testDIAlignR_script1.RData"))
#' @export
script1 <- function(dataPath, outFile = "DIAlignR", params = paramsDIAlignR(), oswMerged = TRUE,
                    runs = NULL, applyFun=lapply){
  fileInfo <- getRunNames(dataPath, oswMerged, params)
  fileInfo <- updateFileInfo(fileInfo, runs)
  runs <- rownames(fileInfo)
  message("Following runs will be aligned:")
  print(fileInfo[, "runName"], sep = "\n")

  start_time <- Sys.time()
  if(params[["transitionIntensity"]]){
    features <- getTransitions(fileInfo, params[["maxFdrQuery"]], params[["runType"]], applyFun)
  } else{
    features <- getFeatures(fileInfo, params[["maxFdrQuery"]], params[["runType"]], applyFun)
  }
  end_time <- Sys.time()
  message("The execution time for fetching features:")
  print(end_time - start_time)

  message("Calculating global alignments.")
  start_time <- Sys.time()
  refRuns <- data.frame("run" = rownames(fileInfo))
  globalFits <- getGlobalFits(refRuns, features, fileInfo, params[["globalAlignment"]],
                              params[["globalAlignmentFdr"]], params[["globalAlignmentSpan"]], applyFun)
  RSE <- applyFun(globalFits, getRSE, params[["globalAlignment"]])
  globalFits <- applyFun(globalFits, extractFit, params[["globalAlignment"]])
  end_time <- Sys.time()
  message("The execution time for calculating global alignment:")
  print(end_time - start_time)
  save(features, globalFits, RSE, file = file.path(dataPath, paste0(outFile, "script1.RData", sep = "_")))
  print("script1 is done.")
}

#' Performs alignment using script1 output
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2021) + GPL-3
#' Date: 2021-02-20
#' @importFrom data.table data.table setkeyv
#' @inheritParams alignTargetedRuns
#' @return NULL
#'
#' @seealso \code{\link{alignTargetedRuns}}
#' @examples
#' params <- paramsDIAlignR()
#' params[["context"]] <- "experiment-wide"
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' BiocParallel::register(BiocParallel::MulticoreParam(workers = 4, progressbar = TRUE))
#' script1(dataPath, outFile = "testDIAlignR", params = params, applyFun = BiocParallel::bplapply)
#' script2(dataPath, outFile = "testDIAlignR", params = params, applyFun = lapply)
#' file.remove(file.path(dataPath, "testDIAlignR_script1.RData"))
#' @export
script2 <- function(dataPath, outFile = "DIAlignR", params = paramsDIAlignR(), oswMerged = TRUE,
                    runs = NULL, refRun = NULL, applyFun = lapply){
  load(file = file.path(dataPath, paste0(outFile, "script1.RData", sep = "_")))
  #### Check if all parameters make sense.  #########
  params <- checkParams(params)

  #### Get filenames from .osw file and check consistency between osw and mzML files. #################
  fileInfo <- getRunNames(dataPath, oswMerged, params)
  fileInfo <- updateFileInfo(fileInfo, runs)
  runs <- rownames(fileInfo)
  message("Following runs will be aligned:")
  print(fileInfo[, "runName"], sep = "\n")

  #### Get Precursors from the query and respectve chromatogram indices. ######
  # Get all the precursor IDs, transition IDs, Peptide IDs, Peptide Sequence Modified, Charge.
  start_time <- Sys.time()
  precursors <- getPrecursors(fileInfo, oswMerged, params[["runType"]], params[["context"]], params[["maxPeptideFdr"]], params[["level"]])
  if(params[["fractionPercent"]] != 100L){
    idx <- getPrecursorSubset(precursors, params)
    precursors <- precursors[idx[1]:idx[2],]
    setkeyv(precursors, c("peptide_id", "transition_group_id"))
    outFile <- paste(outFile, params[["fraction"]], params[["fractionPercent"]], sep = "_")
  }
  outFile <- paste0(outFile,".tsv")
  end_time <- Sys.time()
  message("The execution time for getting precursors:")
  print(end_time - start_time)

  #### Get Peptide scores, pvalue and qvalues. ######
  # Some peptides may not be found due to using a subset of runs. Appends NA for them.
  # This translates as "Chromatogram indices for peptide ID are missing in NA"
  start_time <- Sys.time()
  peptideIDs <- precursors[, logical(1), keyby = peptide_id]$peptide_id
  peptideScores <- getPeptideScores(fileInfo, peptideIDs, oswMerged, params[["runType"]], params[["context"]])
  peptideScores <- lapply(peptideIDs, function(pep) peptideScores[.(pep)])
  names(peptideScores) <- as.character(peptideIDs)
  end_time <- Sys.time()
  message("The execution time for fetching peptide scores:")
  print(end_time - start_time)

  #### Get reference run for each precursor ########
  start_time <- Sys.time()
  idx <- which(fileInfo$runName == refRun)
  if(length(idx) == 0){
    message("Calculating reference run for each peptide.")
    refRuns <- getRefRun(peptideScores)
  } else{
    run <- rownames(fileInfo)[idx]
    refRuns <- data.table("peptide_id" = peptideIDs, "run" = run, key = "peptide_id")
  }
  end_time <- Sys.time()
  message("The execution time for calculating a reference run:")
  print(end_time - start_time)
  rm(peptideScores)

  #### Collect pointers for each mzML file. #######
  start_time <- Sys.time()
  message("Collecting metadata from mzML files.")
  mzPntrs <- getMZMLpointers(fileInfo)
  message("Metadata is collected from mzML files.")
  end_time <- Sys.time()
  message("The execution time for getting pointers:")
  print(end_time - start_time)

  #### Get chromatogram Indices of precursors across all runs. ############
  message("Collecting chromatogram indices for all precursors.")
  start_time <- Sys.time()
  prec2chromIndex <- getChromatogramIndices(fileInfo, precursors, mzPntrs, applyFun)
  end_time <- Sys.time()
  message("The execution time for getting chromatogram indices:")
  print(end_time - start_time)

  #### Convert features into multi-peptide #####
  message("Building multipeptide.")
  start_time <- Sys.time()
  multipeptide <- getMultipeptide(precursors, features, applyFun)
  message(length(multipeptide), " peptides are in the multipeptide.")
  end_time <- Sys.time()
  message("The execution time for building multipeptide:")
  print(end_time - start_time)
  rm(features)

  # TODO: Check dimensions of multipeptide, PeptideIDs, precursors etc makes sense.
  #### Perform pairwise alignment ###########
  message("Performing reference-based alignment.")
  start_time <- Sys.time()
  num_of_batch <- ceiling(length(multipeptide)/params[["batchSize"]])
  invisible(
    lapply(1:num_of_batch, perBatch, peptideIDs, multipeptide, refRuns, precursors,
           prec2chromIndex, fileInfo, mzPntrs, params, globalFits, RSE, lapply)
  )

  #### Cleanup.  #######
  for(mz in mzPntrs){
    if(is(mz)[1] == "SQLiteConnection") DBI::dbDisconnect(mz)
    if(is(mz)[1] == "mzRpwiz") rm(mz)
  }
  rm(prec2chromIndex, globalFits, refRuns, RSE)

  end_time <- Sys.time() # Report the execution time for hybrid alignment step.
  message("The execution time for alignment:")
  print(end_time - start_time)

  #### Write tables to the disk  #######
  finalTbl <- writeTables(fileInfo, multipeptide, precursors)
  if(params[["transitionIntensity"]]){
    finalTbl[,intensity := sapply(intensity,function(x) paste(round(x, 3), collapse=", "))]
  }
  utils::write.table(finalTbl, file = outFile, sep = "\t", row.names = FALSE, quote = FALSE)
  message("Retention time alignment across runs is done.")
  message(paste0(outFile, " file has been written."))

  #### Write alignment summary  #######
  alignmentStats(finalTbl, params)
}
