#' Progressive alignment for few ids
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2021) + GPL-3
#' Date: 2021-01-20
#' @inheritParams progAlignRuns
#' @return (None)
#' @examples
#' \dontrun{
#' ids <- as.integer(scan(file = "data/ids.txt"))
#' params <- paramsDIAlignR()
#' params[["maxFdrQuery"]] <- 1.0
#' params[["maxPeptideFdr"]] <- 1.0
#' params[["kernelLen"]] <- 11
#' params$batchSize <- 10000L
#' params$globalAlignmentFdr <- 0.001
#' BiocParallel::register(BiocParallel::MulticoreParam(workers = 6, log = FALSE, threshold = "INFO", stop.on.error = TRUE))
#' ropenms <- get_ropenms(condaEnv = "TricEnvr")
#' progAlignRuns2(dataPath = ".", params = params, outFile = "prog.tsv", ropenms = ropenms, ids = ids, applyFun = BiocParallel::bplapply)
#' }
#' @seealso \code{\link{progAlignRuns}}
#' @keywords internal
progAlignRuns2 <- function(dataPath, params, outFile = "DIAlignR.tsv", ropenms, ids= NULL, oswMerged = TRUE,
                          runs = NULL, newickTree = NULL, applyFun = lapply){
  #### Get filenames from .osw file and check consistency between osw and mzML files. #################
  fileInfo <- getRunNames(dataPath, oswMerged)
  fileInfo <- updateFileInfo(fileInfo, runs)
  runs <- rownames(fileInfo)
  message("Following runs will be aligned:")
  print(fileInfo[, "runName", drop=FALSE], sep = "\n")

  #### Get Precursors from the query and respectve chromatogram indices. ######
  # Get all the precursor IDs, transition IDs, Peptide IDs, Peptide Sequence Modified, Charge.
  precursors <- getPrecursors(fileInfo, oswMerged, runType = params[["runType"]],
                              context = params[["context"]], maxPeptideFdr = params[["maxPeptideFdr"]])
  precursors <- dplyr::arrange(precursors, .data$peptide_id, .data$transition_group_id)

  #### Get OpenSWATH peak-groups and their retention times. ##########
  features <- getFeatures(fileInfo, params[["maxFdrQuery"]], params[["runType"]], applyFun)

  #### Precursors for which ids are manually annotated. ##############
  tmp <- applyFun(features, function(df)
    df[df[["m_score"]] <= params[["analyteFDR"]] & df[["peak_group_rank"]] == 1, "transition_group_id"])
  allIDs <- unique(unlist(tmp, recursive = FALSE, use.names = FALSE))
  allIDs <- precursors[precursors[["transition_group_id"]] %in% allIDs, "peptide_id"] # Get respective peptide_id
  set.seed(1)
  allIDs <- sample(unique(allIDs), min(700, length(allIDs)), replace = FALSE) # Selecting 700 ids below analyte FDR fo good global fit.
  allIDs <- unique(c(allIDs, ids)) # ids that we are interested in are also added.
  precursors <- precursors[precursors[["peptide_id"]] %in% allIDs, ]
  rownames(precursors) <- NULL
  message(nrow(precursors), " precursors have features identified in osw files.")

  #### Get Peptide scores, pvalue and qvalues. ######
  peptideIDs <- unique(precursors$peptide_id)
  message(paste0(setdiff(ids, peptideIDs), sep = " "), "peptides are missing.")
  peptideScores <- getPeptideScores(fileInfo, peptideIDs, oswMerged, params[["runType"]], params[["context"]])
  peptideScores <- applyFun(peptideIDs, function(pep) dplyr::filter(peptideScores, .data$peptide_id == pep))
  names(peptideScores) <- as.character(peptideIDs)
  peptideScores <- list2env(peptideScores, hash = TRUE)

  ##### Get distances among runs based on the number of high-quality features. #####
  tmp <- applyFun(features, function(df)
    df[df[["m_score"]] <= params[["analyteFDR"]] & df[["peak_group_rank"]] == 1, "transition_group_id"])
  tmp <- tmp[order(names(tmp), decreasing = FALSE)]
  allIDs <- unique(unlist(tmp, recursive = FALSE, use.names = FALSE))
  allIDs <- sort(allIDs)
  distMat <- length(allIDs) - crossprod(table(stack(tmp)))
  distMat <- stats::dist(distMat, method = "manhattan")

  #### Get the guidance tree. ####
  #start_time <- Sys.time()
  if(is.null(newickTree)){
    tree <- getTree(distMat)
    # Check validity of tree: Names are run and master only.
  } else{
    tree <- ape::read.tree(text = newickTree)
  }

  #### Collect pointers for each mzML file. #######
  message("Collecting metadata from mzML files.")
  mzPntrs <- list2env(getMZMLpointers(fileInfo), hash = TRUE)
  message("Metadata is collected from mzML files.")

  #### Get chromatogram Indices of precursors across all runs. ############
  message("Collecting chromatogram indices for all precursors.")
  prec2chromIndex <- list2env(getChromatogramIndices(fileInfo, precursors, mzPntrs, applyFun), hash = TRUE)

  #### Convert features into multi-peptide #####
  start_time <- Sys.time()
  message("Building multipeptide.")
  multipeptide <- getMultipeptide(precursors, features, applyFun)
  end_time <- Sys.time()
  message("The execution time for building multipeptide:")
  print(end_time - start_time)
  message(length(multipeptide), " peptides are in the multipeptide.")

  #### Get all the child runs through hybrid alignment. ####
  adaptiveRTs <- new.env(hash = TRUE)
  refRuns <- new.env(hash = TRUE)
  features <- list2env(features, hash = TRUE)

  # Traverse up the tree
  start_time <- Sys.time()
  multipeptide <- traverseUp(tree, dataPath, fileInfo, features, mzPntrs, prec2chromIndex, precursors,
             params, adaptiveRTs, refRuns, multipeptide, peptideScores, ropenms, applyFun)
  end_time <- Sys.time() # Report the execution time for hybrid alignment step.
  message("The execution time for creating a master run by alignment:")
  print(end_time - start_time)

  #### Map Ids from the master1 run to all parents. ####
  start_time <- Sys.time()
  multipeptide <- traverseDown(tree, dataPath, fileInfo, multipeptide, prec2chromIndex, mzPntrs, precursors,
               adaptiveRTs, refRuns, params, applyFun)
  end_time <- Sys.time()
  message("The execution time for transfering peaks from root to runs:")
  print(end_time - start_time)

  #### Save features and add master runs to osw #####
  addMasterToOSW(dataPath, tree$node.label, oswMerged)
  filename <- file.path(dataPath, "multipeptide.rds")
  saveRDS(multipeptide, file = filename)

  #### Cleanup.  #######
  rm(mzPntrs)

  #### Write tables to the disk  #######
  finalTbl <- writeTables(fileInfo, multipeptide, precursors)
  utils::write.table(finalTbl, file = outFile, sep = "\t", row.names = FALSE)
  message("Retention time alignment across runs is done.")
  message(paste0(outFile, " file has been written."))

  #### End of function. #####
  alignmentStats(finalTbl, params)
  message("DONE DONE.")
}


# ids <- as.integer(scan(file = "data/ids.txt"))
# params <- paramsDIAlignR()
# params[["maxFdrQuery"]] <- 1.0
# params[["maxPeptideFdr"]] <- 1.0
# params[["kernelLen"]] <- 11
# params$globalAlignmentFdr <- 0.001
# params[["batchSize"]] <- 1000L
# BiocParallel::register(BiocParallel::MulticoreParam(workers = 6, log = FALSE, threshold = "INFO", stop.on.error = TRUE))
# DIAlignR:::alignTargetedRuns2(dataPath = ".", params = params, outFile = "te1.tsv", ids = ids, applyFun = BiocParallel::bplapply)
alignTargetedRuns2 <- function(dataPath, outFile = "DIAlignR.tsv", params, ids= NULL, oswMerged = TRUE, runs = NULL,
                              applyFun = lapply){
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
  precursors <- getPrecursors(fileInfo, oswMerged, params[["runType"]], params[["context"]], params[["maxPeptideFdr"]])
  precursors <- dplyr::arrange(precursors, .data$peptide_id, .data$transition_group_id)

  #### Get OpenSWATH peak-groups and their retention times. ##########
  features <- getFeatures(fileInfo, params[["maxFdrQuery"]], params[["runType"]], applyFun)

  #### Precursors for which ids are manually annotated. ##############
  tmp <- lapply(features, function(df)
    df[df[["m_score"]] <= params[["analyteFDR"]] & df[["peak_group_rank"]] == 1, "transition_group_id"])
  allIDs <- unique(unlist(tmp, recursive = FALSE, use.names = FALSE))
  allIDs <- precursors[precursors[["transition_group_id"]] %in% allIDs, "peptide_id"] # Get respective peptide_id
  set.seed(1)
  allIDs <- sample(unique(allIDs), min(700, length(allIDs)), replace = FALSE) # Selecting 700 ids below analyte FDR fo good global fit.
  allIDs <- unique(c(allIDs, ids)) # ids that we are interested in are also added.
  precursors <- precursors[precursors[["peptide_id"]] %in% allIDs, ]
  message(nrow(precursors), " precursors have features identified in osw files.")

  #### Get Peptide scores, pvalue and qvalues. ######
  peptideIDs <- unique(precursors$peptide_id)
  message(paste0(setdiff(ids, peptideIDs), sep = " "), "are missing.")
  peptideScores <- getPeptideScores(fileInfo, peptideIDs, oswMerged, params[["runType"]], params[["context"]])
  peptideScores <- lapply(peptideIDs, function(pep) dplyr::filter(peptideScores, .data$peptide_id == pep))
  names(peptideScores) <- as.character(peptideIDs)

  #### Collect pointers for each mzML file. #######
  message("Collecting metadata from mzML files.")
  mzPntrs <- getMZMLpointers(fileInfo)
  message("Metadata is collected from mzML files.")

  #### Get chromatogram Indices of precursors across all runs. ############
  message("Collecting chromatogram indices for all precursors.")
  prec2chromIndex <- getChromatogramIndices(fileInfo, precursors, mzPntrs, applyFun)

  #### Convert features into multi-peptide #####
  message("Building multipeptide.")
  multipeptide <- getMultipeptide(precursors, features, applyFun)
  message(length(multipeptide), " peptides are in the multipeptide.")

  #### Get reference run for each precursor ########
  message("Calculating reference run for each peptide.")
  refRuns <- getRefRun(peptideScores, applyFun)

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
  utils::write.table(finalTbl, file = outFile, sep = "\t", row.names = FALSE, quote = FALSE)
  message("Retention time alignment across runs is done.")
  message(paste0(outFile, " file has been written."))

  #### Write alignment summary  #######
  alignmentStats(finalTbl, params)
}

