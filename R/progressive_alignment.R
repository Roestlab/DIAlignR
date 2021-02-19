#' Peptide quantification through progressive alignment
#'
#' This function expects osw and xics directories at dataPath. It first reads osw files and fetches
#'  chromatogram indices for each analyte. To perform alignment, first a crude guide-tree is built which
#'  can also be provided with newickTree parameter. As we traverse from the leaf-nodes to the root node,
#'  runs are aligned pairwise. The root node is named master1 that has average of all fragment ion chromatograms
#'  and identified peak-groups. These features are propagated back to leaf nodes and finally aligned
#'  features are written in the output file.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-07-10
#' @inheritParams checkParams
#' @inheritParams alignTargetedRuns
#' @param dataPath (string) path to xics and osw directory.
#' @param outFile (string) name of the output file.
#' @param ropenms (pyopenms module) get this python module through \code{\link{get_ropenms}}. Required only for chrom.mzML files.
#' @param oswMerged (logical) TRUE for experiment-wide FDR and FALSE for run-specific FDR by pyprophet.
#' @param runs (string) names of xics file without extension.
#' @param newickTree (string) guidance tree in newick format. Look up \code{\link{getTree}}.
#' @return (None)
#' @seealso \code{\link{alignTargetedRuns}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' params <- paramsDIAlignR()
#' params[["context"]] <- "experiment-wide"
#' \dontrun{
#' ropenms <- get_ropenms(condaEnv = "envName")
#' progAlignRuns(dataPath, params = params, outFile = "test3", ropenms = ropenms)
#' # Removing aligned vectors
#' file.remove(list.files(dataPath, pattern = "*_av.rds", full.names = TRUE))
#' # Removing temporarily created master chromatograms
#' file.remove(list.files(file.path(dataPath, "xics"), pattern = "^master[0-9]+\\.chrom\\.mzML$", full.names = TRUE))
#' file.remove(file.path(dataPath, "test3.temp.RData"))
#' file.remove(file.path(dataPath, "master.merged.osw"))
#' }
#' @export
progAlignRuns <- function(dataPath, params, outFile = "DIAlignR", ropenms = NULL, oswMerged = TRUE,
                          runs = NULL, newickTree = NULL, applyFun = lapply){
  if(params[["chromFile"]] == "mzML"){
    if(is.null(ropenms)) stop("ropenms is required to write chrom.mzML files.")
  }

  #### Check if all parameters make sense.  #########
  params <- checkParams(params)

  #### Get filenames from .osw file and check consistency between osw and mzML files. #################
  fileInfo <- getRunNames(dataPath, oswMerged, params)
  fileInfo <- updateFileInfo(fileInfo, runs)
  runs <- rownames(fileInfo)
  message("Following runs will be aligned:")
  print(fileInfo[, "runName", drop=FALSE], sep = "\n")

  #### Get Precursors from the query and respectve chromatogram indices. ######
  # Get all the precursor IDs, transition IDs, Peptide IDs, Peptide Sequence Modified, Charge.
  precursors <- getPrecursors(fileInfo, oswMerged, runType = params[["runType"]],
                              context = params[["context"]], maxPeptideFdr = params[["maxPeptideFdr"]])
  precursors <- dplyr::arrange(precursors, .data$peptide_id, .data$transition_group_id)

  #### Get Peptide scores, pvalue and qvalues. ######
  peptideIDs <- unique(precursors$peptide_id)
  peptideScores <- getPeptideScores(fileInfo, peptideIDs, oswMerged, params[["runType"]], params[["context"]])
  peptideScores <- applyFun(peptideIDs, function(pep) dplyr::filter(peptideScores, .data$peptide_id == pep))
  names(peptideScores) <- as.character(peptideIDs)
  peptideScores <- list2env(peptideScores, hash = TRUE)

  #### Get OpenSWATH peak-groups and their retention times. ##########
  features <- getFeatures(fileInfo, params[["maxFdrQuery"]], params[["runType"]], applyFun)

  ##### Get distances among runs based on the number of high-quality features. #####
  tmp <- applyFun(features, function(df)
       df[df[["m_score"]] <= params[["analyteFDR"]] & df[["peak_group_rank"]] == 1, "transition_group_id"])
  tmp <- tmp[order(names(tmp), decreasing = FALSE)]
  allIDs <- unique(unlist(tmp, recursive = FALSE, use.names = TRUE))
  allIDs <- sort(allIDs)
  distMat <- length(allIDs) - crossprod(table(utils::stack(tmp)))
  distMat <- stats::dist(distMat, method = "manhattan")

  #### Get the guidance tree. ####
  start_time <- Sys.time()
  if(is.null(newickTree)){
    tree <- getTree(distMat) # Check validity of tree: Names are run and master only.
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
  multipeptide <- traverseDown(tree, dataPath, fileInfo, multipeptide, prec2chromIndex, mzPntrs,
                                precursors, adaptiveRTs, refRuns, params, applyFun)
  end_time <- Sys.time()
  message("The execution time for transfering peaks from root to runs:")
  print(end_time - start_time)

  #### Save features and add master runs to osw #####
  addMasterToOSW(dataPath, tree$node.label, oswMerged)
  save(multipeptide, fileInfo, peptideScores, refRuns, adaptiveRTs,
       file = file.path(dataPath, paste0(outFile,".temp.RData", sep = "")))

  #### Cleanup.  #######
  for(mz in names(mzPntrs)){
    if(is(mzPntrs[[mz]])[1] == "SQLiteConnection") DBI::dbDisconnect(mzPntrs[[mz]])
  }
  rm(mzPntrs, prec2chromIndex, peptideScores, features)

  #### Write tables to the disk  #######
  finalTbl <- writeTables(fileInfo, multipeptide, precursors)
  filename <- paste0(outFile, ".tsv", sep = "")
  utils::write.table(finalTbl, file = filename, sep = "\t", row.names = FALSE)
  message("Retention time alignment across runs is done.")
  message(paste0(filename, " file has been written."))

  #### End of function. #####
  alignmentStats(finalTbl, params)
  message("DONE DONE.")
}
