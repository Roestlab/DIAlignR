#' Peptide quantification through progressive alignment
#'
#' This function expects osw and mzml directories at dataPath. It first reads osw files and fetches
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
#' @param dataPath (string) path to mzml and osw directory.
#' @param outFile (string) name of the output file.
#' @param ropenms (pyopenms module) get this python module through \code{\link{get_ropenms}}.
#' @param oswMerged (logical) TRUE for experiment-wide FDR and FALSE for run-specific FDR by pyprophet.
#' @param runs (string) names of mzml file without extension.
#' @param newickTree (string) guidance tree in newick format. Look up \code{\link{getTree}}.
#' @return (None)
#' @seealso \code{\link{alignTargetedRuns}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' params <- paramsDIAlignR()
#' \dontrun{
#' ropenms <- get_ropenms(condaEnv = "envName")
#' progAlignRuns(dataPath, params = params, outFile = "test3.tsv", ropenms = ropenms)
#' # Removing aligned vectors
#' file.remove(list.files(dataPath, pattern = "*_av.rds", full.names = TRUE))
#' # Removing temporarily created master chromatograms
#' file.remove(list.files(file.path(dataPath, "mzml"), pattern = "^master[0-9]+\\.chrom\\.mzML$", full.names = TRUE))
#' }
#' @export
progAlignRuns <- function(dataPath, params, outFile = "DIAlignR.tsv", ropenms, oswMerged = TRUE,
                          runs = NULL, newickTree = NULL){
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

  #### Get OpenSWATH peak-groups and their retention times. ##########
  features <- list2env(getFeatures(fileInfo, params[["maxFdrQuery"]], params[["runType"]]), hash = FALSE)

  #### Precursors for which features are identified. ##############
  tmp <- lapply(features, function(df)
    df[df[["m_score"]] <= params[["analyteFDR"]] & df[["peak_group_rank"]] == 1, "transition_group_id"])
  allIDs <- unique(unlist(tmp, recursive = FALSE, use.names = FALSE))
  precursors <- precursors[precursors[["transition_group_id"]] %in% allIDs, ]
  message(nrow(precursors), " precursors have features identified in osw files.")

  # Get distances among runs based on number of common IDs.
  distMat <- length(allIDs) - crossprod(table(stack(tmp)))
  distMat <- stats::dist(distMat, method = "manhattan")

  #### Collect pointers for each mzML file. #######
  message("Collecting metadata from mzML files.")
  mzPntrs <- list2env(getMZMLpointers(fileInfo), hash = FALSE)
  message("Metadata is collected from mzML files.")

  #### Get chromatogram Indices of precursors across all runs. ############
  message("Collecting chromatogram indices for all precursors.")
  prec2chromIndex <- list2env(getChromatogramIndices(fileInfo, precursors, mzPntrs), hash = FALSE)

  #### Convert features into multi-peptide. #####
  message("Building multipeptide.")
  multipeptide <-  list2env(getMultipeptide(precursors, features))
  message(length(multipeptide), " precursors are in the multipeptide")

  #### Get the guidance tree. ####
  start_time <- Sys.time()
  if(is.null(newickTree)){
    tree <- getTree(distMat)
    # Check validity of tree: Names are run and master only.
  } else{
    tree <- ape::read.tree(text = newickTree)
  }

  #### Get all the child runs through hybrid alignment. ####
  adaptiveRTs <- new.env(hash = TRUE)
  refRuns <- new.env(hash = TRUE)
  # Traverse up the tree
  traverseUp(tree, dataPath, fileInfo, features, mzPntrs, prec2chromIndex, precursors,
             params, adaptiveRTs, refRuns, ropenms)
  end_time <- Sys.time() # Report the execution time for hybrid alignment step.
  message("The execution time for alignment:")
  print(end_time - start_time)

  #### Convert features into multi-peptide #####
  message("Building multipeptide.")
  multipeptide <- list2env(getMultipeptide(precursors, features), hash = TRUE)
  message(length(multipeptide), " precursors are in the multipeptide")

  #### Map Ids from the master1 run to all parents. ####
  analytes <- precursors$transition_group_id
  traverseDown(tree, dataPath, fileInfo, multipeptide, prec2chromIndex, mzPntrs, analytes,
               adaptiveRTs, refRuns, params)

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
