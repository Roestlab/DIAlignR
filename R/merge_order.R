#' Number of descendants
#'
#' Get the number of descendants for each node in the tree.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#' @importFrom ape Ntip
#' @param tree (phylo) a phylogenetic tree.
#' @return (numeric)
#' @seealso \code{\link[ape]{Ntip}, \link{getNodeIDs}}
#' @keywords internal
#' @examples
#' m <- matrix(c(0,1,2,3, 1,0,1.5,1.5, 2,1.5,0,1, 3,1.5,1,0), byrow = TRUE,
#'             ncol = 4, dimnames = list(c("run1", "run2", "run3", "run4"),
#'                                       c("run1", "run2", "run3", "run4")))
#' distMat <- as.dist(m, diag = FALSE, upper = FALSE)
#' \dontrun{
#' tree <- getTree(distMat)
#' getNodeIDs(tree)
#' nrDesc(tree)
#' }
nrDesc <- function(tree) {
  res <- numeric(max(tree$edge))
  res[1:Ntip(tree)] <- 1L
  for (i in postorder(tree)) {
    tmp <- tree$edge[i,1]
    res[tmp] <- res[tmp] + res[tree$edge[i, 2]]
  }
  res
}


# nr_desc(k)
# phangorn::Descendants(k, node = "6", type = "children")
# getMRCA(k, c("1", "4"))

#' Create a phylogenetic tree
#'
#' Builds a phylogenetic tree from the distance matrix using UPGMA algorithm.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-05-31
#' @import phangorn ape
#' @param distMat (dist) a pairwise distance matrix.
#' @return (phylo) a tree.
#' @seealso \code{\link[phangorn]{upgma}, \link{getNodeIDs}}
#' @keywords internal
#' @examples
#' m <- matrix(c(0,1,2,3, 1,0,1.5,1.5, 2,1.5,0,1, 3,1.5,1,0), byrow = TRUE,
#'             ncol = 4, dimnames = list(c("run1", "run2", "run3", "run4"),
#'                                       c("run1", "run2", "run3", "run4")))
#' distMat <- as.dist(m, diag = FALSE, upper = FALSE)
#' labels(distMat); length(distMat)
#' \dontrun{
#' getTree(distMat)
#' }
#' tree <- ape::read.tree(text = "(run1:9,(run2:7,run0:2)master2:5)master1;")
#' plot(tree, show.node.label = TRUE)
getTree <- function(distMat){
  tree <- phangorn::upgma(distMat)
  tree <- ape::reorder.phylo(tree, "postorder")
  tree <- ape::makeNodeLabel(tree, method = "number", prefix = "master")
  message("alignment order of runs in Newick format:")
  message(ape::write.tree(tree))
  tree
}


#' Get node IDs from tree
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-05-31
#' @param tree (phylo) a phylogenetic tree.
#' @return (integer) a vector with names being node IDs.
#' @seealso \code{\link{getTree}}
#' @keywords internal
#' @examples
#' m <- matrix(c(0,1,2,3, 1,0,1.5,1.5, 2,1.5,0,1, 3,1.5,1,0), byrow = TRUE,
#'             ncol = 4, dimnames = list(c("run1", "run2", "run3", "run4"),
#'                                       c("run1", "run2", "run3", "run4")))
#' distMat <- as.dist(m)
#' \dontrun{
#' tree <- getTree(distMat)
#' getNodeIDs(tree)
#' }
getNodeIDs <- function(tree){
  nodeIDs <- seq(from = 1L, to = length(tree$tip.label) + tree$Nnode)
  names(nodeIDs) <- c(tree$tip.label, tree$node.label)
  nodeIDs
}


#' Traverses up from leaves to the root
#'
#' @description {
#' While traversing from leaf to root node, at each node a master run is created.
#' Merged features and merged chromatograms from parent runs are estimated. Chromatograms are written on the disk
#' at dataPath/mzml. For each precursor aligned parent time-vectors and corresponding child time-vector
#' are also calculated and written as *_av.rds at dataPath.
#'
#'     Accesors to the new files are added to fileInfo, mzPntrs and prec2chromIndex. Features, reference
#' used for the alignment and adaptiveRTs of global alignments are also added to corresponding environment.
#' }
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-07-01
#' @inheritParams progAlignRuns
#' @param tree (phylo) a phylogenetic tree.
#' @param fileInfo (data-frame) output of \code{\link{getRunNames}}.
#' @param features (list of data-frames) contains features and their properties identified in each run.
#' @param mzPntrs (list) a list of mzRpwiz.
#' @param prec2chromIndex (list) a list of dataframes having following columns: \cr
#' transition_group_id: it is PRECURSOR.ID from osw file. \cr
#' chromatogramIndex: index of chromatogram in mzML file.
#' @param precursors (data-frame) atleast two columns transition_group_id and transition_ids are required.
#' @param adaptiveRTs (environment) For each descendant-pair, it contains the window around the aligned
#'  retention time, within which features with m-score below aligned FDR are considered for quantification.
#' @param refRuns (environment) For each descendant-pair, the reference run is indicated by 1 or 2 for all the peptides.
#' @return (None)
#' @seealso \code{\link{getTree}, \link{getNodeRun}}
#' @keywords internal
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' params <- paramsDIAlignR()
#' fileInfo <- getRunNames(dataPath = dataPath)
#' mzPntrs <- list2env(getMZMLpointers(fileInfo))
#' features <- list2env(getFeatures(fileInfo, maxFdrQuery = params[["maxFdrQuery"]], runType = params[["runType"]]))
#' precursors <- getPrecursors(fileInfo, oswMerged = TRUE, runType = params[["runType"]],
#'  context = "experiment-wide", maxPeptideFdr = params[["maxPeptideFdr"]])
#' precursors <- dplyr::arrange(precursors, .data$peptide_id, .data$transition_group_id)
#' peptideIDs <- unique(precursors$peptide_id)
#' peptideScores <- getPeptideScores(fileInfo, peptideIDs, oswMerged = TRUE, params[["runType"]], params[["context"]])
#' peptideScores <- lapply(peptideIDs, function(pep) dplyr::filter(peptideScores, .data$peptide_id == pep))
#' names(peptideScores) <- as.character(peptideIDs)
#' peptideScores <- list2env(peptideScores)
#' multipeptide <- getMultipeptide(precursors, features)
#' prec2chromIndex <- list2env(getChromatogramIndices(fileInfo, precursors, mzPntrs))
#' adaptiveRTs <- new.env()
#' refRuns <- new.env()
#' tree <- ape::read.tree(text = "(run1:9,(run2:7,run0:2)master2:5)master1;")
#' tree <- ape::reorder.phylo(tree, "postorder")
#' \dontrun{
#' ropenms <- get_ropenms(condaEnv = "envName", useConda=TRUE)
#' multipeptide <- traverseUp(tree, dataPath, fileInfo, features, mzPntrs, prec2chromIndex, precursors, params,
#'  adaptiveRTs, refRuns, multipeptide, peptideScores, ropenms)
#' rm(mzPntrs)
#' # Cleanup
#' file.remove(list.files(dataPath, pattern = "*_av.rds", full.names = TRUE))
#' file.remove(list.files(file.path(dataPath, "mzml"), pattern = "^master[0-9]+\\.chrom\\.mzML$", full.names = TRUE))
#' }
traverseUp <- function(tree, dataPath, fileInfo, features, mzPntrs, prec2chromIndex, precursors,
                       params, adaptiveRTs, refRuns, multipeptide, peptideScores, ropenms, applyFun = lapply){
  vertices <- getNodeIDs(tree)
  ord <- tree$edge[,2] # Traversal order
  num_merge <- length(ord)/2
  merge_start <- 2*(1:num_merge)-1

  # Traverse to the root of all runs.
  for(i in merge_start){
    runA <- names(vertices)[vertices == ord[i]]
    runB <- names(vertices)[vertices == ord[i+1]]
    mergeName <- names(vertices)[vertices == ape::getMRCA(tree, c(ord[i], ord[i+1]))]
    message(runA, " + ", runB, " = ", mergeName)
    multipeptide <- getNodeRun(runA, runB, mergeName, dataPath, fileInfo, features, mzPntrs, prec2chromIndex,
               precursors, params, adaptiveRTs, refRuns, multipeptide, peptideScores, ropenms, applyFun)
  }

  assign("temp", fileInfo, envir = parent.frame(n = 1))
  with(parent.frame(n = 1), fileInfo <- temp)
  message("Created all master runs.")
  multipeptide
}


#' Traverses down from the root to leaves
#'
#' Features of the root node are propagated to all leaves node. Aligned features are set/added in the
#' multipeptide environment.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-07-01
#' @inheritParams traverseUp
#' @param multipeptide (environment) contains multiple data-frames that are collection of features
#'  associated with analytes. This is an output of \code{\link{getMultipeptide}}.
#' @param analytes (integer) this vector contains transition_group_id from precursors. It must be of
#' the same length as of multipeptide.
#' @return (None)
#' @seealso \code{\link{traverseUp}, \link{alignToMaster}}
#' @keywords internal
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' params <- paramsDIAlignR()
#' fileInfo <- getRunNames(dataPath = dataPath)
#' mzPntrs <- list2env(getMZMLpointers(fileInfo))
#' features <- list2env(getFeatures(fileInfo, maxFdrQuery = params[["maxFdrQuery"]], runType = params[["runType"]]))
#' precursors <- getPrecursors(fileInfo, oswMerged = TRUE, runType = params[["runType"]],
#'  context = "experiment-wide", maxPeptideFdr = params[["maxPeptideFdr"]])
#' precursors <- dplyr::arrange(precursors, .data$peptide_id, .data$transition_group_id)
#' peptideIDs <- unique(precursors$peptide_id)
#' peptideScores <- getPeptideScores(fileInfo, peptideIDs, oswMerged = TRUE, params[["runType"]], params[["context"]])
#' peptideScores <- lapply(peptideIDs, function(pep) dplyr::filter(peptideScores, .data$peptide_id == pep))
#' names(peptideScores) <- as.character(peptideIDs)
#' prec2chromIndex <- list2env(getChromatogramIndices(fileInfo, precursors, mzPntrs))
#' multipeptide <- getMultipeptide(precursors, features)
#' adaptiveRTs <- new.env()
#' refRuns <- new.env()
#' tree <- ape::read.tree(text = "(run1:9,(run2:7,run0:2)master2:5)master1;")
#' tree <- ape::reorder.phylo(tree, "postorder")
#' \dontrun{
#' ropenms <- get_ropenms(condaEnv = "envName", useConda=TRUE)
#' multipeptide <- traverseUp(tree, dataPath, fileInfo, features, mzPntrs, prec2chromIndex, precursors, params,
#'  adaptiveRTs, refRuns, multipeptide, peptideScores, ropenms)
#' multipeptide <- getMultipeptide(precursors, features)
#' multipeptide <- traverseDown(tree, dataPath, fileInfo, multipeptide, prec2chromIndex, mzPntrs, precursors,
#'  adaptiveRTs, refRuns, params)
#' # Cleanup
#' rm(mzPntrs)
#' file.remove(list.files(dataPath, pattern = "*_av.rds", full.names = TRUE))
#' file.remove(list.files(file.path(dataPath, "mzml"), pattern = "^master[0-9]+\\.chrom\\.mzML$", full.names = TRUE))
#' }
traverseDown <- function(tree, dataPath, fileInfo, multipeptide, prec2chromIndex, mzPntrs, precursors,
                         adaptiveRTs, refRuns, params, applyFun = lapply){
  vertices <- getNodeIDs(tree)
  ord <- rev(tree$edge[,1])
  num_merge <- length(ord)/2
  junctions <- 2*(1:num_merge)-1

  # set alignment rank for the master1 run.
  master1 <- names(vertices)[vertices == ord[1]]
  peptideIDs <- unique(precursors$peptide_id)
  newMP <- lapply(seq_along(peptideIDs), function(i){
    ##### Set alignment rank in the master1 #####
    peptide <- peptideIDs[i]
    df <- multipeptide[[i]]
    df.other <- df[df$run != master1,]
    df.ref <- df[df$run == master1,]
    refIdx <- which(df.ref[["peak_group_rank"]] == 1)
    refIdx <- refIdx[which.min(df.ref$m_score[refIdx])]
    df.ref[["alignment_rank"]][refIdx] <- 1L

    ##### Set alignment rank for other precursors #####
    idx <- which(precursors$peptide_id == peptide)
    analytes <- precursors[idx, "transition_group_id"]

    if(length(analytes)>1){
      ##### Get XIC_group from reference run. if missing, return unaligned features #####
      chromIndices <- prec2chromIndex[[master1]][["chromatogramIndex"]][idx]
      if(any(is.na(unlist(chromIndices))) | is.null(unlist(chromIndices))){
        warning("Chromatogram indices for peptide ", peptide, " are missing in ", fileInfo[master1, "runName"])
        message("Skipping peptide ", peptide, " across all runs.")
        return(df)
      } else {
        XICs <- lapply(chromIndices, function(iM) extractXIC_group(mz = mzPntrs[[master1]], chromIndices = iM))
        XICs.s <- lapply(XICs, smoothXICs, type = params[["XICfilter"]], kernelLen = params[["kernelLen"]],
                         polyOrd = params[["polyOrd"]])
        names(XICs.s) <- names(XICs) <- as.character(analytes)
      }
      if(params[["smoothPeakArea"]]) XICs <- XICs.s
      df.ref <- setOtherPrecursors(df.ref, XICs, analytes, params)
    }

    ##### Add the new dataframe to the multipeptide #####
    dplyr::bind_rows(df.ref, df.other)
  })
  names(newMP) <- names(multipeptide)

  # Traverse from root to leaf node.
  for(i in junctions){
    pair <- phangorn::Descendants(tree, node = ord[i], type = "children")
    runA <- names(vertices)[vertices == pair[1]]
    runB <- names(vertices)[vertices == pair[2]]
    master <- names(vertices)[vertices == ord[i]]
    message("Mapping peaks from ", master, " to ",  runA, " and ", runB, ".")

    # Get parents to master aligned time vectors.
    filename <- file.path(dataPath, paste0(master, "_av.rds"))
    alignedVecs <- readRDS(file = filename)

    adaptiveRT <- max(adaptiveRTs[[paste(runA, runB, sep = "_")]],
                      adaptiveRTs[[paste(runB, runA, sep = "_")]])

    # Map master to runA
    refA <- refRuns[[master]][,1]
    newMP <- alignToMaster(master, runA, alignedVecs, refA, adaptiveRT, newMP,
                           prec2chromIndex, mzPntrs, fileInfo, precursors, params, applyFun)

    # Map master to runB
    refB <- as.integer(!(refA-1))+1L
    newMP <- alignToMaster(master, runB, alignedVecs, refB, adaptiveRT, newMP,
                           prec2chromIndex, mzPntrs, fileInfo, precursors, params, applyFun)
  }

  # Done
  message("master1 run has been propagated to all parents.")
  newMP
}


#' Align descendants to master
#'
#' During traverse-down, parent runs are aligned to the master/child run. This function performs the
#' alignment by already saved aligned parent-child time vectors. For the aligned peaks, alignment_rank is
#' set to 1 in multipeptide environment.
#'
#' refRun is flipped if eXp is runB instead of runA.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-07-19
#' @import tibble
#' @inherit alignToMaster params
#' @inheritParams traverseDown
#' @inheritParams getChildXICs
#' @param ref (string) name of the descendant run. Must be in the rownames of fileInfo.
#' @param eXp (string) name of one of the parent run. Must be in the rownames of fileInfo.
#' @param alignedVecs (list of dataframes) Each dataframe contains aligned parents time-vectors and
#'  resulting child/master time vector for a precursor. This is the second element of
#'   \code{\link{getChildXICs}} output.
#' @return (None)
#' @seealso \code{\link{traverseUp}, \link{traverseDown}, \link{setAlignmentRank}}
#' @keywords internal
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' params <- paramsDIAlignR()
#' fileInfo <- DIAlignR::getRunNames(dataPath = dataPath)
#' mzPntrs <- list2env(getMZMLpointers(fileInfo))
#' features <- list2env(getFeatures(fileInfo, maxFdrQuery = 0.05, runType = "DIA_proteomics"))
#' precursors <- getPrecursors(fileInfo, oswMerged = TRUE, runType = params[["runType"]],
#'  context = "experiment-wide", maxPeptideFdr = params[["maxPeptideFdr"]])
#' precursors <- dplyr::arrange(precursors, .data$peptide_id, .data$transition_group_id)
#' peptideIDs <- unique(precursors$peptide_id)
#' peptideScores <- getPeptideScores(fileInfo, peptideIDs, oswMerged = TRUE, params[["runType"]], params[["context"]])
#' peptideScores <- lapply(peptideIDs, function(pep) dplyr::filter(peptideScores, .data$peptide_id == pep))
#' names(peptideScores) <- as.character(peptideIDs)
#' prec2chromIndex <- list2env(getChromatogramIndices(fileInfo, precursors, mzPntrs))
#' multipeptide <- getMultipeptide(precursors, features)
#' prec2chromIndex <- list2env(getChromatogramIndices(fileInfo, precursors, mzPntrs))
#' adaptiveRTs <- new.env()
#' refRuns <- new.env()
#' tree <- ape::reorder.phylo(ape::read.tree(text = "(run1:7,run2:2)master1;"), "postorder")
#' \dontrun{
#' ropenms <- get_ropenms(condaEnv = "envName", useConda=TRUE)
#' multipeptide <- traverseUp(tree, dataPath, fileInfo, features, mzPntrs, prec2chromIndex, precursors, params,
#'  adaptiveRTs, refRuns, multipeptide, peptideScores, ropenms)
#' multipeptide <- getMultipeptide(precursors, features)
#' alignedVecs <- readRDS(file = file.path(dataPath, "master1_av.rds"))
#' adaptiveRT <- (adaptiveRTs[["run1_run2"]] + adaptiveRTs[["run2_run1"]])/2
#' multipeptide[["14383"]]$alignment_rank[multipeptide[["14383"]]$run == "master1"] <- 1L
#' multipeptide <- alignToMaster(ref = "master1", eXp = "run1", alignedVecs, refRuns[["master1"]][,1], adaptiveRT,
#'  multipeptide, prec2chromIndex, mzPntrs, fileInfo, precursors, params)
#' # Cleanup
#' rm(mzPntrs)
#' file.remove(file.path(dataPath, "master1_av.rds"))
#' file.remove(file.path(dataPath, "mzml", "master1.chrom.mzML"))
#' }
alignToMaster <- function(ref, eXp, alignedVecs, refRun, adaptiveRT, multipeptide, prec2chromIndex,
                          mzPntrs, fileInfo, precursors, params, applyFun = lapply){
  peptideIDs <- unique(precursors$peptide_id)

  # Aign each peptide to its parent
  num_of_batch <- ceiling(length(peptideIDs)/params[["batchSize"]])
  mp <- lapply(1:num_of_batch, function(iBatch){
    batchSize <- params[["batchSize"]]
    strt <- ((iBatch-1)*batchSize+1)
    stp <- min((iBatch*batchSize), length(multipeptide))
    ##### Get XICs for the batch from the experiment run #####
    XICs <- lapply(strt:stp, function(rownum){
      ##### Get transition_group_id for that peptideID #####
      idx <- which(precursors$peptide_id == peptideIDs[rownum])
      analytes <- precursors[idx, "transition_group_id"]
      ##### Get XIC_group from runA and runB. If missing, add NULL #####
      chromIndices <- prec2chromIndex[[eXp]][["chromatogramIndex"]][idx]
      nope <- any(is.na(unlist(chromIndices))) | is.null(unlist(chromIndices))
      if(nope) return(NULL)
      xics <- lapply(chromIndices, function(i1) mzR::chromatograms(mzPntrs[[eXp]], i1))
      names(xics) <- as.character(analytes)
      xics
    })

    # Align each peptide in multipeptide to its parent
    result <- lapply(strt:stp, function(i){
      peptide <- peptideIDs[i]
      idx <- (i - (iBatch-1)*batchSize)
      alignedVec <- alignedVecs[[i]]
      df <- multipeptide[[i]]

      if(is.null(alignedVec)){
        # Try to  set alignment rank without chromatogram
        return(df)
      } else {
        tAligned <- alignedVec[, c(5L, 2+refRun[i])]
        colnames(tAligned) <- c("tAligned.ref", "tAligned.eXp")
      }

      ##### Get XIC_group from reference run. if missing, return unaligned features #####
      XICs.eXp <- XICs[[idx]]
      if(is.null(XICs.eXp)){
        warning("Chromatogram indices for peptide ", peptide, " are missing in ", fileInfo[eXp, "runName"])
        message("Skipping peptide ", peptide, " across all runs.")
        return(df)
      } else {
        XICs.eXp.s <- lapply(XICs.eXp, smoothXICs, type = params[["XICfilter"]], kernelLen = params[["kernelLen"]],
                             polyOrd = params[["polyOrd"]])
        names(XICs.eXp.s) <- names(XICs.eXp)
      }
      if(params[["smoothPeakArea"]]) XICs.eXp <- XICs.eXp.s
      analytes <- as.integer(names(XICs.eXp))

      # Update alignment rank for the eXp.
      df.other <- df[df$run != eXp,]
      tAligned <- list(tAligned[,1], tAligned[,2])
      df.eXp <- setAlignmentRank(df, ref, eXp, tAligned, XICs.eXp, params, adaptiveRT)
      df.eXp <- setOtherPrecursors(df.eXp, XICs.eXp, analytes, params)
      tibble::remove_rownames(dplyr::bind_rows(df.other, df.eXp))
    })
    result
  })
  mp <- unlist(mp, recursive = FALSE)
  names(mp) <- names(multipeptide)
  message(eXp, " has been aligned to ", ref, ".")
  mp
}
