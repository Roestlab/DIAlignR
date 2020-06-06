## get the number of descendants for each tip or node:
nr_desc <- function(tree) {
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
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-05-31
#' @import phangorn ape
#' @param distMat (dist) pairwise distance matrix.
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
getTree <- function(distMat){
  tree <- phangorn::upgma(distMat)
  tree <- ape::reorder.phylo(tree, "postorder")
  tree <- ape::makeNodeLabel(tree, method = "number", prefix = "master")
  message("alignment order of runs is as:")
  plot(tree, show.node.label = TRUE)
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
#' distMat <- as.dist(m, diag = FALSE, upper = FALSE)
#' \dontrun{
#' tree <- getTree(distMat)
#' getNodeIDs(tree)
#' }
getNodeIDs <- function(tree){
  nodeIDs <- seq(from = 1L, to = length(tree$tip.label) + tree$Nnode)
  names(nodeIDs) <- c(tree$tip.label, tree$node.label)
  nodeIDs
}

traverseTree <- function(){
  traversalOrder <- tree$edge[,2]
  num_merge <- length(traversalOrder)/2
  merge_start <- 2*(1:num_merge)-1
  for(i in merge_start){
    # Create new list in Feature

    # Get XIC from i XICs.ref
    # Get XIC from i+1 XICs.eXp
    # Get alignedIndices
    # Align to multipeptide
  }
}
