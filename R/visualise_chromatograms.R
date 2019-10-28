#' Plot Extracted-ion chromatogram group.
#'
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot ggtitle geom_vline geom_line theme theme_bw aes
#' @export
plotXICgroup <- function(XIC_group, peakAnnot = NULL, Title =NULL){
  df <- do.call("cbind", XIC_group)
  df <- df[,!duplicated(colnames(df))]
  colnames(df) <- c("time", paste("V", 1:(ncol(df)-1), sep=""))
  df <- gather(df, key = "Transition", value = "Intensity", -time)
  g <- ggplot(df, aes(time, Intensity, col=Transition)) + geom_line(show.legend = FALSE) + theme_bw()
  if(!is.null(Title)) g <- g + ggtitle(paste0(Title)) + theme(plot.title = element_text(hjust = 0.5))
  if(!is.null(peakAnnot)){
    g <- g + geom_vline(xintercept=peakAnnot, lty="dotted", size = 0.4)
  }
  return(g)
}

#' Plot an aligned XIC-group.
#'
#' @importFrom zoo na.locf
plotSingleAlignedChrom <- function(XIC_group, idx, peakAnnot = NULL){
  intensity <- list()
  # Update intensities with aligned time indices.
  for(k in 1:length(XIC_group)){
    mutateInt <- XIC_group[[k]][idx, 2]
    mutateInt <- na.locf(na.locf(mutateInt, na.rm = FALSE),fromLast = TRUE)
    intensity[[k]] <- mutateInt
  }
  df <- do.call("cbind", intensity)
  Index <- 1:nrow(df)
  df <- cbind(Index, as.data.frame(df))
  df <- gather(df, key = "Transition", value = "Intensity", -Index)
  # Plot chromatogram
  g <- ggplot(df, aes(Index, Intensity, col=Transition)) + geom_line(show.legend = FALSE) + theme_bw()
  if(!is.null(peakAnnot)){
    g <- g + geom_vline(xintercept=peakAnnot, lty="dotted", size = 0.4)
  }
  return(g)}

#' Plot Extracted-ion chromatogram group for a specific peptide.
#'
#' @export
plotPeptideXICs <- function(peptides, runs, dataPath = ".", maxFdrQuery = 1.0,
                            SgolayFiltOrd = 4, SgolayFiltLen = 9,
                            query = NULL, oswMerged = TRUE, nameCutPattern = "(.*)(/)(.*)",
                            peakAnnot = NULL, Title =NULL){
  XICs <- getXICs(peptides, runs, dataPath , maxFdrQuery,
                   SgolayFiltOrd, SgolayFiltLen,
                   query, oswMerged, nameCutPattern)
  plotXICgroup(XICs[[1]][[1]], peakAnnot, Title)
}

#' Plot Extracted-ion chromatogram group for a specific peptide from MRM run.
#'
#' @export
plotMRMPeptideXICs <- function(peptides, runs, dataPath = ".", maxFdrQuery = 1.0,
                            SgolayFiltOrd = 2, SgolayFiltLen = 3,
                            query = NULL, oswMerged = FALSE, nameCutPattern = "(.*)(/)(.*)",
                            peakAnnot = NULL, Title =NULL){
  XICs <- getMRMXICs(peptides, runs, dataPath , maxFdrQuery,
                  SgolayFiltOrd, SgolayFiltLen,
                  query, oswMerged, nameCutPattern)
  plotXICgroup(XICs[[1]][[1]], peakAnnot, Title)
}


#' Plot aligned XICs group for a specific peptide.
#' AlignObjOutput is the output from getAlignObjs fucntion.
#'
#' @importFrom ggplot2 geom_vline xlab scale_y_continuous
#' @importFrom gridExtra grid.arrange
#' @importFrom scales scientific_format
#' @export
  plotAlignedPeptides <- function(AlignObjOutput, plotType = "All", annotatePeak = FALSE){
  Alignobj <- AlignObjOutput[[1]][[1]]
  XICs.ref <- AlignObjOutput[[1]][[2]]
  XICs.eXp <- AlignObjOutput[[1]][[3]]
  refPeakLabel <- AlignObjOutput[[1]][[4]]
  AlignedIndices <- cbind(Alignobj@indexA_aligned,
                          Alignobj@indexB_aligned,
                          Alignobj@score)
  colnames(AlignedIndices) <- c("indexAligned.ref", "indexAligned.eXp", "score")
  AlignedIndices[, 1:2][AlignedIndices[, 1:2] == 0] <- NA
  peptide <- names(AlignObjOutput)[1]
  refRun <- names(AlignObjOutput[[1]])[2]
  eXpRun <- names(AlignObjOutput[[1]])[3]
  t.ref <- mapIdxToTime(XICs.ref[[1]][[1]], AlignedIndices[,"indexAligned.ref"])
  t.eXp <- mapIdxToTime(XICs.eXp[[1]][[1]], AlignedIndices[,"indexAligned.eXp"])

  ###################### Plot unaligned chromatogram ######################################
  pTL <- plotXICgroup(XICs.ref) + scale_y_continuous(labels = scientific_format(digits = 1))
  if(annotatePeak){
    pTL <- pTL +
      geom_vline(xintercept=refPeakLabel$RT[1], lty="dotted", size = 0.3) +
      geom_vline(xintercept=refPeakLabel$leftWidth[1], lty="dashed", size = 0.1) +
      geom_vline(xintercept=refPeakLabel$rightWidth[1], lty="dashed", size = 0.1)
  }

  pBL <- plotXICgroup(XICs.eXp) + scale_y_continuous(labels = scientific_format(digits = 1))
  if(annotatePeak){
    pBL <- pBL +
      geom_vline(xintercept=t.eXp[which.min(abs(t.ref - refPeakLabel$RT[1]))], lty="dotted", size = 0.3) +
      geom_vline(xintercept=t.eXp[which.min(abs(t.ref - refPeakLabel$leftWidth[1]))], lty="dashed", size = 0.1) +
      geom_vline(xintercept=t.eXp[which.min(abs(t.ref - refPeakLabel$rightWidth[1]))], lty="dashed", size = 0.1)
  }

  ###################### Plot aligned chromatogram ######################################
  pTR <- plotSingleAlignedChrom(XICs.ref, idx = AlignedIndices[,"indexAligned.ref"]) +
    scale_y_continuous(labels = scientific_format(digits = 1)) + xlab("ref Index")
  if(annotatePeak){
    pTR <- pTR +
      geom_vline(xintercept=which.min(abs(t.ref - refPeakLabel$RT[1])),
                 lty="dotted", size = 0.3) +
      geom_vline(xintercept=which.min(abs(t.ref - refPeakLabel$leftWidth[1])),
                 lty="dashed", size = 0.1) +
      geom_vline(xintercept=which.min(abs(t.ref - refPeakLabel$rightWidth[1])),
                 lty="dashed", size = 0.1)
  }

  pBR <- plotSingleAlignedChrom(XICs.eXp, idx = AlignedIndices[,"indexAligned.eXp"]) +
    scale_y_continuous(labels = scientific_format(digits = 1)) + xlab("eXp Index")
  if(annotatePeak){
    pBR <- pBR +
      geom_vline(xintercept=which.min(abs(t.ref - refPeakLabel$RT[1])),
               lty="dotted", size = 0.3) +
      geom_vline(xintercept=which.min(abs(t.ref - refPeakLabel$leftWidth[1])),
               lty="dashed", size = 0.1) +
      geom_vline(xintercept=which.min(abs(t.ref - refPeakLabel$rightWidth[1])),
               lty="dashed", size = 0.1)
  }

  if(plotType == "onlyAligned"){
    grid.arrange(pTR, pBR, nrow=2, ncol=1, top = paste0(peptide,"\n", "ref: ", refRun,
                                                        "\n", "eXp: ", eXpRun ))
  } else if(plotType == "onlyUnaligned"){
    grid.arrange(pTL, pBL, nrow=2, ncol=1, top = paste0(peptide,"\n", "ref: ", refRun,
                                                        "\n", "eXp: ", eXpRun ))
  } else{
    grid.arrange(pTL, pTR, pBL, pBR, nrow=2, ncol=2, top = paste0(peptide,"\n", "ref: ", refRun,
                                                        "\n", "eXp: ", eXpRun ))
  }
}

#' Plot aligned path through the similarity matrix.
#' Reference run has indices on X-axis, eXp run has them on Y-axis.
#'
#' @importFrom lattice levelplot
#' @importFrom latticeExtra as.layer
#' @export
plotAlignemntPath <- function(AlignObjOutput){
  Alignobj <- AlignObjOutput[[1]][[1]]
  peptide <- names(AlignObjOutput)[1]
  s <- Alignobj@s
  Path <- Alignobj@path[2:nrow(Alignobj@path), 2:ncol(Alignobj@path)]
  levelplot(s, axes = TRUE, xlab = "ref index", ylab = "eXp index",
            main = paste0("Hybrid alignment through the similarity matrix\n for ",
                          peptide), fontsize = 7) +
    as.layer(levelplot(Path, col.regions = c("transparent", "green"),
                                     alpha = 1, axes = FALSE))
}
