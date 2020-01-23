## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("time", "Transition"))

#' Plot Extracted-ion chromatogram group.
#'
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot ggtitle geom_vline geom_line theme theme_bw aes element_text
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-13
#'
#' @param XIC_group (list) It is a list of dataframe which has two columns. First column is for time
#'  and second column indicates intensity.
#' @param peakAnnot (numeric) Peak-apex time.
#' @param Title (logical) TRUE: name of the list will be displayed as title.
#' @return A plot to the current device.
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' runs <- c("hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
#'  "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
#' XICs <- getXICs(analytes = "QFNNTDIVLLEDFQK_3", runs = runs, dataPath = dataPath,
#'  XICfilter = "none")
#' plotXICgroup(XICs[["hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"]][[1]])
#'
#' XICs <- getXICs(analytes = "14299_QFNNTDIVLLEDFQK/3", runs = runs, dataPath = dataPath,
#'        XICfilter = "sgolay", SgolayFiltOrd = 4, SgolayFiltLen = 13, analyteInGroupLabel = TRUE)
#' plotXICgroup(XICs[["hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"]][[1]])
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

#' Plot extracted-ion chromatogram.
#'
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-13
#'
#' @param analyte (string) An analyte is as PRECURSOR.GROUP_LABEL or as PEPTIDE.MODIFIED_SEQUENCE and PRECURSOR.CHARGE from osw file.
#' @param run (string) Name of a mzml file without extension.
#' @param dataPath (char) Path to mzml and osw directory.
#' @param maxFdrQuery (numeric) A numeric value between 0 and 1. It is used to filter features from osw file which have SCORE_MS2.QVALUE less than itself.
#' @param XICfilter (string) This must be one of the strings "sgolay", "none".
#' @param SgolayFiltOrd (integer) It defines the polynomial order of filer.
#' @param SgolayFiltLen (integer) Must be an odd number. It defines the length of filter.
#' @param runType (char) This must be one of the strings "DIA_proteomics", "DIA_Metabolomics".
#' @param oswMerged (logical) TRUE for experiment-wide FDR and FALSE for run-specific FDR by pyprophet.
#' @param nameCutPattern (string) regex expression to fetch mzML file name from RUN.FILENAME columns of osw files.
#' @param analyteInGroupLabel (logical) TRUE for getting analytes as PRECURSOR.GROUP_LABEL from osw file.
#' @param peakAnnot (numeric) Peak-apex time.
#' @param Title (logical) TRUE: name of the list will be displayed as title.
#'
#' @return A plot to the current device.
#'
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' run <- "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"
#' plotAnalyteXICs(analyte = "QFNNTDIVLLEDFQK_3", run, dataPath = dataPath, XICfilter = "none")
#' plotAnalyteXICs(analyte = "14299_QFNNTDIVLLEDFQK/3", run, dataPath = dataPath,
#' XICfilter = "sgolay", analyteInGroupLabel = TRUE)
#' @export
plotAnalyteXICs <- function(analyte, run, dataPath = ".", maxFdrQuery = 1.0,
                            XICfilter = "sgolay", SgolayFiltOrd = 4, SgolayFiltLen = 9,
                            runType = "DIA_proteomics", oswMerged = TRUE, nameCutPattern = "(.*)(/)(.*)",
                            analyteInGroupLabel = FALSE, peakAnnot = NULL, Title = NULL){
  if((length(run) != 1) | (length(analyte) != 1)){
    return(stop("One analyte and single run are needed."))
  }
  XICs <- getXICs(analytes = analyte, runs = run, dataPath = dataPath, maxFdrQuery = maxFdrQuery,
                  XICfilter = XICfilter, SgolayFiltOrd = SgolayFiltOrd, SgolayFiltLen = SgolayFiltLen,
                  runType = runType, oswMerged = oswMerged, nameCutPattern = nameCutPattern, analyteInGroupLabel = analyteInGroupLabel)
  plotXICgroup(XICs[[run]][[analyte]], peakAnnot, Title)
}

#' Plot an aligned XIC-group.
#'
#' @details
#' x-axis cannot have the same time-values, therefore, x-axis is indecized.
#'
#' @importFrom zoo na.locf
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-13
#' @param XIC_group (list) It is a list of dataframe which has two columns. First column is for time
#'  and second column indicates intensity.
#' @param idx (integer) Indices of aligned chromatograms.
#' @param peakAnnot (numeric) Peak-apex time.
#' @return A plot to the current device.
#' @keywords internal
plotSingleAlignedChrom <- function(XIC_group, idx, peakAnnot = NULL){
  intensity <- list()
  # Update intensities with aligned time indices.
  for(k in seq_along(XIC_group)){
    mutateInt <- XIC_group[[k]][idx, 2]
    mutateInt <- na.locf(na.locf(mutateInt, na.rm = FALSE),fromLast = TRUE)
    intensity[[k]] <- mutateInt
  }
  #TODO: interpolate mutateT so that it can be plotted on x-axis.
  mutateT <- mapIdxToTime(XIC_group[[1]][, "time"], idx)
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

#' Plot aligned XICs group for a specific peptide.
#'
#' @description
#' AlignObj is the output from getAlignObjs fucntion. This function prepares ggplot objects from AlignObj.
#'
#' @importFrom ggplot2 geom_vline xlab scale_y_continuous
#' @importFrom scales scientific_format
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-13
#'
#' @param AlignObj (S4 object)
#' @param XICs.ref (list) List of extracted ion chromatograms (dataframe) from reference run. The dataframe has two columns: first column is for time
#'  and second column indicates intensity.
#' @param XICs.eXp (list) List of extracted ion chromatograms (dataframe) from experiment run.The dataframe has two columns: first column is for time
#'  and second column indicates intensity.
#' @param refPeakLabel (numeric vector) It contains peak apex, left width and right width.
#' @param annotatePeak (logical) TRUE: Peak boundaries and apex will be highlighted.
#' @return A plot to the current device.
#'
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' runs <- c("hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
#'  "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
#' AlignObjOutput <- getAlignObjs(analytes = "QFNNTDIVLLEDFQK_3", runs, dataPath = dataPath)
#' AlignObj <- AlignObjOutput[["QFNNTDIVLLEDFQK_3"]][[1]]
#' XICs.ref <- AlignObjOutput[["QFNNTDIVLLEDFQK_3"]][[2]]
#' XICs.eXp <- AlignObjOutput[["QFNNTDIVLLEDFQK_3"]][[3]]
#' refPeakLabel <- AlignObjOutput[["QFNNTDIVLLEDFQK_3"]][[4]]
#' @keywords internal
getAlignedFigs <- function(AlignObj, XICs.ref, XICs.eXp, refPeakLabel,
                               annotatePeak = FALSE){
  AlignedIndices <- cbind(AlignObj@indexA_aligned, AlignObj@indexB_aligned,
                          AlignObj@score)
  colnames(AlignedIndices) <- c("indexAligned.ref", "indexAligned.eXp", "score")
  AlignedIndices <- AlignedIndices[(AlignedIndices[,"indexAligned.ref"] != 0L), ]
  AlignedIndices[, 1:2][AlignedIndices[, 1:2] == 0] <- NA
  t.ref <- XICs.ref[[1]][["time"]]
  t.eXp <- mapIdxToTime(XICs.eXp[[1]][["time"]], AlignedIndices[,"indexAligned.eXp"])
  ###################### Plot unaligned chromatogram ######################################
  prefU <- plotXICgroup(XICs.ref) + scale_y_continuous(labels = scientific_format(digits = 1)) + xlab("ref time")
  if(annotatePeak){
    prefU <- prefU +
      geom_vline(xintercept=refPeakLabel$RT[1], lty="dotted", size = 0.3) +
      geom_vline(xintercept=refPeakLabel$leftWidth[1], lty="dashed", size = 0.1) +
      geom_vline(xintercept=refPeakLabel$rightWidth[1], lty="dashed", size = 0.1)
  }

  peXpU <- plotXICgroup(XICs.eXp) + scale_y_continuous(labels = scientific_format(digits = 1)) + xlab("eXp time")
  if(annotatePeak){
    peXpU <- peXpU +
      geom_vline(xintercept=t.eXp[which.min(abs(t.ref - refPeakLabel$RT[1]))], lty="dotted", size = 0.3) +
      geom_vline(xintercept=t.eXp[which.min(abs(t.ref - refPeakLabel$leftWidth[1]))], lty="dashed", size = 0.1) +
      geom_vline(xintercept=t.eXp[which.min(abs(t.ref - refPeakLabel$rightWidth[1]))], lty="dashed", size = 0.1)
  }

  ###################### Plot aligned chromatogram ######################################
  peXpA <- plotSingleAlignedChrom(XICs.eXp, idx = AlignedIndices[,"indexAligned.eXp"]) +
    scale_y_continuous(labels = scientific_format(digits = 1)) + xlab("eXp Aligned index")
  if(annotatePeak){
    peXpA <- peXpA +
      geom_vline(xintercept=which.min(abs(t.ref - refPeakLabel$RT[1])),
                 lty="dotted", size = 0.3) +
      geom_vline(xintercept=which.min(abs(t.ref - refPeakLabel$leftWidth[1])),
                 lty="dashed", size = 0.1) +
      geom_vline(xintercept=which.min(abs(t.ref - refPeakLabel$rightWidth[1])),
                 lty="dashed", size = 0.1)
  }

  ###################### return ggplot objects ######################################
  figs <- list("prefU" = prefU,"peXpU" = peXpU, "peXpA" = peXpA)
  figs
}

#' Plot aligned XICs group for a specific peptide.
#' AlignObjOutput is the output from getAlignObjs fucntion.
#'
#' @importFrom gridExtra grid.arrange
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-13
#'
#' @param AlignObjOutput (list) The list contains AlignObj, raw XICs for reference and experiment, and reference-peak label.
#' @param plotType This must be one of the strings "All", "onlyUnaligned" and "onlyAligned".
#' @param DrawAlignR (logical) TRUE: ggplot objects will be returned.
#' @param annotatePeak (logical) TRUE: Peak boundaries and apex will be highlighted.
#' @param saveFigs (logical) TRUE: Figures will be saved in AlignedAnalytes.pdf .
#' @return A plot to the current device.
#'
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' runs <- c("hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
#'  "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
#' AlignObjOutput <- getAlignObjs(analytes = "QFNNTDIVLLEDFQK_3", runs, dataPath = dataPath)
#' plotAlignedAnalytes(AlignObjOutput)
#' @export
plotAlignedAnalytes <- function(AlignObjOutput, plotType = "All", DrawAlignR = FALSE,
                                annotatePeak = FALSE, saveFigs = FALSE){
  if((length(AlignObjOutput) > 1) | saveFigs){
    grDevices::pdf("AlignedAnalytes.pdf")
  }
  for(i in seq_along(AlignObjOutput)){
    if(is.null(AlignObjOutput[[i]])){
      next
    }
    AlignObj <- AlignObjOutput[[i]][[1]]
    XICs.ref <- AlignObjOutput[[i]][[2]]
    XICs.eXp <- AlignObjOutput[[i]][[3]]
    refPeakLabel <- AlignObjOutput[[i]][[4]]
    analyte <- names(AlignObjOutput)[i]
    refRun <- names(AlignObjOutput[[i]])[2]
    eXpRun <- names(AlignObjOutput[[i]])[3]
    figs <- getAlignedFigs(AlignObj, XICs.ref, XICs.eXp, refPeakLabel, annotatePeak)

    if(DrawAlignR){
      return(figs)}

    if(plotType == "onlyAligned"){
      grid.arrange(figs[["prefU"]], figs[["peXpA"]], nrow=2, ncol=1,
                   top = paste0(analyte,"\n", "ref: ", refRun, "\n", "eXp: ", eXpRun ))
    } else if(plotType == "onlyUnaligned"){
      grid.arrange(figs[["prefU"]], figs[["peXpU"]], nrow=2, ncol=1,
                   top = paste0(analyte,"\n", "ref: ", refRun, "\n", "eXp: ", eXpRun ))
    } else{
      grid.arrange(figs[["peXpU"]], figs[["prefU"]], figs[["peXpA"]],
                   nrow=3, ncol=1, top = paste0(analyte,"\n", "ref: ", refRun, "\n", "eXp: ", eXpRun ))
    }
  }
  if((length(AlignObjOutput) > 1) | saveFigs){
    grDevices::dev.off()
  }
}


#' Visualize alignment path through similarity matrix
#'
#' Plot aligned path through the similarity matrix. Reference run has indices on X-axis, eXp run has them on Y-axis.
#' In getAlignObjs function, objType must be set to medium.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-13
#' @param AlignObjOutput (list) The list contains AlignObj, raw XICs for reference and experiment, and reference-peak label.
#' @return A plot to the current device.
#'
#' @examples
#' library(lattice)
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' runs <- c("hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
#'  "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
#' AlignObjOutput <- getAlignObjs(analytes = "QFNNTDIVLLEDFQK_3", runs, dataPath = dataPath,
#'  objType = "medium")
#' plotAlignmentPath(AlignObjOutput)
#' @export
plotAlignmentPath <- function(AlignObjOutput){
  Alignobj <- AlignObjOutput[[1]][[1]]
  analyte <- names(AlignObjOutput)[1]
  s <- Alignobj@s
  Path <- Alignobj@path[2:nrow(Alignobj@path), 2:ncol(Alignobj@path)]
  lattice::levelplot(s, axes = TRUE, xlab = "ref index", ylab = "eXp index",
            main = paste0("Hybrid alignment through the similarity matrix\n for ",
                          analyte), fontsize = 7) +
    latticeExtra::as.layer(lattice::levelplot(Path, col.regions = c("transparent", "green"),
                                     alpha = 1, axes = FALSE))
}
