#' Plot Extracted-ion chromatogram group.
#'
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot ggtitle geom_vline geom_line theme theme_bw aes element_text xlab ylab
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-13
#' @importFrom rlang .data
#' @param XIC_group (list) It is a list of dataframe which has two columns. First column is for time
#'  and second column indicates intensity.
#' @param peakAnnot (numeric) Peak-apex time.
#' @param Title (logical) TRUE: name of the list will be displayed as title.
#' @return A plot to the current device.
#' @seealso \code{\link{plotAnalyteXICs}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' runs <- c("hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
#'  "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
#' XICs <- getXICs(analytes = 4618L, runs = runs, dataPath = dataPath, oswMerged = TRUE)
#' plotXICgroup(XICs[["hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]])
#' XICs <- smoothXICs(XICs[["hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]],
#'   type = "sgolay", kernelLen = 13, polyOrd = 4)
#' plotXICgroup(XICs, Title = "Precursor 4618 \n
#'  run hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt")
#' @export
plotXICgroup <- function(XIC_group, peakAnnot = NULL, Title =NULL){
  df <- do.call("cbind", lapply(XIC_group, `[`, i =, j = 2))
  df <- data.frame(cbind(XIC_group[[1]][,1],  df))
  colnames(df) <- c("time", paste("V", 1:(ncol(df)-1), sep=""))
  df <- gather(df, key = "Transition", value = "Intensity", -.data$time)
  g <- ggplot(df, aes(.data$time, .data$Intensity, col=.data$Transition)) +
    geom_line(show.legend = FALSE) + xlab("time") + ylab("intensity")+ theme_bw()
  if(!is.null(Title)) g <- g + ggtitle(paste0(Title)) + theme(plot.title = element_text(hjust = 0.5))
  if(!is.null(peakAnnot)){
    g <- g + geom_vline(xintercept=peakAnnot, lty="dotted", size = 0.4)
  }
  return(g)
}


#' Plot extracted-ion chromatogram.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-13
#'
#' @param analyte (integer) an analyte is a PRECURSOR.ID from the osw file.
#' @param run (string) Name of a mzml file without extension.
#' @param dataPath (string) path to mzml and osw directory.
#' @param maxFdrQuery (numeric) A numeric value between 0 and 1. It is used to filter features from osw file which have SCORE_MS2.QVALUE less than itself.
#' @param XICfilter (string) must be either sgolay, boxcar, gaussian, loess or none.
#' @param polyOrd (integer) order of the polynomial to be fit in the kernel.
#' @param kernelLen (integer) number of data-points to consider in the kernel.
#' @param runType (char) This must be one of the strings "DIA_proteomics", "DIA_Metabolomics".
#' @param oswMerged (logical) TRUE for experiment-wide FDR and FALSE for run-specific FDR by pyprophet.
#' @param peakAnnot (numeric) Peak-apex time.
#' @param Title (logical) TRUE: name of the list will be displayed as title.
#'
#' @return A plot to the current device.
#' @seealso \code{\link{plotXICgroup}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' run <- "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"
#' plotAnalyteXICs(analyte = 2474L, run, dataPath = dataPath, oswMerged = TRUE, XICfilter = "none")
#' plotAnalyteXICs(analyte = 2474L, run, dataPath = dataPath, oswMerged = TRUE, XICfilter = "sgolay")
#' @export
plotAnalyteXICs <- function(analyte, run, dataPath = ".", maxFdrQuery = 1.0,
                            XICfilter = "sgolay", polyOrd = 4, kernelLen = 9,
                            runType = "DIA_proteomics", oswMerged = TRUE,
                            peakAnnot = NULL, Title = NULL){
  if((length(run) != 1) | (length(analyte) != 1)){
    return(stop("One analyte and single run are needed."))
  }
  XICs <- getXICs(analytes = analyte, runs = run, dataPath = dataPath, maxFdrQuery = maxFdrQuery,
                  runType = runType, oswMerged = oswMerged)
  XICs <- smoothXICs(XICs[[run]][[1]], type=XICfilter, kernelLen = kernelLen, polyOrd = polyOrd)
  plotXICgroup(XICs, peakAnnot, Title)
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
    mutateInt <- na.locf(na.locf(na.approx(mutateInt, na.rm = FALSE), na.rm = FALSE), fromLast = TRUE)
    intensity[[k]] <- mutateInt
  }
  mutateT <- mapIdxToTime(XIC_group[[1]][, "time"], idx)

  # Extrapolate time
  if(sum(is.na(mutateT)) > 0){
    df <- data.frame(x = which(!is.na(mutateT)), y = mutateT[!is.na(mutateT)] )
    fit <- lm(y ~ x, df)
    df <- data.frame(x = which(is.na(mutateT)), y = NA)
    df$y <- predict(fit, newdata = df)
    for(i in nrow(df)) mutateT[df$x[i]]  <- df$y[i]
  }

  df <- do.call("cbind", intensity)
  df <- cbind(mutateT, as.data.frame(df))
  df <- gather(df, key = "Transition", value = "Intensity", -mutateT)
  # Plot chromatogram
  g <- ggplot(df, aes(mutateT, .data$Intensity, col=.data$Transition)) + geom_line(show.legend = FALSE) + theme_bw()
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
#' @keywords internal
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' runs <- c("hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
#'  "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
#' AlignObjOutput <- getAlignObjs(analytes = 4618L, runs, dataPath = dataPath)
#' AlignObj <- AlignObjOutput[[2]][["4618"]][[1]][["AlignObj"]]
#' XICs.ref <- AlignObjOutput[[2]][["4618"]][[1]][["ref"]]
#' XICs.eXp <- AlignObjOutput[[2]][["4618"]][[1]][["eXp"]]
#' refPeakLabel <- AlignObjOutput[[2]][["4618"]][[1]][["peak"]]
#' \dontrun{
#' getAlignedFigs(AlignObj, XICs.ref, XICs.eXp, refPeakLabel)
#' }
getAlignedFigs <- function(AlignObj, XICs.ref, XICs.eXp, refPeakLabel,
                               annotatePeak = FALSE){
  AlignedIndices <- cbind(slot(AlignObj, "indexA_aligned"),
                          slot(AlignObj, "indexB_aligned"))
  colnames(AlignedIndices) <- c("indexAligned.ref", "indexAligned.eXp")
  # Do not include gaps in reference run.
  AlignedIndices <- AlignedIndices[(AlignedIndices[,"indexAligned.ref"] != 0L), ]
  AlignedIndices[, 1:2][AlignedIndices[, 1:2] == 0] <- NA
  t.ref <- XICs.ref[[1]][, "time"]
  t.eXp <- mapIdxToTime(XICs.eXp[[1]][, "time"], AlignedIndices[,"indexAligned.eXp"])
  ###################### Plot unaligned chromatogram ######################################
  prefU <- plotXICgroup(XICs.ref) + scale_y_continuous(labels = scales::scientific_format(digits = 1)) + xlab("ref time")
  if(annotatePeak){
    prefU <- prefU +
      geom_vline(xintercept=refPeakLabel$RT[1], lty="dotted", size = 0.3) +
      geom_vline(xintercept=refPeakLabel$leftWidth[1], lty="dashed", size = 0.1) +
      geom_vline(xintercept=refPeakLabel$rightWidth[1], lty="dashed", size = 0.1)
  }

  peXpU <- plotXICgroup(XICs.eXp) + scale_y_continuous(labels = scales::scientific_format(digits = 1)) + xlab("eXp unaligned time")
  if(annotatePeak){
    peXpU <- peXpU +
      geom_vline(xintercept=t.eXp[which.min(abs(t.ref - refPeakLabel$RT[1]))], lty="dotted", size = 0.3) +
      geom_vline(xintercept=t.eXp[which.min(abs(t.ref - refPeakLabel$leftWidth[1]))], lty="dashed", size = 0.1) +
      geom_vline(xintercept=t.eXp[which.min(abs(t.ref - refPeakLabel$rightWidth[1]))], lty="dashed", size = 0.1)
  }

  ###################### Plot aligned chromatogram ######################################
  peXpA <- plotSingleAlignedChrom(XICs.eXp, idx = AlignedIndices[,"indexAligned.eXp"]) +
    scale_y_continuous(labels = scales::scientific_format(digits = 1)) + xlab("eXp aligned index")
  if(annotatePeak){
    peXpA <- peXpA +
      geom_vline(xintercept=t.eXp[which.min(abs(t.ref - refPeakLabel$RT[1]))], lty="dotted", size = 0.3) +
      geom_vline(xintercept=t.eXp[which.min(abs(t.ref - refPeakLabel$leftWidth[1]))], lty="dashed", size = 0.1) +
      geom_vline(xintercept=t.eXp[which.min(abs(t.ref - refPeakLabel$rightWidth[1]))], lty="dashed", size = 0.1)
  }

  ###################### return ggplot objects ######################################
  figs <- list("prefU" = prefU,"peXpU" = peXpU, "peXpA" = peXpA)
  figs
}


#' Plot aligned XICs group for a specific peptide.
#' AlignObjOutput is the output from getAlignObjs fucntion.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-13
#'
#' @param AlignObjOutput (list) list contains fileInfo, AlignObj, raw XICs for reference and experiment, and reference-peak label.
#' @param plotType (string) must be one of the strings "All", "onlyUnaligned" and "onlyAligned".
#' @param outFile (string) name of the output pdf file.
#' @param annotatePeak (logical) TRUE: Peak boundaries and apex will be highlighted.
#' @param saveFigs (logical) TRUE: Figures will be saved in AlignedAnalytes.pdf .
#' @return A plot to the current device.
#'
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' runs <- c("hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
#'  "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
#' AlignObjOutput <- getAlignObjs(analytes = 4618L, runs, dataPath = dataPath)
#' plotAlignedAnalytes(AlignObjOutput)
#' @export
plotAlignedAnalytes <- function(AlignObjOutput, plotType = "All", outFile = "AlignedAnalytes.pdf",
                                annotatePeak = FALSE, saveFigs = FALSE){
  if((length(AlignObjOutput[[2]][[1]]) > 1) | length(AlignObjOutput[[2]]) > 1 | saveFigs){
    grDevices::pdf(outFile)
  }
  # Get fileInfo (output of getRunNames)
  vec <- AlignObjOutput[[1]][,"runName"]
  names(vec) <- rownames(AlignObjOutput[[1]])
  # Iterate over each precursor, check if it is NULL.
  for(i in seq_along(AlignObjOutput[[2]])){
    if(is.null(AlignObjOutput[[2]][[i]])){
      next
    }
    pairs <- names(AlignObjOutput[[2]][[i]])
    for(pair in pairs){
      refRun <- strsplit(pair, split = "_")[[1]][1]
      eXpRun <- strsplit(pair, split = "_")[[1]][2]
      AlignObj <- AlignObjOutput[[2]][[i]][[pair]][["AlignObj"]]
      XICs.ref <- AlignObjOutput[[2]][[i]][[pair]][["ref"]]
      XICs.eXp <- AlignObjOutput[[2]][[i]][[pair]][["eXp"]]
      refPeakLabel <- AlignObjOutput[[2]][[i]][[pair]][["peak" ]]

      analyte <- names(AlignObjOutput[[2]])[i]
      figs <- getAlignedFigs(AlignObj, XICs.ref, XICs.eXp, refPeakLabel, annotatePeak)
      if(plotType == "onlyAligned"){
        gridExtra::grid.arrange(figs[["prefU"]], figs[["peXpA"]], nrow=2, ncol=1,
                     top = paste0(analyte,"\n", "ref: ", refRun, "\n", "eXp: ", eXpRun ))
      } else if(plotType == "onlyUnaligned"){
        gridExtra::grid.arrange(figs[["prefU"]], figs[["peXpU"]], nrow=2, ncol=1,
                     top = paste0(analyte,"\n", "ref: ", refRun, "\n", "eXp: ", eXpRun ))
      } else{
        gridExtra::grid.arrange(figs[["prefU"]], figs[["peXpA"]], figs[["peXpU"]],
                     nrow=3, ncol=1, top = paste0(analyte,"\n", "ref: ", vec[[refRun]], "\n", "eXp: ", vec[[eXpRun]] ))
      }
    }
  }
  if((length(AlignObjOutput[[2]][[1]]) > 1) | length(AlignObjOutput[[2]]) > 1 | saveFigs){
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
#' AlignObjOutput <- getAlignObjs(analytes = 4618L, runs, dataPath = dataPath, objType = "medium")
#' plotAlignmentPath(AlignObjOutput)
#' @export
plotAlignmentPath <- function(AlignObjOutput){
  vec <- AlignObjOutput[[1]][,"runName"]
  names(vec) <- rownames(AlignObjOutput[[1]])
  pair <- names(AlignObjOutput[[2]][[1]])
  ref <- strsplit(pair, split = "_")[[1]][1]
  eXp <- strsplit(pair, split = "_")[[1]][2]
  AlignObj <- AlignObjOutput[[2]][[1]][[pair]][["AlignObj"]]
  analyte <- names(AlignObjOutput[[2]])[1]

  s <- slot(AlignObj, "s")
  Path <- slot(AlignObj, "path")
  Path <- Path[2:nrow(Path), 2:ncol(Path)]
  lattice::levelplot(s, axes = TRUE, xlab = vec[[ref]], ylab = vec[[eXp]],
            main = paste0("Hybrid alignment through the similarity matrix\n for ",
                          analyte), fontsize = 7) +
    latticeExtra::as.layer(lattice::levelplot(Path, col.regions = c("transparent", "green"),
                                     alpha = 1, axes = FALSE))
}
