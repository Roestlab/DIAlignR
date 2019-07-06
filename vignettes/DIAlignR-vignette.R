## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----installDIAlignR, eval=FALSE-----------------------------------------
#  require(devtools)
#  install_github("Roestlab/DIAlignR")

## ----loadDIAlignR--------------------------------------------------------
library(DIAlignR)

## ----loadChroms, eval=FALSE----------------------------------------------
#  library(mzR)
#  library(signal)
#  TargetPeptides <- read.table("500Peptide4Alignment.csv", sep = ",", header = T)
#  temp <- list.files(pattern="*.mzML", recursive = TRUE)
#  for(filename in temp){
#    # This makes sure that order of extracted MS2 chromatograms is same for each run.
#    mz <- openMSfile(filename, backend = "pwiz")
#    chromHead <- chromatogramHeader(mz)
#    filename <- gsub("(.*)(/hroest_)(.*)(_SW.chrom.mzML)", replacement = "\\3", filename)
#    chromatogramIndices <- chromHead$chromatogramIndex[match(TargetPeptides$transition_name, chromHead$chromatogramId)]
#    TargetPeptides[filename] <- chromatogramIndices
#    transition_group_ids <- unique(TargetPeptides$transition_group_id)
#    ChromsExtractedPerRun <- sapply(transition_group_ids, function(id){
#      chromIndices <- chromatogramIndices[TargetPeptides$transition_group_id == id]
#      # ChromsExtracted <- lapply(1:length(chromIndices), function(i) chromatograms(mz, chromIndices[i]))
#      ChromsExtracted <- lapply(1:length(chromIndices), function(i) {
#        rawChrom <- chromatograms(mz, chromIndices[i])
#        rawChrom[,2] <- sgolayfilt(rawChrom[,2], p = 4, n = 9) # To smooth chromatograms, use Savitzky-Golay filter
#        return(rawChrom)
#      } )
#      return(ChromsExtracted)
#    })
#    names(ChromsExtractedPerRun) <- transition_group_ids
#    rm(mz)
#    saveRDS(ChromsExtractedPerRun, paste0(filename, "_ChromSelected.rds"))
#  }
#  write.table(TargetPeptides, file = "TargetPeptidesWchromIndex.csv", sep = ",")
#  
#  # Load chromatograms of all runs
#  temp <- list.files(pattern = "*_ChromSelected.rds")
#  StrepChroms1 <- list()
#  for(i in 1:length(temp)){
#    StrepChroms1[[i]] <- readRDS(temp[i])
#  }
#  temp <- sapply(temp, strsplit, split = "_ChromSelected.rds", USE.NAMES = FALSE)
#  names(StrepChroms1) <- temp

## ---- eval=FALSE---------------------------------------------------------
#  runs <- names(StrepChroms);
#  peptides <- names(StrepChroms[["run1"]])

## ----VisualizeAlignment, fig.width=6, fig.align='center', fig.height=4, eval=TRUE, fig.show = "hold"----
library(lattice)
library(ggplot2)
library(reshape2)

plotChromatogram <- function(data, run, peptide, peakAnnot, printTitle =TRUE){
  df <- do.call("cbind", data[[run]][[peptide]])
  df <- df[,!duplicated(colnames(df))]
  df <- melt(df, id.vars="time", value.name = "Intensity")
  g <- ggplot(df, aes(time, Intensity, col=variable)) + geom_line(show.legend = FALSE) + theme_bw()
  if(printTitle) g <- g + ggtitle(paste0(run, ", ",peptide)) + theme(plot.title = element_text(hjust = 0.5))
  g <- g + geom_vline(xintercept=peakAnnot[peptide, run], lty="dotted", size = 0.4)
  return(g)
}

plotChromatogram(StrepChroms, "run3", "3505_LATWYSEMK/2", StrepAnnot)
plotChromatogram(StrepChroms, "run4", "3505_LATWYSEMK/2", StrepAnnot)

## ----globalFit, eval=TRUE------------------------------------------------
run_pair <- "run1_run2"
loess.fit <- getLOESSfit(run_pair, peptides, oswOutStrep, 0.15)
StrepAnnot <- as.data.frame(StrepAnnot)
predict.run2 <- predict(loess.fit, data.frame(RUN1 = StrepAnnot[, "run1"]))
Err_global <- predict.run2 - StrepAnnot[,"run2"]
sum(abs(Err_global))/20
MappedTimeGlobal <- getPepPeakCorp(featureTable =  StrepAnnot, pairName = run_pair, alignType = "global", oswFeatureList = oswOutStrep, spanvalue = 0.15)
Err_global <- StrepAnnot[,"run2"] - MappedTimeGlobal
sum(abs(Err_global))/20

## ----localFit, eval=TRUE-------------------------------------------------
simMeasure <- "dotProductMasked"
pair <- "run1_run2"
MappedTimeLocal <- getPepPeakCorp(StrepAnnot, pair, StrepChroms[["run1"]], StrepChroms[["run2"]], "local", oswOutStrep)
Err_local <- StrepAnnot[,"run2"] - MappedTimeLocal
sum(abs(Err_local))/20

## ----hybridAlignParam----------------------------------------------------
samplingTime <-3.4 # In example dataset, all points are acquired at 3.4 second interval.
meanRSE <- 7.4 # Defines base upper limit of sampling space within that local alignment is done. 
pair <- "run1_run2"
MappedTimeHybrid <- getPepPeakCorp(StrepAnnot, pair, StrepChroms[["run1"]], StrepChroms[["run2"]], "hybrid", oswOutStrep, expRSE = meanRSE, samplingTime = samplingTime)
Err_hybrid <- StrepAnnot[,"run2"] - MappedTimeHybrid
sum(abs(Err_hybrid))/20

## ----plotErr, fig.width=6, fig.align='center', fig.height=6, fig.show='hold', eval=TRUE----
plotErrorCurve <- function(x, clr = "black", SameGraph = FALSE, xmax = 120, ...){
    x <- x[!is.na(x)]
    breaks = seq(0, xmax, by=0.5)
    duration.cut = cut(x, breaks, right = FALSE) 
    duration.freq = table(duration.cut)
    cumfreq0 = c(0, cumsum(duration.freq))
    if(SameGraph == TRUE){lines(breaks, cumfreq0/length(x), col = clr, ...)}
    else{plot(breaks, cumfreq0/length(x), col = clr, type = "l", ...)}
}
plotErrorCurve(abs(Err_global), "blue",xmax = 60, xlab = "Retention time difference (in sec)", ylab = "Cumulative fraction of peptides")
plotErrorCurve(abs(Err_local), "darkgreen", SameGraph = TRUE, xlab = "Retention time difference (in sec)", ylab = "Cumulative fraction of peptides")
plotErrorCurve(abs(Err_hybrid), "red", SameGraph = TRUE, xlab = "Retention time difference (in sec)", ylab = "Cumulative fraction of peptides")
grid(); legend("bottomright",  c("global", "local", "hybrid"), col = c("blue", "darkgreen", "red"), lwd = 1.5, box.lwd = 1)

## ----unalignedPeak, eval = FALSE-----------------------------------------
#  
#  levelplot(s, axes = TRUE, xlab = "run1 index", ylab = "run2 index")
#  Path <- getAlignmentPath(AlignedIndices[[1]], s)
#  levelplot(s, axes = TRUE, xlab = "run1 index", ylab = "run2 index", main = paste0("Alignment path through the ", simMeasure, " similarity matrix\n for ", peptide)) + latticeExtra::as.layer(levelplot(Path, col.regions = c("transparent", "green"), alpha = 1, axes = FALSE))

## ----sessionInfo, eval=TRUE----------------------------------------------
devtools::session_info()

