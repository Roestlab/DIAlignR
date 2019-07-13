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

## ---- echo = FALSE, warning=FALSE----------------------------------------
library(lattice)
library(ggplot2)
library(reshape2)
library(zoo)

plotChromatogram <- function(data, run, peptide, peakAnnot, printTitle =TRUE){
  df <- do.call("cbind", data[[run]][[peptide]])
  df <- df[,!duplicated(colnames(df))]
  df <- melt(df, id.vars="time", value.name = "Intensity")
  g <- ggplot(df, aes(time, Intensity, col=variable)) + geom_line(show.legend = FALSE) + theme_bw()
  if(printTitle) g <- g + ggtitle(paste0(run, ", ",peptide)) + theme(plot.title = element_text(hjust = 0.5))
  g <- g + geom_vline(xintercept=peakAnnot[peptide, run], lty="dotted", size = 0.4)
  return(g)
}
plotSingleAlignedChrom <- function(data, run, peptide, idx, t, printTitle = TRUE){
  intensity <- list()
  for(k in 1:length(data[[run]][[peptide]])){
      mutateInt <- data[[run]][[peptide]][[k]][idx, 2]
      mutateInt <- na.locf(na.locf(mutateInt, na.rm = FALSE),fromLast = TRUE)
      intensity[[k]] <- mutateInt
  }
  df <- do.call("cbind", intensity)
  Index <- 1:nrow(df)
  df <- cbind(Index, as.data.frame(df))
  df <- melt(df, id.vars="Index", value.name = "Intensity")

  g <- ggplot(df, aes(Index, Intensity, col=variable)) + geom_line(show.legend = FALSE) + theme_bw()
  if(printTitle) g <- g + ggtitle(paste0(run, ", ",peptide)) + theme(plot.title = element_text(hjust = 0.5))
  return(g)}

## ----VisualizeAlignment, fig.width=6, fig.align='center', fig.height=4, eval=TRUE, fig.show = "hold"----
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

## ------------------------------------------------------------------------
plotErrorCurve <- function(x, clr = "black", SameGraph = FALSE, xmax = 120, ...){
    x <- x[!is.na(x)]
    breaks = seq(0, xmax, by=0.5)
    duration.cut = cut(x, breaks, right = FALSE) 
    duration.freq = table(duration.cut)
    cumfreq0 = c(0, cumsum(duration.freq))
    if(SameGraph == TRUE){lines(breaks, cumfreq0/length(x), col = clr, ...)}
    else{plot(breaks, cumfreq0/length(x), col = clr, type = "l", ...)}
}

## ----plotErr, fig.width=6, fig.align='center', fig.height=6, fig.show='hold', eval=TRUE----
plotErrorCurve(abs(Err_global), "blue",xmax = 60, xlab = "Retention time difference (in sec)", ylab = "Cumulative fraction of peptides")
plotErrorCurve(abs(Err_local), "darkgreen", SameGraph = TRUE, xlab = "Retention time difference (in sec)", ylab = "Cumulative fraction of peptides")
plotErrorCurve(abs(Err_hybrid), "red", SameGraph = TRUE, xlab = "Retention time difference (in sec)", ylab = "Cumulative fraction of peptides")
grid(); legend("bottomright",  c("global", "local", "hybrid"), col = c("blue", "darkgreen", "red"), lwd = 1.5, box.lwd = 1)

## ---- fig.width=6, fig.align='center', fig.height=6, echo=FALSE----------
run_pair <- c("run1", "run2")
RUN1 <- oswOutStrep[[run_pair[1]]]; RUN2 <- oswOutStrep[[run_pair[2]]]
cmp <- intersect(RUN1[,1], RUN2[,1]) # First column corresponds to transition_group_record
RUN1 <- RUN1[which(RUN1[,1] %in% cmp), ]
RUN2 <- RUN2[which(RUN2[,1] %in% cmp), ]
RUN1 <- RUN1[match(cmp, RUN1[,1]),]
RUN2 <- RUN2[match(cmp, RUN2[,1]),]
RUNS_RT <- data.frame( "transition_group_record" = RUN1[,1], "RUN1" = RUN1$RT, "RUN2" = RUN2$RT)
RUNS_RT <- RUNS_RT[order(RUNS_RT$RUN1), ]
testPeptides <-intersect(cmp, peptides)
plot(RUNS_RT[,"RUN1"], RUNS_RT[,"RUN2"], xlim = c(2820,2980), ylim = c(2820,2980), xlab = "Run1 retention time (sec)", ylab = "Run2 retention time (sec)", pch=4, cex = 0.4, main = "Non-linear global alignment")
# For testing we want to avoid validation peptides getting used in the fit.
Loess.fit <- loess(RUN2 ~ RUN1, data = RUNS_RT, subset = !transition_group_record %in% testPeptides,
                   span = 0.15, control=loess.control(surface="direct"))
lines(RUNS_RT[,"RUN1"], predict(Loess.fit, newdata=RUNS_RT[,"RUN1"]), col = "blue")
points(StrepAnnot$run1[c(18,19)], StrepAnnot$run2[c(18,19)], col = "red", pch = 19)
lines(RUNS_RT[,"RUN1"], predict(Loess.fit, newdata=RUNS_RT[,"RUN1"])-2*Loess.fit$s, lty = "dashed")
lines(RUNS_RT[,"RUN1"], predict(Loess.fit, newdata=RUNS_RT[,"RUN1"])+2*Loess.fit$s, lty = "dashed")
grid(); legend("topleft", c("LOESS fit", "2*RSE"), lty= c("solid", "dashed"), col = c("blue", "black"))

## ---- fig.width=6, fig.align='center', fig.height=4, fig.show='hold'-----
plotChromatogram(StrepChroms, "run1", "15605_YFMPVHGEYR/3", StrepAnnot)
plotChromatogram(StrepChroms, "run2", "15605_YFMPVHGEYR/3", StrepAnnot)

## ---- echo=FALSE---------------------------------------------------------
peptide <- "15605_YFMPVHGEYR/3"
simMeasure <- "dotProductMasked"
intensityListA <- lapply(StrepChroms[["run1"]][[peptide]], `[[`, 2) # Extracting intensity values
intensityListB <- lapply(StrepChroms[["run2"]][[peptide]], `[[`, 2) # Extracting intensity values
tA <- StrepChroms[["run1"]][[peptide]][[1]][["time"]] # Extracting time component
tB <- StrepChroms[["run2"]][[peptide]][[1]][["time"]] # Extracting time component
localObj <- alignChromatogramsCpp(intensityListA, intensityListB, "local", tA, tB, "mean", simMeasure)

## ----localAlign, fig.width=6, fig.align='center', fig.height=4-----------
s <- localObj@s
Path <- localObj@path[2:nrow(localObj@path), 2:ncol(localObj@path)]
levelplot(s, axes = TRUE, xlab = "run1 index", ylab = "run2 index")
levelplot(s, axes = TRUE, xlab = "run1 index", ylab = "run2 index", main = paste0("Local alignment through the similarity matrix\n for ", peptide), fontsize = 7) + latticeExtra::as.layer(levelplot(Path, col.regions = c("transparent", "green"), alpha = 1, axes = FALSE))

## ---- fig.width=6, fig.align='center', fig.height=4, fig.show='hold'-----
pTR <- plotSingleAlignedChrom(StrepChroms, "run1", peptide, localObj@indexA_aligned, tA, TRUE) + geom_vline(xintercept=which.min(abs(tA - StrepAnnot[peptide, "run1"])), lty="dotted", size = 0.4)
pBR <- plotSingleAlignedChrom(StrepChroms, "run2", peptide, localObj@indexB_aligned, tB, TRUE) + geom_vline(xintercept=which.min(abs(tB - StrepAnnot[peptide, "run2"])), lty="dotted", size = 0.4)
pBR <- pBR + geom_vline(xintercept=which.min(abs(tA - StrepAnnot[peptide, "run1"])), lty="dashed", size = 0.4, color = "red")
pTR
pBR

## ----hybridAlign, fig.width=6, fig.align='center', fig.height=4----------
hybridObj <- getAlignedObj(peptide, pair, StrepChroms[["run1"]], StrepChroms[["run2"]], "hybrid", oswOutStrep, samplingTime = 3.4)
s <- hybridObj@s
Path <- hybridObj@path[2:nrow(hybridObj@path), 2:ncol(hybridObj@path)]
levelplot(s, axes = TRUE, xlab = "run1 index", ylab = "run2 index")
levelplot(s, axes = TRUE, xlab = "run1 index", ylab = "run2 index", main = paste0("Hybrid alignment through the similarity matrix\n for ", peptide), fontsize = 7) + latticeExtra::as.layer(levelplot(Path, col.regions = c("transparent", "green"), alpha = 1, axes = FALSE))

## ----sessionInfo, eval=TRUE----------------------------------------------
devtools::session_info()

