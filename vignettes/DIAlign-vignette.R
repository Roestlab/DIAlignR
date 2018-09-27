## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----installDIAlign------------------------------------------------------
require(devtools)
install_github("Roestlab/DIAlign")
library(DIAlign)

## ----globalFit-----------------------------------------------------------
run_pair <- c("run1", "run2")
loess.fit <- getLOESSfit(run_pair, peptides, oswOutStrep, 0.15)
StrepAnnot <- as.data.frame(StrepAnnot)
predict.run2 <- predict(loess.fit, data.frame(RUN1 = StrepAnnot[, run_pair[1]]))
Err <- predict.run2 - StrepAnnot[,run_pair[2]]

## ----plotGlobal----------------------------------------------------------
plotErrorCurve <- function(x, clr = "black", SameGraph = FALSE, xmax = 120, ...){
    x <- x[!is.na(x)]
    breaks = seq(0, xmax, by=0.5)
    duration.cut = cut(x, breaks, right = FALSE) 
    duration.freq = table(duration.cut)
    cumfreq0 = c(0, cumsum(duration.freq))
    if(SameGraph == TRUE){lines(breaks, cumfreq0/length(x), col = clr, ...)}
    else{plot(breaks, cumfreq0/length(x), col = clr, type = "l", ...)}
}
plotErrorCurve(abs(Err), "blue", xlab = "Retention time difference (in sec)", ylab = "Cumulative fraction of peptides")

## ----localFit------------------------------------------------------------
gapQuantile <- 0.5; goFactor <- 1/8; geFactor <- 40
simMeasure <- "dotProductMasked"
run_pair <- c("run1", "run2")
Err <- matrix(NA, nrow = length(peptides), ncol = 1)
rownames(Err) <- peptides
for(peptide in peptides){
  s <- getSimilarityMatrix(StrepChroms, peptide, run_pair[1], run_pair[2], type = simMeasure)
  gapPenalty <- getGapPenalty(s, gapQuantile, type = simMeasure)
  Alignobj <- getAffineAlignObj(s, go = gapPenalty*goFactor, ge = gapPenalty*geFactor)
  AlignedIndices <- getAlignment(Alignobj)
  tA <- StrepChroms[[run_pair[1]]][[peptide]][[1]][["time"]]
  tB <- StrepChroms[[run_pair[2]]][[peptide]][[1]][["time"]]
  tA.aligned <- mapIdxToTime(tA, AlignedIndices[[1]][,"indexA_aligned"])
  tB.aligned <- mapIdxToTime(tB, AlignedIndices[[1]][,"indexB_aligned"])
  predictTime <- tB.aligned[which.min(abs(tA.aligned - StrepAnnot[peptide, run_pair[1]]))]
  deltaT <- predictTime - StrepAnnot[peptide, run_pair[2]]
  Err[peptide, 1] <- deltaT
}

## ----plotLocal-----------------------------------------------------------
plotErrorCurve(abs(Err), "red", xlab = "Retention time difference (in sec)", ylab = "Cumulative fraction of peptides")

## ------------------------------------------------------------------------
library(lattice)
levelplot(s, axes = TRUE, xlab = "run1 index", ylab = "run2 index")
Path <- getAlignmentPath(AlignedIndices[[1]], s)
levelplot(s, axes = TRUE, xlab = "run1 index", ylab = "run2 index", main = paste0("Alignment path through the ", simMeasure, " similarity matrix\n for ", peptide)) + latticeExtra::as.layer(levelplot(Path, col.regions = c("transparent", "green"), alpha = 1, axes = FALSE))

## ------------------------------------------------------------------------
library(ggplot2)
library(reshape2)
plotChromatogram <- function(data, run, peptide, StrepAnnot, printTitle =TRUE){
  df <- do.call("cbind", data[[run]][[peptide]])
  df <- df[,!duplicated(colnames(df))]
  df <- melt(df, id.vars="time", value.name = "Intensity")
  g <- ggplot(df, aes(time, Intensity, col=variable)) + geom_line(show.legend = FALSE) + theme_bw()
  if(printTitle) g <- g + ggtitle(paste0(run, ", ",peptide)) + theme(plot.title = element_text(hjust = 0.5))
  g <- g + geom_vline(xintercept=StrepAnnot[peptide, run], lty="dotted", size = 0.4)
  return(g)
}

## ---- fig.show='hold'----------------------------------------------------
plot(1:10)
plot(10:1)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(mtcars, 10))

