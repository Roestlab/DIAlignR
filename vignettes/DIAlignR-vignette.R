## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----installDIAlignR, eval=FALSE----------------------------------------------
#  if(!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("DIAlignR")

## ----loadDIAlignR-------------------------------------------------------------
library(DIAlignR)

## ----getDataPath--------------------------------------------------------------
dataPath <- system.file("extdata", package = "DIAlignR")

## ---- results=FALSE, message=FALSE, warning=FALSE-----------------------------
runs <- c("hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
          "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
params <- paramsDIAlignR()
params[["context"]] <- "experiment-wide"
# For specific runs provide their names.
alignTargetedRuns(dataPath = dataPath, outFile = "test.csv", runs = runs, oswMerged = TRUE, params = params)
# For all the analytes in all runs, keep them as NULL.
alignTargetedRuns(dataPath = dataPath, outFile = "test.csv", runs = NULL, oswMerged = TRUE, params = params)

## ---- message=FALSE-----------------------------------------------------------
runs <- c("hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
          "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
AlignObjLight <- getAlignObjs(analytes = 4618L, runs = runs, dataPath = dataPath, objType	= "light", params = params)
# First element contains names of runs, spectra files, chromatogram files and feature files.
AlignObjLight[[1]][, c("runName", "spectraFile")]
obj <- AlignObjLight[[2]][["4618"]][[1]][["AlignObj"]]
slotNames(obj)
names(as.list(obj))
AlignObjMedium <- getAlignObjs(analytes = 4618L, runs = runs, dataPath = dataPath, objType	= "medium", params = params)
obj <- AlignObjMedium[[2]][["4618"]][[1]][["AlignObj"]]
slotNames(obj)

## ---- fig.width=6, fig.align='center', fig.height=6, message=FALSE------------
runs <- c("hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
 "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
AlignObj <- getAlignObjs(analytes = 4618L, runs = runs, dataPath = dataPath, params = params)
plotAlignedAnalytes(AlignObj, annotatePeak = TRUE)

## ---- fig.width=5, fig.align='center', fig.height=5, message=FALSE------------
library(lattice)
runs <- c("hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
 "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
AlignObjOutput <- getAlignObjs(analytes = 4618L, runs = runs, params = params, dataPath = dataPath, objType = "medium")
plotAlignmentPath(AlignObjOutput)

## -----------------------------------------------------------------------------
sessionInfo()

