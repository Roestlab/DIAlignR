## ----setup, include = FALSE------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----installDIAlignR, eval=FALSE-------------------------------------------
#  if(!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("DIAlignR")

## ----loadDIAlignR----------------------------------------------------------
library(DIAlignR)

## ----getDataPath-----------------------------------------------------------
dataPath <- system.file("extdata", package = "DIAlignR")

## ---- results=FALSE, message=FALSE, warning=FALSE--------------------------
runs <- c("hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
          "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
# For specific runs provide their names.
intensityTbl <- alignTargetedRuns(dataPath = dataPath, runs = runs, analyteInGroupLabel = TRUE)
# For specific analytes provide their names.
intensityTbl <- alignTargetedRuns(dataPath = dataPath, runs = NULL,
                                  analytes = c("QFNNTDIVLLEDFQK_3"), analyteInGroupLabel = FALSE)
# For all the analytes in all runs, keep them as NULL.
intensityTbl <- alignTargetedRuns(dataPath = dataPath, runs = NULL, analytes = NULL,
                                  analyteInGroupLabel = TRUE)

## ---- message=FALSE--------------------------------------------------------
runs <- c("hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
          "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
AlignObjLight <- getAlignObjs(analytes = "QFNNTDIVLLEDFQK_3", runs = runs, dataPath = dataPath, objType	= "light")
obj <- AlignObjLight[["QFNNTDIVLLEDFQK_3"]][[1]]
slotNames(obj)
names(as.list(obj))
AlignObjMedium <- getAlignObjs(analytes = "QFNNTDIVLLEDFQK_3", runs = runs, dataPath = dataPath, objType	= "medium")
slotNames(AlignObjMedium[["QFNNTDIVLLEDFQK_3"]][[1]])

## ---- fig.width=6, fig.align='center', fig.height=6, message=FALSE---------
runs <- c("hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
 "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
AlignObj <- getAlignObjs(analytes = "QFNNTDIVLLEDFQK_3", runs = runs, dataPath = dataPath)
plotAlignedAnalytes(AlignObj, annotatePeak = TRUE)

## ---- fig.width=5, fig.align='center', fig.height=5, message=FALSE---------
library(lattice)
runs <- c("hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
 "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
AlignObjOutput <- getAlignObjs(analytes = "QFNNTDIVLLEDFQK_3", runs = runs,
                               dataPath = dataPath, objType = "medium")
plotAlignmentPath(AlignObjOutput)

## --------------------------------------------------------------------------
sessionInfo()

