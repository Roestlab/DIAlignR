#' Extract XICs of all transitions requested in chromIndices.
#'
#' @return A list of data-frames. Each data frame has elution time and intensity of fragment-ion XIC.
extractXIC_group <- function(mz, chromIndices, XICfilter = "sgolay", SgolayFiltOrd = 4, SgolayFiltLen = 9){
  XIC_group <- lapply(1:length(chromIndices), function(i) {
    rawChrom <- mzR::chromatograms(mz, chromIndices[i])
    # Savitzky-Golay filter to smooth chromatograms, filter order p = 3, filter length n = 13
    if(XICfilter == "sgolay"){
      rawChrom[,2] <- signal::sgolayfilt(rawChrom[,2], p = SgolayFiltOrd, n = SgolayFiltLen)
    }
    return(rawChrom)
  })
  return(XIC_group)
}

#' Extract XICs of all transitions requested in chromIndices.
#'
#' @return A list of list of data-frames. Each data frame has elution time and intensity of fragment-ion XIC.
getXICs4AlignObj <- function(dataPath, runs, oswFiles, analytes, XICfilter = "sgolay",
                             SgolayFiltOrd = 4, SgolayFiltLen = 9){
  mzPntrs <- getMZMLpointers(dataPath, runs)
  XICs <- vector("list", length(runs))
  names(XICs) <- names(runs)
  for(i in seq_along(runs)){
    runname = names(runs)[i]
    message("Fetching XICs from run ", runs[[runname]])
    XICs[[i]] <- lapply(1:length(analytes), function(j){
      analyte <- analytes[j]
      chromIndices <- selectChromIndices(oswFiles, runname = runname, analyte = analyte)
      if(is.null(chromIndices)){
        warning("Chromatogram indices for ", analyte, " are missing in ", runs[[runname]])
        message("Skipping ", analyte)
        XIC_group <- NULL
      } else {
        XIC_group <- extractXIC_group(mzPntrs[[runname]], chromIndices, XICfilter, SgolayFiltOrd, SgolayFiltLen)
      }
      XIC_group
    })
    names(XICs[[i]]) <- analytes
  }
  rm(mzPntrs)
  XICs
}

#' Extract XICs of all transitions requested in chromIndices.
#'
#' @return A list of list of data-frames. Each data frame has elution time and intensity of fragment-ion XIC.
getAlignObj <- function(XICs.ref, XICs.eXp, Loess.fit, adaptiveRT, samplingTime,
                        normalization, simType, goFactor, geFactor,
                        cosAngleThresh, OverlapAlignment,
                        dotProdThresh, gapQuantile, hardConstrain,
                        samples4gradient, objType = "light"){
  # Set up constraints for penalizing similarity matrix
  noBeef <- ceiling(adaptiveRT/samplingTime)
  tVec.ref <- XICs.ref[[1]][["time"]] # Extracting time component
  tVec.eXp <- XICs.eXp[[1]][["time"]] # Extracting time component
  B1p <- predict(Loess.fit, tVec.ref[1])
  B2p <- predict(Loess.fit, tVec.ref[length(tVec.ref)])
  # Perform dynamic programming for chromatogram alignment
  intensityList.ref <- lapply(XICs.ref, `[[`, 2) # Extracting intensity values
  intensityList.eXp <- lapply(XICs.eXp, `[[`, 2) # Extracting intensity values
  AlignObj <- alignChromatogramsCpp(intensityList.ref, intensityList.eXp,
                                    alignType = "hybrid", tVec.ref, tVec.eXp,
                                    normalization = normalization, simType = simType,
                                    B1p = B1p, B2p = B2p, noBeef = noBeef,
                                    goFactor = goFactor, geFactor = geFactor,
                                    cosAngleThresh = cosAngleThresh, OverlapAlignment = OverlapAlignment,
                                    dotProdThresh = dotProdThresh, gapQuantile = gapQuantile,
                                    hardConstrain = hardConstrain, samples4gradient = samples4gradient,
                                    objType = objType)
  AlignObj
}

#' Extract XICs of all transitions requested in chromIndices.
#'
#' @return A list of list of data-frames. Each data frame has elution time and intensity of fragment-ion XIC.
#' @export
getMappedRT <- function(refRT, XICs.ref, XICs.eXp, Loess.fit, alignType, adaptiveRT, samplingTime,
                        normalization, simMeasure, goFactor, geFactor, cosAngleThresh,
                        OverlapAlignment, dotProdThresh, gapQuantile, hardConstrain,
                        samples4gradient, objType = "light"){
  AlignObj <- getAlignObj(XICs.ref, XICs.eXp, Loess.fit, adaptiveRT, samplingTime,
                          normalization, simType = simMeasure, goFactor, geFactor,
                          cosAngleThresh, OverlapAlignment,
                          dotProdThresh, gapQuantile, hardConstrain, samples4gradient, objType)
  tVec.ref <- XICs.ref[[1]][["time"]] # Extracting time component
  tVec.eXp <- XICs.eXp[[1]][["time"]] # Extracting time component
  eXpRT <- mappedRTfromAlignObj(refRT, tVec.ref, tVec.eXp, AlignObj)
  eXpRT
}


#' Get XICs for a list of peptides
#'
#' @return A list of list. Each list contains XICs for that run.
#' @importFrom dplyr %>%
#' @export
getXICs <- function(analytes, runs, dataPath = ".", maxFdrQuery = 1.0, XICfilter = "sgolay",
                    SgolayFiltOrd = 4, SgolayFiltLen = 9, runType = "DIA_proteomics",
                    oswMerged = TRUE, nameCutPattern = "(.*)(/)(.*)", analyteInGroupLabel = FALSE){
  if( (SgolayFiltLen %% 2) != 1){
    print("SgolayFiltLen can only be odd number")
    return(NULL)
  }
  # Get filenames from .merged.osw file and check if names are consistent between osw and mzML files.
  filenames <- getRunNames(dataPath, oswMerged, nameCutPattern)
  filenames <- filenames[filenames$runs %in% runs,]

  # Get Chromatogram indices for each peptide in each run.
  oswFiles <- getOswFiles(dataPath, filenames, maxFdrQuery = maxFdrQuery, analyteFDR = 1.00,
                         oswMerged = oswMerged, analytes = analytes, runType = runType,
                         analyteInGroupLabel = analyteInGroupLabel)
  refAnalytes <- getAnalytesName(oswFiles, commonAnalytes = FALSE)
  analytesFound <- intersect(analytes, refAnalytes)
  analytesNotFound <- setdiff(analytes, analytesFound)
  if(length(analytesNotFound)>0){
    message("Analytes ", paste(analytesNotFound, ", "), "not found.")
  }

  ####################### Get XICs ##########################################
  runs <- filenames$runs
  names(runs) <- rownames(filenames)
  # Get Chromatogram for each peptide in each run.
  message("Fetching Extracted-ion chromatograms from runs")
  XICs <- getXICs4AlignObj(dataPath, runs, oswFiles, analytesFound, XICfilter,
                           SgolayFiltOrd, SgolayFiltLen)
  XICs
}
