#' Extract XICs of all transitions requested in chromIndices.
#'
#' Extracts XICs using mz object. Generally Savitzky–Golay filter is used, however, filter can be turned-off as well.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + MIT
#' Date: 2019-12-13
#' @param mz (mzRpwiz object)
#' @param chromIndices (vector of Integers) Indices of chromatograms to be extracted.
#' @param XICfilter (string) This must be one of the strings "sgolay", "none".
#' @param SgolayFiltOrd (integer) It defines the polynomial order of filer.
#' @param SgolayFiltLen (integer) Must be an odd number. It defines the length of filter.
#' @return A list of data-frames. Each data frame has elution time and intensity of fragment-ion XIC.
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' mzmlName <- paste0(dataPath,"/mzml/hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML")
#' mz <- mzR::openMSfile(mzmlName, backend = "pwiz")
#' chromIndices <- c(37L, 38L, 39L, 40L, 41L, 42L)
#' XIC_group <- extractXIC_group(mz, chromIndices, SgolayFiltOrd = 4, SgolayFiltLen = 13)
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

#' Extract XICs of all analytes from oswFiles
#'
#' For all the analytes requested, it fetches chromatogram indices from oswFiles and
#' extract chromatograms from mzML files.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + MIT
#' Date: 2019-12-13
#' @param dataPath (char) path to mzml and osw directory.
#' @param runs (vector of string) names of mzML files without extension. Names of the vector must be a combination of "run" and an iteger e.g. "run2".
#' @param oswFiles (list of data-frames) it is output from getOswFiles function.
#' @param analytes (string) analyte is as PRECURSOR.GROUP_LABEL or as PEPTIDE.MODIFIED_SEQUENCE and PRECURSOR.CHARGE from osw file.
#' @param XICfilter (string) this must be one of the strings "sgolay", "none".
#' @param SgolayFiltOrd (integer) it defines the polynomial order of filer.
#' @param SgolayFiltLen (integer) must be an odd number. It defines the length of filter.
#' @return A list of list of data-frames. Each data frame has elution time and intensity of fragment-ion XIC.
#'
#' @seealso \code{\link{getOswFiles}, \link{getRunNames}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' filenames <- DIAlignR::getRunNames(dataPath = dataPath)
#' runs <- c("run1" = "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
#'           "run0" =  "hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt")
#' oswFiles <- DIAlignR::getOswFiles(dataPath, filenames)
#' analytes <- "QFNNTDIVLLEDFQK_3"
#' XICs <- getXICs4AlignObj(dataPath, runs, oswFiles, analytes)
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
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + MIT
#' Date: 2019-12-13
#'
#' @param XICs.ref List of extracted ion chromatograms from reference run.
#' @param XICs.eXp List of extracted ion chromatograms from experiment run.
#' @param Loess.fit LOESS fit object between reference and experiment run.
#' @param adaptiveRT (numeric) Similarity matrix is not penalized within adaptive RT.
#' @param samplingTime (numeric) Time difference between two data-points in each chromatogram. For hybrid and local alignment, samples are assumed to be equally time-spaced.
#' @param normalization (character) Must be selected from "mean", "l2".
#' @param simType (string) Must be selected from dotProduct, cosineAngle,
#' cosine2Angle, dotProductMasked, euclideanDist, covariance and correlation.
#' @param goFactor (numeric) Penalty for introducing first gap in alignment. This value is multiplied by base gap-penalty.
#' @param geFactor (numeric) Penalty for introducing subsequent gaps in alignment. This value is multiplied by base gap-penalty.
#' @param cosAngleThresh (numeric) In simType = dotProductMasked mode, angular similarity should be higher than cosAngleThresh otherwise similarity is forced to zero.
#' @param OverlapAlignment (logical) An input for alignment with free end-gaps. False: Global alignment, True: overlap alignment.
#' @param dotProdThresh (numeric) In simType = dotProductMasked mode, values in similarity matrix higher than dotProdThresh quantile are checked for angular similarity.
#' @param gapQuantile (numeric) Must be between 0 and 1. This is used to calculate base gap-penalty from similarity distribution.
#' @param hardConstrain (logical) If FALSE; indices farther from noBeef distance are filled with distance from linear fit line.
#' @param samples4gradient (numeric) This parameter modulates penalization of masked indices.
#' @param objType (char) Must be selected from light, medium and heavy.
#' @return A S4 object. Three most-important slots are:
#' \item{indexA_aligned}{(integer) aligned indices of reference run.}
#' \item{indexB_aligned}{(integer) aligned indices of experiment run.}
#' \item{score}{(numeric) cumulative score of alignment.}
#'
#' @seealso \code{\link{alignChromatogramsCpp}}
#' @examples
#' data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
#' data(oswFiles_DIAlignR, package="DIAlignR")
#' XICs.ref <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["run1"]][["14299_QFNNTDIVLLEDFQK/3"]]
#' XICs.eXp <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["run2"]][["14299_QFNNTDIVLLEDFQK/3"]]
#' Loess.fit <- getLOESSfit(oswFiles_DIAlignR, ref = "run1", eXp = "run2", maxFdrLoess = 0.05, spanvalue = 0.1)
#' AlignObj <- getAlignObj(XICs.ref, XICs.eXp, Loess.fit, adaptiveRT = 77.82315, samplingTime = 3.414,
#' normalization = "mean", simType = "dotProductMasked", goFactor = 0.125, geFactor = 40,
#' cosAngleThresh = 0.3, OverlapAlignment = TRUE, dotProdThresh = 0.96, gapQuantile = 0.5, hardConstrain = FALSE,
#' samples4gradient = 100, objType = "light")
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

#' Get mapping of reference RT on experiment run.
#'
#' This function aligns XICs of reference and experiment runs. Using alignment, it maps retention time from refernce run on experiment run.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + MIT
#' Date: 2019-12-13
#' @param XICs.ref List of extracted ion chromatograms from reference run.
#' @param XICs.eXp List of extracted ion chromatograms from experiment run.
#' @param Loess.fit LOESS fit object between reference and experiment run.
#' @param alignType Available alignment methods are "global", "local" and "hybrid".
#' @param adaptiveRT (numeric) Similarity matrix is not penalized within adaptive RT.
#' @param samplingTime (numeric) Time difference between two data-points in each chromatogram. For hybrid and local alignment, samples are assumed to be equally time-spaced.
#' @param normalization (character) Must be selected from "mean", "l2".
#' @param simMeasure (string) Must be selected from dotProduct, cosineAngle,
#' cosine2Angle, dotProductMasked, euclideanDist, covariance and correlation.
#' @param goFactor (numeric) Penalty for introducing first gap in alignment. This value is multiplied by base gap-penalty.
#' @param geFactor (numeric) Penalty for introducing subsequent gaps in alignment. This value is multiplied by base gap-penalty.
#' @param cosAngleThresh (numeric) In simType = dotProductMasked mode, angular similarity should be higher than cosAngleThresh otherwise similarity is forced to zero.
#' @param OverlapAlignment (logical) An input for alignment with free end-gaps. False: Global alignment, True: overlap alignment.
#' @param dotProdThresh (numeric) In simType = dotProductMasked mode, values in similarity matrix higher than dotProdThresh quantile are checked for angular similarity.
#' @param gapQuantile (numeric) Must be between 0 and 1. This is used to calculate base gap-penalty from similarity distribution.
#' @param hardConstrain (logical) If FALSE; indices farther from noBeef distance are filled with distance from linear fit line.
#' @param samples4gradient (numeric) This parameter modulates penalization of masked indices.
#' @param objType (char) Must be selected from light, medium and heavy.
#' @return (numeric)
#' @seealso \code{\link{alignChromatogramsCpp}}
#' @examples
#' data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
#' data(oswFiles_DIAlignR, package="DIAlignR")
#' XICs.ref <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["run1"]][["14299_QFNNTDIVLLEDFQK/3"]]
#' XICs.eXp <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["run2"]][["14299_QFNNTDIVLLEDFQK/3"]]
#' Loess.fit <- getLOESSfit(oswFiles_DIAlignR, ref = "run2", eXp = "run0", maxFdrLoess = 0.05, spanvalue = 0.1)
#' adaptiveRT <- 77.82315 #3.5*Loess.fit$s
#' getMappedRT(refRT = 5238.35, XICs.ref, XICs.eXp, Loess.fit, alignType = "hybrid",
#'  adaptiveRT = adaptiveRT, samplingTime = 3.414, normalization = "mean", simMeasure = "dotProductMasked",
#'  goFactor = 0.125, geFactor = 40, cosAngleThresh = 0.3, OverlapAlignment = TRUE,
#'  dotProdThresh = 0.96, gapQuantile = 0.5, hardConstrain = FALSE, samples4gradient = 100)
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


#' Get XICs of all analytes
#'
#' For all the analytes requested in runs, it first creates oswFiles, then, fetches chromatogram indices from oswFiles and
#' extract chromatograms from mzML files.
#'
#' @importFrom dplyr %>%
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + MIT
#' Date: 2019-12-13
#'
#' @param analytes (string) An analyte is as PRECURSOR.GROUP_LABEL or as PEPTIDE.MODIFIED_SEQUENCE and PRECURSOR.CHARGE from osw file.
#' @param dataPath (char) Path to mzml and osw directory.
#' @param runs (A vector of string) Names of mzml file without extension. Vector must have names as shown in the example.
#' @param maxFdrQuery (numeric) A numeric value between 0 and 1. It is used to filter features from osw file which have SCORE_MS2.QVALUE less than itself.
#' @param XICfilter (string) This must be one of the strings "sgolay", "none".
#' @param SgolayFiltOrd (integer) It defines the polynomial order of filer.
#' @param SgolayFiltLen (integer) Must be an odd number. It defines the length of filter.
#' @param runType (char) This must be one of the strings "DIA_proteomics", "DIA_Metabolomics".
#' @param oswMerged (logical) TRUE for experiment-wide FDR and FALSE for run-specific FDR by pyprophet.
#' @param nameCutPattern (string) regex expression to fetch mzML file name from RUN.FILENAME columns of osw files.
#' @param analyteInGroupLabel (logical) TRUE for getting analytes as PRECURSOR.GROUP_LABEL from osw file.
#' @return A list of list. Each list contains XIC-group for that run. XIC-group is a list of dataframe that has elution time and intensity of fragment-ion XIC.
#'
#' @seealso \code{\link{getOswFiles}, \link{getRunNames}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' runs <- c("hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt", "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
#' XICs <- getXICs(analytes = c("QFNNTDIVLLEDFQK_3"), runs = runs, dataPath = dataPath)
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
  names(XICs) <- filenames$runs
  XICs
}
