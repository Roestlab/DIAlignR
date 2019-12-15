#' Outputs intensities for each analyte from aligned Targeted-MS runs
#'
#' This function expects osw and mzml directories at dataPath. It first reads osw files and fetches chromatogram indices for each analyte.
#' It then align XICs of each analyte to its reference XICs. Best peak, which has lowest m-score, about the aligned retention time is picked for quantification.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + MIT
#' Date: 2019-12-14
#' @param dataPath (char) Path to mzml and osw directory.
#' @param alignType Available alignment methods are "global", "local" and "hybrid".
#' @param analyteInGroupLabel (logical) TRUE for getting analytes as PRECURSOR.GROUP_LABEL from osw file.
#' @param oswMerged (logical) TRUE for experiment-wide FDR and FALSE for run-specific FDR by pyprophet.
#' @param runs (A vector of string) Names of mzml file without extension.
#' @param analytes (vector of strings) transition_group_ids for which features are to be extracted. analyteInGroupLabel must be set according the pattern used here.
#' @param nameCutPattern (string) regex expression to fetch mzML file name from RUN.FILENAME columns of osw files.
#' @param maxFdrQuery (numeric) A numeric value between 0 and 1. It is used to filter features from osw file which have SCORE_MS2.QVALUE less than itself.
#' @param maxFdrLoess (numeric) A numeric value between 0 and 1. Features should have m-score lower than this value for participation in LOESS fit.
#' @param analyteFDR (numeric) only analytes that have m-score less than this, will be included in the output.
#' @param spanvalue (numeric) Spanvalue for LOESS fit. For targeted proteomics 0.1 could be used.
#' @param runType (char) This must be one of the strings "DIA_proteomics", "DIA_Metabolomics".
#' @param normalization (character) Must be selected from "mean", "l2".
#' @param simMeasure (string) Must be selected from dotProduct, cosineAngle,
#' cosine2Angle, dotProductMasked, euclideanDist, covariance and correlation.
#' @param XICfilter (string) This must be one of the strings "sgolay", "none".
#' @param SgolayFiltOrd (integer) It defines the polynomial order of filer.
#' @param SgolayFiltLen (integer) Must be an odd number. It defines the length of filter.
#' @param goFactor (numeric) Penalty for introducing first gap in alignment. This value is multiplied by base gap-penalty.
#' @param geFactor (numeric) Penalty for introducing subsequent gaps in alignment. This value is multiplied by base gap-penalty.
#' @param cosAngleThresh (numeric) In simType = dotProductMasked mode, angular similarity should be higher than cosAngleThresh otherwise similarity is forced to zero.
#' @param OverlapAlignment (logical) An input for alignment with free end-gaps. False: Global alignment, True: overlap alignment.
#' @param dotProdThresh (numeric) In simType = dotProductMasked mode, values in similarity matrix higher than dotProdThresh quantile are checked for angular similarity.
#' @param gapQuantile (numeric) Must be between 0 and 1. This is used to calculate base gap-penalty from similarity distribution.
#' @param hardConstrain (logical) If FALSE; indices farther from noBeef distance are filled with distance from linear fit line.
#' @param samples4gradient (numeric) This parameter modulates penalization of masked indices.
#' @param samplingTime (numeric) Time difference between two data-points in each chromatogram. For hybrid and local alignment, samples are assumed to be equally time-spaced.
#' @param RSEdistFactor (numeric) This defines how much distance in the unit of rse remains a noBeef zone.
#' @param saveFiles (logical) Must be selected from light, medium and heavy.
#' @return Two tables of intensity and rention times for every analyte in each run.
#' @seealso \code{\link{getRunNames}, \link{getOswFiles}, \link{getAnalytesName}, \link{getMappedRT}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' runs <- c("hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt", "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
#' intensityTbl <- alignTargetedRuns(dataPath, runs = runs, analytes = c("QFNNTDIVLLEDFQK_3"), analyteInGroupLabel = FALSE)
#' intensityTbl <- alignTargetedRuns(dataPath, runs = runs, analytes = c("14299_QFNNTDIVLLEDFQK/3"), analyteInGroupLabel = TRUE)
#' @importFrom dplyr %>%
#' @references Gupta S, Ahadi S, Zhou W, Röst H. "DIAlignR Provides Precise Retention Time Alignment Across Distant Runs in DIA and Targeted Proteomics." Mol Cell Proteomics. 2019 Apr;18(4):806-817. doi: https://doi.org/10.1074/mcp.TIR118.001132 Epub 2019 Jan 31.
#'
#' @export
alignTargetedRuns <- function(dataPath, alignType = "hybrid", analyteInGroupLabel = FALSE, oswMerged = TRUE,
                              runs = NULL, analytes = NULL, nameCutPattern = "(.*)(/)(.*)",
                         maxFdrQuery = 0.05, maxFdrLoess = 0.01, analyteFDR = 0.01,
                         spanvalue = 0.1, runType = "DIA_Proteomics",
                         normalization = "mean", simMeasure = "dotProductMasked",
                         XICfilter = "sgolay", SgolayFiltOrd = 4, SgolayFiltLen = 9,
                         goFactor = 0.125, geFactor = 40,
                         cosAngleThresh = 0.3, OverlapAlignment = TRUE,
                         dotProdThresh = 0.96, gapQuantile = 0.5,
                         hardConstrain = FALSE, samples4gradient = 100,
                         samplingTime = 3.4,  RSEdistFactor = 3.5, saveFiles = FALSE){
  # Check if filter length is odd for Savitzky-Golay filter.
  if( (SgolayFiltLen %% 2) != 1){
    return(stop("SgolayFiltLen can only be odd number"))
  }

  # Get filenames from .merged.osw file and check if names are consistent between osw and mzML files.
  filenames <- getRunNames(dataPath, oswMerged, nameCutPattern)
  if(!is.null(runs)){
    filenames <- filenames[filenames$runs %in% runs,]
    missingRun <- setdiff(runs, filenames$runs)
    if(length(missingRun) != 0){
      return(stop(missingRun, " runs are not found."))
    }
  }
  message("Following runs will be aligned:")
  print(filenames[, "runs"], sep = "\n")

  ######### Get Precursors from the query and respectve chromatogram indices. ######
  oswFiles <- getOswFiles(dataPath, filenames,  maxFdrQuery = maxFdrQuery, analyteFDR = analyteFDR,
                          oswMerged = oswMerged, analytes = NULL, runType = runType, analyteInGroupLabel = analyteInGroupLabel)

  refAnalytes <- getAnalytesName(oswFiles, analyteFDR, commonAnalytes = FALSE)
  if(!is.null(analytes)){
    analytesFound <- intersect(analytes, refAnalytes)
    analytesNotFound <- setdiff(analytes, analytesFound)
    if(length(analytesNotFound)>0){
      message(paste(analytesNotFound, "not found."))
    }
    refAnalytes <- analytesFound
  }

  ######### Collect pointers for each mzML file. #######
  runs <- filenames$runs
  names(runs) <- rownames(filenames)
  # Collect all the pointers for each mzML file.
  message("Collecting metadata from mzML files.")
  mzPntrs <- getMZMLpointers(dataPath, runs)
  message("Metadata is collected from mzML files.")

  ######### Initilize output tables. #######
  rtTbl <- matrix(NA, nrow = length(refAnalytes), ncol = length(runs))
  intesityTbl <- matrix(NA, nrow = length(refAnalytes), ncol = length(runs))
  lwTbl <- matrix(NA, nrow = length(refAnalytes), ncol = length(runs))
  rwTbl <- matrix(NA, nrow = length(refAnalytes), ncol = length(runs))
  rownames(rtTbl) <- refAnalytes; colnames(rtTbl) <- names(runs)
  rownames(intesityTbl) <- refAnalytes; colnames(intesityTbl) <- names(runs)
  rownames(lwTbl) <- refAnalytes; colnames(lwTbl) <- names(runs)
  rownames(rwTbl) <- refAnalytes; colnames(rwTbl) <- names(runs)

  ######### Container to save loess fits.  #######
  loessFits <- list()
  #alignedTables <- performRefAlignment(alignType, ...)

  message("Performing reference-based alignment.")
  start_time <- Sys.time()
  for(analyteIdx in seq_along(refAnalytes)){
    analyte <- refAnalytes[analyteIdx]
    # Select reference run based on m-score
    refRunIdx <- getRefRun(oswFiles, analyte)
    refPeak <- oswFiles[[refRunIdx]] %>%
      dplyr::filter(transition_group_id == analyte & peak_group_rank == 1) %>%
      dplyr::select(leftWidth, RT, rightWidth, Intensity)

    # Get XIC_group from reference run. if missing, go to next analyte.
    ref <- names(runs)[refRunIdx]
    exps <- setdiff(names(runs), ref)
    chromIndices <- selectChromIndices(oswFiles, runname = ref, analyte = analyte)
    if(is.null(chromIndices)){
      warning("Chromatogram indices for ", analyte, " are missing in ", runs[ref])
      message("Skipping ", analyte)
      next
    } else {
      XICs.ref <- extractXIC_group(mz = mzPntrs[[ref]], chromIndices = chromIndices,
                                   XICfilter = XICfilter, SgolayFiltOrd = SgolayFiltOrd,
                                   SgolayFiltLen = SgolayFiltLen)
    }

    # Align all runs to reference run
    for(eXp in exps){
      # Get XIC_group from experiment run
      chromIndices <- selectChromIndices(oswFiles, runname = eXp, analyte = analyte)
      if(!is.null(chromIndices)){
        XICs.eXp <- extractXIC_group(mzPntrs[[eXp]], chromIndices)
        # Get the loess fit for hybrid alignment
        pair <- paste(ref, eXp, sep = "_")
        if(any(pair %in% names(loessFits))){
          Loess.fit <- loessFits[[pair]]
        } else{
          Loess.fit <- getGlobalAlignment(oswFiles, ref, eXp, maxFdrLoess, spanvalue, fitType = "loess")
          loessFits[[pair]] <- Loess.fit
        }
        # Set up constraints for penalizing similarity matrix
        adaptiveRT <- RSEdistFactor*Loess.fit$s
        # Get retention time in experiment run mapped to reference run retention time.
        eXpRT <- getMappedRT(refPeak$RT, XICs.ref, XICs.eXp, Loess.fit, alignType, adaptiveRT, samplingTime,
                             normalization, simMeasure, goFactor, geFactor, cosAngleThresh,
                             OverlapAlignment, dotProdThresh, gapQuantile, hardConstrain,
                             samples4gradient)
        eXp_feature <- pickNearestFeature(eXpRT, analyte, oswFiles, runname = eXp,
                                          adaptiveRT = adaptiveRT, featureFDR = 0.05)
        if(!is.null(eXp_feature)){
          # A feature is found. Use this feature for quantification.
          lwTbl[analyteIdx, eXp] <- eXp_feature[["leftWidth"]]
          rtTbl[analyteIdx, eXp] <- eXp_feature[["RT"]]
          rwTbl[analyteIdx, eXp] <- eXp_feature[["rightWidth"]]
          intesityTbl[analyteIdx, eXp] <- eXp_feature[["Intensity"]]
        } else {
          # Feature is not found.}
        }
      } else {
        warning("Chromatogram indices for ", analyte, " are missing in ", runs[eXp])
        next
      }
    }
    # Get the feature from reference run
    lwTbl[analyteIdx, refRunIdx] <- refPeak[["leftWidth"]]
    rtTbl[analyteIdx, refRunIdx] <- refPeak[["RT"]]
    rwTbl[analyteIdx, refRunIdx] <- refPeak[["rightWidth"]]
    intesityTbl[analyteIdx, refRunIdx] <- refPeak[["Intensity"]]
  }
  ######### Cleanup.  #######
  rm(mzPntrs)
  # Report the execution time for hybrid alignment step.
  end_time <- Sys.time()
  message("Execution time for alignment = ", end_time - start_time)

  colnames(rtTbl) <- unname(runs[colnames(rtTbl)])
  colnames(intesityTbl) <- unname(runs[colnames(intesityTbl)])
  if(saveFiles){
    write.table(rtTbl,file = "rtTbl.csv", col.names = NA, sep = ",")
    write.table(intesityTbl,file = "intesityTbl.csv", col.names = NA, sep = ",")
    print("Data matrix is available in the current directory")
    return(1)
  } else {
    return(intesityTbl)
  }
}

#' AlignObj for analytes between a pair of runs
#'
#' This function expects osw and mzml directories at dataPath. It first reads osw files and fetches chromatogram indices for each requested analyte.
#' It then align XICs of each analyte to its reference XICs. AlignObj is returned which contains aligned indices and cumulative score along the alignment path.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + MIT
#' Date: 2019-12-14
#' @param analytes (vector of strings) transition_group_ids for which features are to be extracted. analyteInGroupLabel must be set according the pattern used here.
#' @param runs (A vector of string) Names of mzml file without extension.
#' @param dataPath (char) Path to mzml and osw directory.
#' @param alignType Available alignment methods are "global", "local" and "hybrid".
#' @param runType (char) This must be one of the strings "DIA_proteomics", "DIA_Metabolomics".
#' @param refRun (string)
#' @param analyteInGroupLabel (logical) TRUE for getting analytes as PRECURSOR.GROUP_LABEL from osw file.
#' @param oswMerged (logical) TRUE for experiment-wide FDR and FALSE for run-specific FDR by pyprophet.
#' @param nameCutPattern (string) regex expression to fetch mzML file name from RUN.FILENAME columns of osw files.
#' @param maxFdrQuery (numeric) A numeric value between 0 and 1. It is used to filter features from osw file which have SCORE_MS2.QVALUE less than itself.
#' @param maxFdrLoess (numeric) A numeric value between 0 and 1. Features should have m-score lower than this value for participation in LOESS fit.
#' @param analyteFDR (numeric) only analytes that have m-score less than this, will be included in the output.
#' @param spanvalue (numeric) Spanvalue for LOESS fit. For targeted proteomics 0.1 could be used.
#' @param normalization (character) Must be selected from "mean", "l2".
#' @param simMeasure (string) Must be selected from dotProduct, cosineAngle,
#' cosine2Angle, dotProductMasked, euclideanDist, covariance and correlation.
#' @param XICfilter (string) This must be one of the strings "sgolay", "none".
#' @param SgolayFiltOrd (integer) It defines the polynomial order of filer.
#' @param SgolayFiltLen (integer) Must be an odd number. It defines the length of filter.
#' @param goFactor (numeric) Penalty for introducing first gap in alignment. This value is multiplied by base gap-penalty.
#' @param geFactor (numeric) Penalty for introducing subsequent gaps in alignment. This value is multiplied by base gap-penalty.
#' @param cosAngleThresh (numeric) In simType = dotProductMasked mode, angular similarity should be higher than cosAngleThresh otherwise similarity is forced to zero.
#' @param OverlapAlignment (logical) An input for alignment with free end-gaps. False: Global alignment, True: overlap alignment.
#' @param dotProdThresh (numeric) In simType = dotProductMasked mode, values in similarity matrix higher than dotProdThresh quantile are checked for angular similarity.
#' @param gapQuantile (numeric) Must be between 0 and 1. This is used to calculate base gap-penalty from similarity distribution.
#' @param hardConstrain (logical) If FALSE; indices farther from noBeef distance are filled with distance from linear fit line.
#' @param samples4gradient (numeric) This parameter modulates penalization of masked indices.
#' @param samplingTime (numeric) Time difference between two data-points in each chromatogram. For hybrid and local alignment, samples are assumed to be equally time-spaced.
#' @param RSEdistFactor (numeric) This defines how much distance in the unit of rse remains a noBeef zone.
#' @param objType (char) Must be selected from light, medium and heavy.
#' @return A list of AlignObj. Each AlignObj is an S4 object. Three most-important slots are:
#' \item{indexA_aligned}{(integer) aligned indices of reference run.}
#' \item{indexB_aligned}{(integer) aligned indices of experiment run.}
#' \item{score}{(numeric) cumulative score of alignment.}
#' @seealso \code{\link{plotAlignedAnalytes}, \link{getRunNames}, \link{getOswFiles}, \link{getXICs4AlignObj}, \link{getAlignObj}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' runs <- c("hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt", "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
#' AlignObjOutput <- getAlignObjs(analytes = "QFNNTDIVLLEDFQK_3", runs, dataPath = dataPath)
#' plotAlignedAnalytes(AlignObjOutput)
#'
#' @references Gupta S, Ahadi S, Zhou W, Röst H. "DIAlignR Provides Precise Retention Time Alignment Across Distant Runs in DIA and Targeted Proteomics." Mol Cell Proteomics. 2019 Apr;18(4):806-817. doi: https://doi.org/10.1074/mcp.TIR118.001132 Epub 2019 Jan 31.
#'
#' @export
getAlignObjs <- function(analytes, runs, dataPath = ".", alignType = "hybrid",
                         runType = "DIA_Proteomics", refRun = NULL,
                         analyteInGroupLabel = FALSE, oswMerged = TRUE, nameCutPattern = "(.*)(/)(.*)",
                         maxFdrQuery = 0.05, maxFdrLoess = 0.01, analyteFDR = 1.00, spanvalue = 0.1,
                         normalization = "mean", simMeasure = "dotProductMasked",
                         XICfilter = "sgolay", SgolayFiltOrd = 4, SgolayFiltLen = 9,
                         goFactor = 0.125, geFactor = 40,
                         cosAngleThresh = 0.3, OverlapAlignment = TRUE,
                         dotProdThresh = 0.96, gapQuantile = 0.5,
                         hardConstrain = FALSE, samples4gradient = 100,
                         samplingTime = 3.4,  RSEdistFactor = 3.5, objType = "light"){
  if(length(runs) != 2){
    print("For pairwise alignment, two runs are required.")
    return(NULL)
  }

  if( (SgolayFiltLen %% 2) != 1){
    print("SgolayFiltLen can only be odd number")
    return(NULL)
  }
  ##### Get filenames from osw files and check if names are consistent between osw and mzML files. ######
  filenames <- getRunNames(dataPath, oswMerged, nameCutPattern)
  filenames <- filenames[filenames$runs %in% runs,]
  missingRun <- setdiff(runs, filenames$runs)
  if(length(missingRun) != 0){
    return(stop(missingRun, " runs are not found."))
  }

  message("Following runs will be aligned:")
  print(filenames[, "runs"], sep = "\n")

  ######### Get Precursors from the query and respectve chromatogram indices. ######
  oswFiles <- getOswFiles(dataPath, filenames, maxFdrQuery = maxFdrQuery, analyteFDR = analyteFDR,
                          oswMerged = oswMerged, analytes = NULL, runType = runType, analyteInGroupLabel = analyteInGroupLabel)

  # Report analytes that are not found
  refAnalytes <- getAnalytesName(oswFiles, analyteFDR, commonAnalytes = FALSE)
  analytesFound <- intersect(analytes, refAnalytes)
  analytesNotFound <- setdiff(analytes, analytesFound)
  if(length(analytesNotFound)>0){
    message(paste(analytesNotFound, "not found."))
  }
  analytes <- analytesFound

  ####################### Get XICs ##########################################
  runs <- filenames$runs
  names(runs) <- rownames(filenames)
  # Get Chromatogram for each peptide in each run.
  message("Fetching Extracted-ion chromatograms from runs")
  XICs <- getXICs4AlignObj(dataPath, runs, oswFiles, analytes, XICfilter = XICfilter,
                           SgolayFiltOrd = SgolayFiltOrd, SgolayFiltLen = SgolayFiltLen)

  ####################### Perfrom alignment ##########################################
  AlignObjs <- vector("list", length(analytes))
  names(AlignObjs) <- analytes
  loessFits <- list()
  print("Perfroming alignment")
  for(analyteIdx in seq_along(analytes)){
    analyte <- analytes[analyteIdx]
    # Select reference run based on m-score
    if(is.null(refRun)){
      refRunIdx <- getRefRun(oswFiles, analyte)
    } else{
      refRunIdx <- which(filenames$runs == refRun)
    }

    # Get XIC_group from reference run
    ref <- names(runs)[refRunIdx]
    exps <- setdiff(names(runs), ref)
    XICs.ref <- XICs[[ref]][[analyte]]

    # Align experiment run to reference run
    for(eXp in exps){
      # Get XIC_group from experiment run
      XICs.eXp <- XICs[[eXp]][[analyte]]
      if(!is.null(XICs.eXp)){
        # Get the loess fit for hybrid alignment
        pair <- paste(ref, eXp, sep = "_")
        if(any(pair %in% names(loessFits))){
          Loess.fit <- loessFits[[pair]]
        } else{
          Loess.fit <- getGlobalAlignment(oswFiles, ref, eXp, maxFdrLoess, spanvalue, fitType = "loess")
          loessFits[[pair]] <- Loess.fit
        }
        adaptiveRT <-  RSEdistFactor*Loess.fit$s # Residual Standard Error
        # Fetch alignment object between XICs.ref and XICs.eXp
        AlignObj <- getAlignObj(XICs.ref, XICs.eXp, Loess.fit, adaptiveRT = adaptiveRT, samplingTime,
                                normalization, simType = simMeasure, goFactor, geFactor,
                                cosAngleThresh, OverlapAlignment,
                                dotProdThresh, gapQuantile, hardConstrain, samples4gradient,
                                objType)
        AlignObjs[[analyte]] <- list()
        # Attach AlignObj for the analyte.
          AlignObjs[[analyte]][[pair]] <- AlignObj
        # Attach intensities of reference XICs.
        AlignObjs[[analyte]][[runs[ref]]] <- XICs.ref
        # Attach intensities of experiment XICs.
        AlignObjs[[analyte]][[runs[eXp]]] <- XICs.eXp
        # Attach peak boundaries to the object.
        AlignObjs[[analyte]][[paste0(pair, "_pk")]] <- oswFiles[[refRunIdx]] %>%
          dplyr::filter(transition_group_id == analyte & peak_group_rank == 1) %>%
          dplyr::select(leftWidth, RT, rightWidth) %>%
          as.vector()
      }
      else {AlignObjs[[analyte]] <- NULL}
    }
  }

  ####################### Return AlignedObjs ##########################################
  message("Alignment done. Returning AlignObjs")
  AlignObjs
}
