#' Align XICs of precursors across multiple Targeted-MS runs and outputs quantitative data matrix.
#'
#' @examples
#' dataPath <- "data/testData2"
#' analytes <- c("AAAGPVADLF_2", "YQYSDQGIDY_2")
#' runs <- c("170407_AM_BD-ZH12_Spleen_W_10%_DIA_#2_400-650mz_msms42", "170413_AM_BD-ZH12_BoneMarrow_W_10%_DIA_#2_400-650mz_msms35")
#' alignTargetedruns(dataPath, runs = runs, analytes = analytes, oswMerged = FALSE, SgolayFiltOrd = 4, SgolayFiltLen = 5)
#' alignTargetedruns(dataPath, runType = "DIA_Metabolomics", samplingTime = 0.8229)
#' @return Saves intensity table in the current directory.
#' @importFrom dplyr %>%
#' @export
alignTargetedruns <- function(dataPath, alignType = "hybrid", oswMerged = TRUE,
                              runs = NULL, analytes = NULL, nameCutPattern = "(.*)(/)(.*)",
                         maxFdrQuery = 0.05, maxFdrLoess = 0.01, analyteFDR = 0.01,
                         spanvalue = 0.1, runType = "DIA_Proteomics",
                         normalization = "mean", simMeasure = "dotProductMasked",
                         SgolayFiltOrd = 4, SgolayFiltLen = 9,
                         goFactor = 0.125, geFactor = 40,
                         cosAngleThresh = 0.3, OverlapAlignment = TRUE,
                         dotProdThresh = 0.96, gapQuantile = 0.5,
                         hardConstrain = FALSE, samples4gradient = 100,
                         samplingTime = 3.4,  RSEdistFactor = 3.5){
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
  oswFiles <- getOswFiles(dataPath, filenames, maxFdrQuery, analyteFDR,
                          oswMerged, analytes = NULL, runType = runType)

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
      XICs.ref <- extractXIC_group(mzPntrs[[ref]], chromIndices, SgolayFiltOrd, SgolayFiltLen)
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
          Loess.fit <- getLOESSfit(oswFiles, ref, eXp, maxFdrLoess, spanvalue)
          loessFits[[pair]] <- Loess.fit
        }
        # Set up constraints for penalizing similarity matrix
        rse <- Loess.fit$s
        adaptiveRT <- RSEdistFactor*rse
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
  end_time <- Sys.time()
  # Report the execution time for hybrid alignment step.
  print(end_time - start_time)
  colnames(rtTbl) <- runs[colnames(rtTbl)]
  write.table(rtTbl,file = "rtTbl.csv", col.names = NA, sep = ",")
  colnames(intesityTbl) <- runs[colnames(intesityTbl)]
  write.table(intesityTbl,file = "intesityTbl.csv", col.names = NA, sep = ",")

  ######### Cleanup.  #######
  rm(mzPntrs)
  rm(oswFiles, loessFits)
  rm(rtTbl, intesityTbl, lwTbl, rwTbl)
  print("Data matrix is available in the current directory")
  return(1)
}


#' AlignObj for analytes between a pair
#'
#' @return A list of AlignObj. Each AlignObj contains alignment path, similarity matrix and related parameters.
#' @export
getAlignObjs <- function(analytes, runs, dataPath = ".", alignType = "hybrid",
                         runType = "DIA_Proteomics", refRun = NULL,
                         oswMerged = TRUE, nameCutPattern = "(.*)(/)(.*)",
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
  oswFiles <- getOswFiles(dataPath, filenames, maxFdrQuery, analyteFDR,
                          oswMerged, analytes = NULL, runType)

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
          Loess.fit <- getLOESSfit(oswFiles, ref, eXp, maxFdrLoess, spanvalue)
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
