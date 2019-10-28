#' Align XICs of precursors across multiple MRM runs and outputs quantitative data matrix.
#'
#' In MRM runs, since XICs have less interference, higher quantile value from similarity
#' matrix is needed to get base-gap-penalty.
#' @return Saves intensity table in the current directory.
#' @importFrom dplyr %>%
#' @export
alignMRMruns <- function(dataPath, alignType = "local", oswMerged = FALSE, nameCutPattern = "(.*)(/)(.*)",
                              maxFdrQuery = 0.05, maxFdrLoess = 0.01, spanvalue = 0.1, runType = "MRM_Proteomics",
                              normalization = "mean", simMeasure = "dotProductMasked",
                              SgolayFiltOrd = 2, SgolayFiltLen = 3,
                              goFactor = 0.125, geFactor = 400,
                              cosAngleThresh = 0.3, OverlapAlignment = TRUE,
                              dotProdThresh = 0.96, gapQuantile = 0.98,
                              hardConstrain = FALSE, samples4gradient = 100,
                              expRSE = 8.0, samplingTime = 3.4,  RSEdistFactor = 3.5){
  # Check if filter length is odd for Savitzky-Golay filter.
  if( (SgolayFiltLen %% 2) != 1){
    print("SgolayFiltLen can only be odd number")
    return(NULL)
  }

  # Check if names are consistent between osw and mzML files. Fetch run names.
  filenames <- getRunNames(dataPath, oswMerged, nameCutPattern)
  print("Following runs will be aligned:")
  rownames(filenames) <- paste0("run", 0:(nrow(filenames)-1), "")
  print(filenames[, "runs"])

  # TODO: Make this part in a separate function. Use Environment for pass-by-referene.
  # Get Precursors from the query and respectve chromatogram indices.
  oswFiles <- getMRMoswFiles(filenames, dataPath, peptides = NULL,  query = NULL,
                 oswMerged = oswMerged, maxFdrQuery)

  peptides <- c()
  for(x in oswFiles){
    peptides <- x %>% .$transition_group_id %>% dplyr::union(peptides)
  }

  runs <- filenames$runs
  names(runs) <- rownames(filenames)
  # Collect all the pointers for each mzML file.
  print("Collecting metadata from mzML files.")
  mzPntrs <- list()
  for(mzMLindex in 1:length(runs)){
    run <- names(runs)[mzMLindex]
    filename <- file.path(dataPath, "mzml", paste0(runs[run], ".chrom.mzML"))
    mzPntrs[[mzMLindex]] <- mzR::openMSfile(filename, backend = "pwiz")
  }
  names(mzPntrs) <- names(runs)
  print("Metadata is collected from mzML files.")

  # Initilize output tables.
  rtTbl <- matrix(NA, nrow = length(peptides), ncol = length(runs))
  intesityTbl <- matrix(NA, nrow = length(peptides), ncol = length(runs))
  lwTbl <- matrix(NA, nrow = length(peptides), ncol = length(runs))
  rwTbl <- matrix(NA, nrow = length(peptides), ncol = length(runs))
  rownames(rtTbl) <- peptides; colnames(rtTbl) <- names(runs)
  rownames(intesityTbl) <- peptides; colnames(intesityTbl) <- names(runs)
  rownames(lwTbl) <- peptides; colnames(lwTbl) <- names(runs)
  rownames(rwTbl) <- peptides; colnames(rwTbl) <- names(runs)

  print("Performing reference-based alignment.")
  start_time <- Sys.time()
  pdf("MRM.pdf")
  AlignObjs <- list()
  for(pepIdx in 1:length(peptides)){
    peptide <- peptides[pepIdx]
    # Pick a reference run.
    ref <- "run1"
    exps <- setdiff(names(runs), ref)
    chromIndices <- oswFiles[[ref]] %>%
      dplyr::filter(transition_group_id == peptide) %>% .$chromatogramIndex
    chromIndices <- as.integer(strsplit(chromIndices, split = ",")[[1]])
    productMZ.ref <- oswFiles[[ref]] %>%
      dplyr::filter(transition_group_id == peptide) %>% .$productMZ
    productMZ.ref <- as.integer(strsplit(productMZ.ref, split = ",")[[1]])
    XICs.ref <- extractXIC_group(mzPntrs[[ref]], chromIndices, SgolayFiltOrd, SgolayFiltLen)
    # Align all runs to reference run
    for(eXp in exps){
      # Get XIC_group from experiment run
      chromIndices <- oswFiles[[eXp]] %>%
        dplyr::filter(transition_group_id == peptide) %>% .$chromatogramIndex
      if(length(chromIndices) > 0){
        chromIndices <- as.integer(strsplit(chromIndices, split = ",")[[1]])
        # Make sure that order of fragment-ions is same.
        productMZ.eXp <- oswFiles[[eXp]] %>%
          dplyr::filter(transition_group_id == peptide) %>% .$productMZ
        productMZ.eXp <- as.integer(strsplit(productMZ.eXp, split = ",")[[1]])
        if(!all(productMZ.ref == productMZ.eXp)){
          print(paste("Order of transitions is not same across run for peptide"), peptide)
          chromIndices <- chromIndices[match(productMZ.ref, productMZ.eXp)]
        }
        XICs.eXp <- extractXIC_group(mzPntrs[[eXp]], chromIndices, SgolayFiltOrd, SgolayFiltLen)
        # Get the loess fit for hybrid alignment
        pair <- paste(ref,eXp, sep = "_")
        # Perform dynamic programming for chromatogram alignment
        intensityList.ref <- lapply(XICs.ref, `[[`, 2) # Extracting intensity values
        intensityList.eXp <- lapply(XICs.eXp, `[[`, 2) # Extracting intensity values
        tVec.ref <- XICs.ref[[1]][["time"]] # Extracting time component
        tVec.eXp <- XICs.eXp[[1]][["time"]] # Extracting time component
        Alignobj <- alignChromatogramsCpp(intensityList.ref, intensityList.eXp,
                                          alignType = alignType, tVec.ref, tVec.eXp,
                                          normalization = normalization, simType = simMeasure,
                                          goFactor = goFactor, geFactor = geFactor,
                                          gapQuantile = gapQuantile)
        AlignedIndices <- cbind(Alignobj@indexA_aligned,
                                Alignobj@indexB_aligned,
                                Alignobj@score)
        colnames(AlignedIndices) <- c("indexAligned.ref", "indexAligned.eXp", "score")
        AlignedIndices[, 1:2][AlignedIndices[, 1:2] == 0] <- NA
        tAligned.ref <- mapIdxToTime(tVec.ref, AlignedIndices[,"indexAligned.ref"])
        tAligned.eXp <- mapIdxToTime(tVec.eXp, AlignedIndices[,"indexAligned.eXp"])
        # Map retention time from reference to eXp.
        AlignObjs[[1]] <- list()
        AlignObjs[[1]][[1]] <- Alignobj
        AlignObjs[[1]][[runs[ref]]] <- XICs.ref
        AlignObjs[[1]][[runs[eXp]]] <- XICs.eXp
        AlignObjs[[1]][[4]] <- oswFiles[[ref]] %>%
          dplyr::filter(transition_group_id == peptide) %>%
          dplyr::select(leftWidth, RT, rightWidth) %>%
          as.vector()
      }
      else{
        AlignObjs[[1]] <- NULL
      }
      names(AlignObjs)[1] <- peptide
    }
    plotAlignedPeptides(AlignObjs)
  }
  dev.off()
  end_time <- Sys.time()
  end_time - start_time
  rm(mzPntrs)
  rm(oswFiles)
  rm(rtTbl, intesityTbl, lwTbl, rwTbl)
  print("Data matrix is available in the current directory")
  return(1)
}

#' AlignObj for peptides between a pair for MRM runs.
#'
#' In MRM runs, since XICs have less interference, higher quantile value from similarity
#' matrix is needed to get base-gap-penalty.
#' @return A list of AlignObj. Each AlignObj contains alignment path, similarity matrix and related parameters.
#' @export
getMRMAlignObjs <- function(peptides, runs, dataPath = ".", alignType = "local",
                         query = NULL, oswMerged = FALSE, nameCutPattern = "(.*)(/)(.*)",
                         maxFdrQuery = 0.05, maxFdrLoess = 0.01, spanvalue = 0.1,
                         normalization = "mean", simMeasure = "dotProductMasked",
                         SgolayFiltOrd = 2, SgolayFiltLen = 3,
                         goFactor = 0.125, geFactor = 400,
                         cosAngleThresh = 0.3, OverlapAlignment = TRUE,
                         dotProdThresh = 0.96, gapQuantile = 0.98,
                         hardConstrain = FALSE, samples4gradient = 100,
                         expRSE = 8.0, samplingTime = 3.4,  RSEdistFactor = 3.5){
  if(length(runs)!= 2){
    print("For pairwise alignment, two runs are required.")
    return(NULL)
  }

  if( (SgolayFiltLen %% 2) != 1){
    print("SgolayFiltLen can only be odd number")
    return(NULL)
  }
  # Check if names are consistent between osw and mzML files. Fetch run names.
  filenames <- getRunNames(dataPath, oswMerged, nameCutPattern)
  filenames <- filenames[filenames$runs %in% runs,]
  rownames(filenames) <- paste0("run", 0:(length(runs)-1), "")

  # Get Precursors from the query and respective chromatogram indices.
  oswFiles <- getMRMoswFiles(filenames, dataPath, peptides = NULL,  query = NULL,
                             oswMerged = oswMerged, maxFdrQuery)
  PeptidesFound <- c()
  for(x in oswFiles){
    PeptidesFound <- x %>% .$transition_group_id %>% dplyr::union(PeptidesFound)
  }

  PeptidesFound <- intersect(peptides, PeptidesFound)
  # Report peptides that are not found
  PeptidesNotFound <- setdiff(peptides, PeptidesFound)
  if(length(PeptidesNotFound)>0){
    messsage(paste(PeptidesNotFound, "not found."))
  }

  ####################### Get XICs ##########################################
  runs <- filenames$runs
  names(runs) <- rownames(filenames)
  # Get Chromatogram for each peptide in each run.
  XICs <- list()
  print("Fetching Extracted-ion chromatograms from runs")
  for(i in 1:length(runs)){
    run <- names(runs)[i]
    mzmlName <- file.path(dataPath, "mzml", paste0(runs[run], ".chrom.mzML"))
    mz <- tryCatch(mzR::openMSfile(mzmlName, backend = "pwiz"),
                   error = function(cond) {
                     c$message <- paste0(c$message,
                                         "If error includes invalid cvParam accession 1002746, use FileConverter from OpenMS to decompress chromatograms")
                     stop(cond)})
    XICs_run <- lapply(1:length(PeptidesFound), function(j){
      chromIndices <- oswFiles[[i]] %>%
        dplyr::filter(transition_group_id == PeptidesFound[j]) %>% .$chromatogramIndex
      if(length(chromIndices) > 0){
        chromIndices <- as.integer(strsplit(chromIndices, split = ",")[[1]])
        XIC_group <- extractXIC_group(mz, chromIndices, SgolayFiltOrd, SgolayFiltLen)
      } else{
        XIC_group <- list()
      }
      return(XIC_group)
    })
    names(XICs_run) <- PeptidesFound
    XICs[[i]] <- XICs_run
    rm(mz)
    print(paste("Fetched Extracted-ion chromatograms from run", runs[run]))
  }
  names(XICs) <- runs

  ####################### Perfrom alignment ##########################################
  AlignObjs <- list()
  print("Perfroming alignment")
  for(pepIdx in 1:length(peptides)){
    peptide <- peptides[pepIdx]
    # Select reference run based on m-score
    # Get XIC_group from reference run
    ref <- "run1"
    exps <- setdiff(names(runs), ref)
    XICs.ref <- XICs[[runs[ref]]][[peptide]]

    # Align experiment run to reference run
    for(eXp in exps){
      # Get XIC_group from experiment run
      XICs.eXp <- XICs[[runs[eXp]]][[peptide]]
      if(length(XICs.eXp) > 0){
        # Get the loess fit for hybrid alignment
        pair <- paste(ref, eXp, sep = "_")
        # Set up constraints for penalizing similarity matrix
        tVec.ref <- XICs.ref[[1]][["time"]] # Extracting time component
        tVec.eXp <- XICs.eXp[[1]][["time"]] # Extracting time component
        # Perform dynamic programming for chromatogram alignment
        intensityList.ref <- lapply(XICs.ref, `[[`, 2) # Extracting intensity values
        intensityList.eXp <- lapply(XICs.eXp, `[[`, 2) # Extracting intensity values
        Alignobj <- alignChromatogramsCpp(intensityList.ref, intensityList.eXp,
                                          alignType = alignType, tVec.ref, tVec.eXp,
                                          normalization = normalization, simType = simMeasure,
                                          goFactor = goFactor, geFactor = geFactor,
                                          gapQuantile = gapQuantile)
        AlignObjs[[pepIdx]] <- list()
        AlignObjs[[pepIdx]][[1]] <- Alignobj
        AlignObjs[[pepIdx]][[runs[ref]]] <- XICs.ref
        AlignObjs[[pepIdx]][[runs[eXp]]] <- XICs.eXp
        AlignObjs[[pepIdx]][[4]] <- oswFiles[[ref]] %>%
          dplyr::filter(transition_group_id == peptide) %>%
          dplyr::select(leftWidth, RT, rightWidth) %>%
          as.vector()
      }
      else {AlignObjs[[pepIdx]] <- NULL}
    }
  }
  names(AlignObjs) <- peptides
  print("Alignment done. Returning AlignObjs")
  return(AlignObjs)
}
