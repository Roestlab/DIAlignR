#' Selects only chromatogramId, chromatogramIndex columns and convert them into integer.
#'
#' @return Invisible NULL
chromatogramIdAsInteger <- function(chromatogramHeader){
  assign("chromHead", dplyr::mutate(dplyr::select(chromatogramHeader, chromatogramId, chromatogramIndex),
                                    chromatogramId = as.integer(chromatogramId)),
         envir = parent.frame(n = 1))
  invisible(NULL)
}


#' Merge dataframes from OSW and mzML files.
#'
#' @importFrom dplyr %>%
#' @return Update list of data-frame
mergeOswAnalytes_ChromHeader <- function(oswAnalytes, chromHead, analyteFDR =  1.00, runType = "DIA_proteomics"){
  # TODO: Make sure that transition_id has same order across runs. IMO should be specified in query.
  assign("oswAnalytes", dplyr::left_join(oswAnalytes, chromHead,
                                  by = c("transition_id" = "chromatogramId")) %>%
    dplyr::group_by(transition_group_id, peak_group_rank) %>%
    dplyr::mutate(transition_ids = paste0(transition_id, collapse = ","),
                  chromatogramIndex = paste0(chromatogramIndex, collapse = ",")) %>%
    dplyr::ungroup() %>% dplyr::select(-transition_id) %>% dplyr::distinct(),
    envir = parent.frame(n = 1))
  invisible(NULL)
}

#' Get list of peptides and their chromatogram indices.
#'
#' @return A list of data-frames.
#' @export
getOswFiles <- function(dataPath, filenames, maxFdrQuery = 0.05, analyteFDR = 0.01, oswMerged = TRUE,
                         analytes = NULL, runType = "DIA_proteomics", analyteInGroupLabel = FALSE, masterMzExperiment = NULL){
  oswFiles <- list()
  if ( is.null(masterMzExperiment) ){
    message( "A masterMzExperiment object WAS NOT supplied...")
    masterMzExperiment <- list()
  } else {
    message( "A masterMzExperiment object WAS supplied...")
    print( masterMzExperiment )
  }
  for(i in 1:nrow(filenames)){
    run <- rownames(filenames)[i]
    ## Save current run name in master list
    masterMzExperiment[[run]]$run <- run
    # Get a query to search against the osw files.
    if(oswMerged == TRUE){
      oswName <- list.files(path = file.path(dataPath, "osw"), pattern="*merged.osw")
      oswName <- file.path(dataPath, "osw", oswName[1])
    } else{
      oswName <- paste0(file.path(dataPath, "osw", filenames$runs[i]), ".osw")
    }
    # Get transition indices for MS2 fragment-ions.
    tictoc::tic("getOswFiles -> Fetch Analytes Info, MS2 fragment ion indies... Elapsed time")
    oswAnalytes <- fetchAnalytesInfo(oswName, maxFdrQuery, oswMerged, analytes = analytes,
                                     filename = filenames$filename[i], runType, analyteInGroupLabel)
    tictoc::toc()
    # Get chromatogram indices from the header file.
    mzmlName <- file.path(dataPath, "mzml", paste0(filenames$runs[i], ".chrom.mzML"))
    tictoc::tic("getOswFiles -> readChromatogramHeader: get chrom indices... Elapsed time")
    if ( !('mz' %in% names(masterMzExperiment[[run]]$mzExperiment)) ){
      message('An mz object was not found in the masterMZExperiment object supplied...')
      print( (masterMzExperiment[[run]]$mzExperiment) )
      mzExperiment <- readChromatogramHeader(mzmlName)
    } else {
      if ( !('chromHead' %in% names(masterMzExperiment[[run]]$mzExperiment)) ){
        message('A chromHead data.frame was not found in the masterMZExperiment object supplied...')
        mzExperiment <- readChromatogramHeader(mzmlName, masterMzExperiment[[run]]$mz)
      } else {
        message('An mz object and chromHead was found in the masterMZExperiment object supplied, skipping readChromatogramHeader...')
        mzExperiment <- masterMzExperiment[[run]]$mzExperiment
      }
    }
    chromHead <- mzExperiment$chromHead
    chromatogramIdAsInteger(chromHead)
    ## Save mzExperiment list in master list
    masterMzExperiment[[run]]$mzExperiment <- mzExperiment
    tictoc::toc()
    # Merge chromatogram indices with transition indices and save them.
    # Following function merges analytesInfo dataframe with the chromatogram Header.
    mergeOswAnalytes_ChromHeader(oswAnalytes, chromHead, analyteFDR, runType)
    oswFiles[[i]] <- oswAnalytes
    message("Fetched chromatogram indices from ", filenames$filename[i])
  }
  # Assign rownames to the each element of list
  names(oswFiles) <- rownames(filenames)

  masterMzExperiment$oswFiles <- oswFiles
  return(masterMzExperiment)
}
