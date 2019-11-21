#' Selects only chromatogramId, chromatogramIndex columns and convert them into integer.
#'
#' @return Invisible NULL
chromatogramIdAsInteger <- function(chromHead){
  chromHead <- dplyr::select(chromHead, chromatogramId, chromatogramIndex)
  chromHead <- dplyr::mutate(chromHead, chromatogramId = as.integer(chromatogramId))
  chromHead
}


#' Merge dataframes from OSW and mzML files.
#'
#' @importFrom dplyr %>%
#' @return Update list of data-frame
mergeOswAnalytes_ChromHeader <- function(oswAnalytes, chromHead, analyteFDR, runType = "DIA_proteomics"){
  # TODO: Make sure that transition_id has same order across runs. IMO should be specified in query.
  oswAnalytes <- dplyr::left_join(oswAnalytes, chromHead,
                                   by = c("transition_id" = "chromatogramId"))
  oswAnalytes <- oswAnalytes %>%
    dplyr::group_by(transition_group_id, peak_group_rank) %>%
    dplyr::mutate(transition_ids = paste0(transition_id, collapse = ","),
                  chromatogramIndex = paste0(chromatogramIndex, collapse = ",")) %>%
    dplyr::ungroup() %>% dplyr::select(-transition_id) %>% dplyr::distinct()
  peptides <- dplyr::filter(oswAnalytes, m_score < analyteFDR) %>% .$transition_group_id %>%
    dplyr::union(peptides)
  oswAnalytes
}

#' Get list of peptides and their chromatogram indices.
#'
#' @return A list of data-frames.
#' @export
getOswFiles <- function(dataPath, filenames, maxFdrQuery = 0.05, analyteFDR = 0.01, oswMerged,
                         peptides = NULL, runType = "DIA_proteomics"){
  oswFiles <- list()
  for(i in 1:nrow(filenames)){
    run <- rownames(filenames)[i]
    # Get a query to search against the osw files.
    if(oswMerged == TRUE){
      oswName <- list.files(path = file.path(dataPath, "osw"), pattern="*merged.osw")
      oswName <- file.path(dataPath, "osw", oswName[1])
    } else{
      oswName <- paste0(file.path(dataPath, "osw", filenames$runs[i]), ".osw")
    }
    # Get transition indices for MS2 fragment-ions.
    oswAnalytes <- fetchAnalytesInfo(oswName, maxFdrQuery, oswMerged, peptides = NULL,
                                     filename = filenames$filename[i], runType)

    # Get chromatogram indices from the header file.
    mzmlName <- file.path(dataPath, "mzml", paste0(filenames$runs[i], ".chrom.mzML"))
    chromHead <- readChromatogramHeader(mzmlName)
    chromHead <- chromatogramIdAsInteger(chromHead)
    # Merge chromatogram indices with transition indices and save them.
    # Following function merges analytesInfo dataframe with the chromatogram Header.
    oswAnalytes <- mergeOswAnalytes_ChromHeader(oswAnalytes, chromHead, analyteFDR, runType)
    oswFiles[[i]] <- oswAnalytes
    message("Fetched chromatogram indices from ", filenames$filename[i])
  }
  # Assign rownames to the each element of list
  names(oswFiles) <- rownames(filenames)
  return(oswFiles)
}
