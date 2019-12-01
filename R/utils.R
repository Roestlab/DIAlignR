#' Provides the index of reference run.
#'
#' @importFrom dplyr %>%
#' @return An integer
getRefRun <- function(oswFiles, analyte){
  # Select reference run based on m-score
  minMscore <- 1
  refRunIdx <- NULL
  for(runIdx in 1:length(oswFiles)){
    m_score <- oswFiles[[runIdx]] %>%
      dplyr::filter(transition_group_id == analyte & peak_group_rank == 1) %>% .$m_score
    # Check for numeric(0) condition and proceed.
    if(length(m_score) != 0){
      if(m_score < minMscore){
        minMscore <- m_score
        refRunIdx <- runIdx
      }
    }
  }
  refRunIdx
}

#' Provides the index of reference run.
#'
#' @importFrom dplyr %>%
#' @return A vector of Integers
selectChromIndices <- function(oswFiles, runname, analyte){
  # Pick chromatrogram indices from osw table.
  chromIndices <- oswFiles[[runname]] %>%
    dplyr::filter(transition_group_id == analyte) %>% .$chromatogramIndex
  # Check for character(0) condition and proceed.
  if(length(chromIndices) != 0){
    chromIndices <- as.integer(strsplit(chromIndices, split = ",")[[1]])
  } else {
    return(NULL)
  }
  # Select the first row if there are many peak-groups.
  chromIndices
}

