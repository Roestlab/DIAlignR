#' Generate SQL query to fetch information from osw files.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-14
#' @param maxFdrQuery (numeric) value between 0 and 1. It is used to filter features from osw file which have SCORE_MS2.QVALUE less than itself.
#' @param oswMerged (logical) TRUE for experiment-wide FDR and FALSE for run-specific FDR by pyprophet.
#' @param analytes (vector of strings) transition_group_ids for which features are to be extracted. analyteInGroupLabel must be set according the pattern used here.
#' @param filename (string) as mentioned in RUN table of osw files.
#' @param runType (char) This must be one of the strings "DIA_proteomics", "DIA_Metabolomics".
#' @param analyteInGroupLabel (logical) TRUE for getting analytes as PRECURSOR.GROUP_LABEL from osw file.
#' @return SQL query to be searched.
#' @keywords internal
getQuery <- function(maxFdrQuery, oswMerged = TRUE, analytes = NULL,
                     filename = NULL, runType = "DIA_Proteomics", analyteInGroupLabel = FALSE){
  if(is.null(analytes)){
    selectAnalytes <- ""
  } else{
    selectAnalytes <- paste0(" AND transition_group_id IN ('", paste(analytes, collapse="','"),"')")
  }

  if(oswMerged){
    matchFilename <- paste0(" AND RUN.FILENAME ='", filename,"'")
  } else{
    matchFilename <- ""
  }

  if(analyteInGroupLabel == TRUE){
    transition_group_id <- " PRECURSOR.GROUP_LABEL AS transition_group_id"
  } else {
    transition_group_id <- " PEPTIDE.MODIFIED_SEQUENCE || '_' || PRECURSOR.CHARGE AS transition_group_id"
  }

  if(runType == "DIA_Metabolomics"){
    query <- paste0("SELECT RUN.ID AS id_run,
    COMPOUND.ID AS id_compound,
    COMPOUND.COMPOUND_NAME || '_' || COMPOUND.ADDUCTS AS transition_group_id,
    TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID AS transition_id,
    RUN.ID AS run_id,
    RUN.FILENAME AS filename,
    FEATURE.EXP_RT AS RT,
    FEATURE.DELTA_RT AS delta_rt,
    PRECURSOR.LIBRARY_RT AS assay_RT,
    FEATURE.ID AS id,
    COMPOUND.SUM_FORMULA AS sum_formula,
    COMPOUND.COMPOUND_NAME AS compound_name,
    COMPOUND.ADDUCTS AS Adducts,
    PRECURSOR.CHARGE AS Charge,
    PRECURSOR.PRECURSOR_MZ AS mz,
    FEATURE_MS2.AREA_INTENSITY AS Intensity,
    FEATURE.LEFT_WIDTH AS leftWidth,
    FEATURE.RIGHT_WIDTH AS rightWidth,
    SCORE_MS2.RANK AS peak_group_rank,
    SCORE_MS2.QVALUE AS m_score
    FROM PRECURSOR
    INNER JOIN PRECURSOR_COMPOUND_MAPPING ON PRECURSOR.ID = PRECURSOR_COMPOUND_MAPPING.PRECURSOR_ID
    INNER JOIN COMPOUND ON PRECURSOR_COMPOUND_MAPPING.COMPOUND_ID = COMPOUND.ID
    INNER JOIN FEATURE ON FEATURE.PRECURSOR_ID = PRECURSOR.ID
    INNER JOIN RUN ON RUN.ID = FEATURE.RUN_ID
    LEFT JOIN TRANSITION_PRECURSOR_MAPPING ON TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID = PRECURSOR.ID
    LEFT JOIN FEATURE_MS1 ON FEATURE_MS1.FEATURE_ID = FEATURE.ID
    LEFT JOIN FEATURE_MS2 ON FEATURE_MS2.FEATURE_ID = FEATURE.ID
    LEFT JOIN SCORE_MS2 ON SCORE_MS2.FEATURE_ID = FEATURE.ID
    WHERE COMPOUND.DECOY = 0 AND SCORE_MS2.QVALUE <  ", maxFdrQuery, selectAnalytes, matchFilename, "
    ORDER BY transition_group_id,
    peak_group_rank;")
  } else if (runType == "MRM_Proteomics"){
    query <- paste0("SELECT PEPTIDE.MODIFIED_SEQUENCE || '_' || PRECURSOR.CHARGE AS transition_group_id,
  RUN.FILENAME AS filename,
  FEATURE.EXP_RT AS RT,
  FEATURE.DELTA_RT AS delta_rt,
  PRECURSOR.LIBRARY_RT AS assay_RT,
  FEATURE_MS2.AREA_INTENSITY AS Intensity,
  FEATURE.LEFT_WIDTH AS leftWidth,
  FEATURE.RIGHT_WIDTH AS rightWidth,
  TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID AS transition_id
  FROM PRECURSOR
  INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR.ID = PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID AND PRECURSOR.DECOY=0
  INNER JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID
  INNER JOIN FEATURE ON FEATURE.PRECURSOR_ID = PRECURSOR.ID
  INNER JOIN RUN ON RUN.ID = FEATURE.RUN_ID
  INNER JOIN TRANSITION_PRECURSOR_MAPPING ON TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID = PRECURSOR.ID
  LEFT JOIN FEATURE_MS2 ON FEATURE_MS2.FEATURE_ID = FEATURE.ID
  ORDER BY transition_group_id;")
  } else{
    query <- paste0("SELECT", transition_group_id,",
  RUN.FILENAME AS filename,
  FEATURE.EXP_RT AS RT,
  FEATURE.DELTA_RT AS delta_rt,
  PRECURSOR.LIBRARY_RT AS assay_RT,
  FEATURE_MS2.AREA_INTENSITY AS Intensity,
  FEATURE.LEFT_WIDTH AS leftWidth,
  FEATURE.RIGHT_WIDTH AS rightWidth,
  SCORE_MS2.RANK AS peak_group_rank,
  SCORE_MS2.QVALUE AS m_score,
  TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID AS transition_id
  FROM PRECURSOR
  INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR.ID = PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID AND PRECURSOR.DECOY=0
  INNER JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID
  INNER JOIN FEATURE ON FEATURE.PRECURSOR_ID = PRECURSOR.ID
  INNER JOIN RUN ON RUN.ID = FEATURE.RUN_ID
  INNER JOIN TRANSITION_PRECURSOR_MAPPING ON TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID = PRECURSOR.ID
  LEFT JOIN FEATURE_MS2 ON FEATURE_MS2.FEATURE_ID = FEATURE.ID
  LEFT JOIN SCORE_MS2 ON SCORE_MS2.FEATURE_ID = FEATURE.ID
  WHERE SCORE_MS2.QVALUE < ", maxFdrQuery, selectAnalytes, matchFilename, "
  ORDER BY transition_group_id,
  peak_group_rank;")
  }
  return(query)
}


#' Generate SQL query to fetch limited information from osw files.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-14
#' @param maxFdrQuery (numeric) value between 0 and 1. It is used to filter features from osw file which have SCORE_MS2.QVALUE less than itself.
#' @param oswMerged (logical) TRUE for experiment-wide FDR and FALSE for run-specific FDR by pyprophet.
#' @param filename (string) as mentioned in RUN table of osw files..
#' @param runType (char) This must be one of the strings "DIA_proteomics", "DIA_Metabolomics".
#' @param analyteInGroupLabel (logical) TRUE for getting analytes as PRECURSOR.GROUP_LABEL from osw file.
#' @return SQL query to be searched.
#' @seealso \code{\link{getOswAnalytes}}
#' @keywords internal
getAnalytesQuery <- function(maxFdrQuery, oswMerged = TRUE, filename = NULL,
                             runType = "DIA_Proteomics", analyteInGroupLabel = FALSE){
  if(oswMerged){
    matchFilename <- paste0(" AND RUN.FILENAME ='", filename,"'")
  } else{
    matchFilename <- ""
  }

  if(analyteInGroupLabel == TRUE){
    transition_group_id <- " PRECURSOR.GROUP_LABEL AS transition_group_id"
  } else {
    transition_group_id <- " PEPTIDE.MODIFIED_SEQUENCE || '_' || PRECURSOR.CHARGE AS transition_group_id"
  }

  if(runType == "DIA_Metabolomics"){
    query <- paste0("SELECT COMPOUND.ID AS compound_id,
    COMPOUND.COMPOUND_NAME || '_' || COMPOUND.ADDUCTS AS transition_group_id,
    RUN.FILENAME AS filename,
    SCORE_MS2.RANK AS peak_group_rank,
    SCORE_MS2.QVALUE AS m_score
    FROM PRECURSOR
    INNER JOIN PRECURSOR_COMPOUND_MAPPING ON PRECURSOR.ID = PRECURSOR_COMPOUND_MAPPING.PRECURSOR_ID
    INNER JOIN COMPOUND ON PRECURSOR_COMPOUND_MAPPING.COMPOUND_ID = COMPOUND.ID
    INNER JOIN FEATURE ON FEATURE.PRECURSOR_ID = PRECURSOR.ID
    INNER JOIN RUN ON RUN.ID = FEATURE.RUN_ID
    LEFT JOIN SCORE_MS2 ON SCORE_MS2.FEATURE_ID = FEATURE.ID
    WHERE COMPOUND.DECOY = 0 AND SCORE_MS2.QVALUE <  ", maxFdrQuery, matchFilename, "
    ORDER BY transition_group_id,
    peak_group_rank;")
  } else if (runType == "MRM_Proteomics"){
    query <- paste0("SELECT PEPTIDE.MODIFIED_SEQUENCE || '_' || PRECURSOR.CHARGE AS transition_group_id,
  RUN.FILENAME AS filename,
  FEATURE.EXP_RT AS RT,
  FROM PRECURSOR
  INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR.ID = PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID AND PRECURSOR.DECOY=0
  INNER JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID
  INNER JOIN FEATURE ON FEATURE.PRECURSOR_ID = PRECURSOR.ID
  INNER JOIN RUN ON RUN.ID = FEATURE.RUN_ID
  ORDER BY transition_group_id;")
  } else{
    query <- paste0("SELECT", transition_group_id,",
  RUN.FILENAME AS filename,
  SCORE_MS2.RANK AS peak_group_rank,
  SCORE_MS2.QVALUE AS m_score,
  TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID AS transition_id
  FROM PRECURSOR
  INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR.ID = PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID AND PRECURSOR.DECOY=0
  INNER JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID
  INNER JOIN FEATURE ON FEATURE.PRECURSOR_ID = PRECURSOR.ID
  INNER JOIN RUN ON RUN.ID = FEATURE.RUN_ID
  INNER JOIN TRANSITION_PRECURSOR_MAPPING ON TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID = PRECURSOR.ID
  LEFT JOIN SCORE_MS2 ON SCORE_MS2.FEATURE_ID = FEATURE.ID
  WHERE SCORE_MS2.QVALUE < ", maxFdrQuery, matchFilename, "
  ORDER BY transition_group_id,
  peak_group_rank;")
  }
  return(query)
}

#  https://stackoverflow.com/questions/10622260/how-do-you-query-an-int-column-for-any-value
#' Get precursor Info
#'
#' For each precursor in the table respective transition ids are fetched.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2020-04-04
#' @param runType (char) This must be one of the strings "DIA_proteomics", "DIA_Metabolomics".
#' @return SQL query to be searched.
#' @seealso \code{\link{fetchPrecursorsInfo}}
#' @keywords internal
getPrecursorsQuery <- function(runType = "DIA_Proteomics"){
  query <- "SELECT PRECURSOR.ID AS transition_group_id,
      TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID AS transition_id,
      PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID AS peptide_id,
      PEPTIDE.MODIFIED_SEQUENCE AS sequence,
      PRECURSOR.CHARGE AS charge,
      PRECURSOR.GROUP_LABEL AS group_label
      FROM PRECURSOR
      INNER JOIN TRANSITION_PRECURSOR_MAPPING ON TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID = PRECURSOR.ID
      INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID = PRECURSOR.ID
      INNER JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID
      ORDER BY transition_group_id, transition_id;"
  query
}

#' Get features from a SQLite file
#'
#' Query is generated to identify features below a FDR cut-off from a run.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2020-04-07
#' @param runType (char) This must be one of the strings "DIA_proteomics", "DIA_Metabolomics".
#' @return SQL query to be searched.
#' @seealso \code{\link{fetchFeaturesFromRun}}
#' @keywords internal
getFeaturesQuery <- function(runType = "DIA_Proteomics"){
  query <- "SELECT PRECURSOR.ID AS transition_group_id,
  FEATURE.EXP_RT AS RT,
  FEATURE_MS2.AREA_INTENSITY AS intensity,
  FEATURE.LEFT_WIDTH AS leftWidth,
  FEATURE.RIGHT_WIDTH AS rightWidth,
  SCORE_MS2.RANK AS peak_group_rank,
  SCORE_MS2.QVALUE AS m_score
  FROM PRECURSOR
  INNER JOIN FEATURE ON FEATURE.PRECURSOR_ID = PRECURSOR.ID
  INNER JOIN RUN ON RUN.ID = FEATURE.RUN_ID
  LEFT JOIN FEATURE_MS2 ON FEATURE_MS2.FEATURE_ID = FEATURE.ID
  LEFT JOIN SCORE_MS2 ON SCORE_MS2.FEATURE_ID = FEATURE.ID
  WHERE RUN.ID = $runID AND SCORE_MS2.QVALUE < $FDR
  ORDER BY transition_group_id, peak_group_rank;"
  query
}


#' Get precursor Info
#'
#' For each precursor in the table respective transition ids are fetched.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-04-04
#' @param analytes (integer) A vector of integer that is searched in PRECURSOR.ID.
#' @param runType (char) This must be one of the strings "DIA_proteomics", "DIA_Metabolomics".
#' @return SQL query to be searched.
#' @seealso \code{\link{fetchPrecursorsInfo}}
#' @keywords internal
getPrecursorsQueryID <- function(analytes, runType = "DIA_Proteomics"){
  selectAnalytes <- paste0(" transition_group_id IN ('", paste(analytes, collapse="','"),"')")

  query <- paste0("SELECT PRECURSOR.ID AS transition_group_id,
      TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID AS transition_id,
      PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID AS peptide_id,
      PEPTIDE.MODIFIED_SEQUENCE AS sequence,
      PRECURSOR.CHARGE AS charge,
      PRECURSOR.GROUP_LABEL AS group_label
      FROM PRECURSOR
      INNER JOIN TRANSITION_PRECURSOR_MAPPING ON TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID = PRECURSOR.ID
      INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID = PRECURSOR.ID
      INNER JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID
      WHERE ", selectAnalytes, "
      ORDER BY transition_group_id, transition_id;")
  query
}
