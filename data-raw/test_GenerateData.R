context("Generate test data")

# save(ascii = TRUE) : File will be a text file instead of binary.
# saved data is compressed, thats why size is different than pryr::object_size
generateDIAlignRdata <- function(){
  dataPath <- system.file("extdata", package = "DIAlignR")
  oswMerged <- TRUE
  filenames <- getRunNames(dataPath, oswMerged)
  oswFiles_DIAlignR <- getFeatures(filenames, maxFdrQuery = 0.05, runType = "DIA_proteomics")
  save(oswFiles_DIAlignR, file = "oswFiles_DIAlignR.rda", version = 2)

  analytes <- "14299_QFNNTDIVLLEDFQK/3"
  runs <- c("run0" = "hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt",
            "run1" = "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
            "run2" = "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
  outData <- getXICs4AlignObj(dataPath, runs, oswFiles_DIAlignR, analytes, XICfilter = "none")
  XIC_QFNNTDIVLLEDFQK_3_DIAlignR <- outData
  save(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, file = "XIC_QFNNTDIVLLEDFQK_3_DIAlignR.rda", version = 2)

  precursors <- getPrecursors(filenames, oswMerged = TRUE, runType = "DIA_proteomics")
  precursors_DIAlignR <- precursors
  save(precursors_DIAlignR, file = "precursors_DIAlignR.rda", version = 2)

  multipeptide <- getMultipeptide(precursors, oswFiles_DIAlignR)
  multipeptide_DIAlignR <- multipeptide
  save(multipeptide_DIAlignR, file = "multipeptide_DIAlignR.rda", version = 2)
}

generateDIAlignRchrom <- function(){
  transition_group_id <- c("19051_KLIVTSEGC[160]FK/2", "19052_KLIVTSEGC[160]FK/3",
  "1182_GLPIVNLLK/2", "997_SGEISLSSWEN/2", "13597_VVAGGELFKESVVVNDK/3", "4080_VITMPAGVELTNNNNVITVK/3",
  "4081_VITM[147]PAGVELTNNNNVITVK/3", "12300_IHFLSPVRPFTLTPGDEEESFIQLITPVR/3", "14299_QFNNTDIVLLEDFQK/3",
  "4731_GEANVELTPELAFK/2", "7351_ANAMGIPSLTVTNVPGSTLSR/3", "8704_HEQIDLFDSPNEGTR/2")
  runs <- c("hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW", "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW",
            "hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW")
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = paste0(runs[1], ".osw"))
  # Generate a query.
  query <- "SELECT PRECURSOR.GROUP_LABEL AS transition_group_id,
  TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID AS transition_id
  FROM PRECURSOR
  INNER JOIN TRANSITION_PRECURSOR_MAPPING ON TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID = PRECURSOR.ID
  WHERE transition_group_id  = ?"
  IDquery <- DBI::dbSendQuery(con, query)
  DBI::dbBind(IDquery, list(transition_group_id))
  NativeIDS <- DBI::dbFetch(IDquery)
  DBI::dbClearResult(IDquery)
  DBI::dbDisconnect(con)
  write.table(NativeIDS, file = "NativeIDs.csv", sep = ",", row.names= FALSE)
  }

generateAlignObj <- function(){
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  data(oswFiles_DIAlignR, package="DIAlignR")
  XICs.ref <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["run1"]][["14299_QFNNTDIVLLEDFQK/3"]]
  XICs.eXp <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["run2"]][["14299_QFNNTDIVLLEDFQK/3"]]
  globalFit <- getGlobalAlignment(oswFiles_DIAlignR, ref = "run1", eXp = "run2",
                                  maxFdrGlobal = 0.05, spanvalue = 0.1)
  alignObj <- getAlignObj(XICs.ref, XICs.eXp, globalFit, alignType = "hybrid", adaptiveRT = 77.82315,
                          normalization = "mean", simType = "dotProductMasked", goFactor = 0.125,
                          geFactor = 40, cosAngleThresh = 0.3, OverlapAlignment = TRUE, dotProdThresh = 0.96,
                          gapQuantile = 0.5, hardConstrain = FALSE, kerLen = 9, samples4gradient = 100, objType = "heavy")
  alignObj_DIAlignR <- alignObj
  save(alignObj_DIAlignR, file = "alignObj_DIAlignR.rda", version = 2)
}

generateMasterXICs <- function(){
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  data(alignObj_DIAlignR, package="DIAlignR")
  XICs.ref <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["run1"]][["14299_QFNNTDIVLLEDFQK/3"]]
  XICs.eXp <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["run2"]][["14299_QFNNTDIVLLEDFQK/3"]]
  alignedIndices <- cbind(alignObj_DIAlignR@indexA_aligned, alignObj_DIAlignR@indexB_aligned)
  colnames(alignedIndices) <- c("indexAligned.ref", "indexAligned.eXp")
  alignedIndices[, 1:2][alignedIndices[, 1:2] == 0] <- NA_integer_
  newXICs <- childXICs(XICs.ref, XICs.eXp, alignedIndices)

  masterXICs_DIAlignR <- newXICs
  save(masterXICs_DIAlignR, file = "masterXICs_DIAlignR.rda", version = 2)
}

generateMergedOsw <- function(){
  transition_group_id <- c("19051_KLIVTSEGC[160]FK/2", "19052_KLIVTSEGC[160]FK/3",
                           "1182_GLPIVNLLK/2", "997_SGEISLSSWEN/2", "13597_VVAGGELFKESVVVNDK/3", "4080_VITMPAGVELTNNNNVITVK/3",
                           "4081_VITM[147]PAGVELTNNNNVITVK/3", "12300_IHFLSPVRPFTLTPGDEEESFIQLITPVR/3", "14299_QFNNTDIVLLEDFQK/3",
                           "4731_GEANVELTPELAFK/2", "7351_ANAMGIPSLTVTNVPGSTLSR/3", "8704_HEQIDLFDSPNEGTR/2")
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = "merged.osw")
  query <- "SELECT PRECURSOR.GROUP_LABEL AS transition_group_id FROM PRECURSOR WHERE PRECURSOR.DECOY = 0"
  IDquery <- DBI::dbSendQuery(con, query)
  IDs <- DBI::dbFetch(IDquery)
  DBI::dbClearResult(IDquery)
  set.seed(10)
  IDs <- sample(IDs$transition_group_id, 300)
  transition_group_id <- union(IDs, transition_group_id)

  query <- "SELECT PRECURSOR.GROUP_LABEL AS transition_group_id FROM PRECURSOR WHERE PRECURSOR.DECOY = 1"
  IDquery <- DBI::dbSendQuery(con, query)
  IDs <- DBI::dbFetch(IDquery)
  DBI::dbClearResult(IDquery)
  set.seed(10)
  IDs <- sample(IDs$transition_group_id, 10)
  transition_group_id <- union(IDs, transition_group_id)
  DBI::dbDisconnect(con)

  TEMP <- data.frame("ID" = transition_group_id)
  query1 <- "DELETE FROM PRECURSOR WHERE PRECURSOR.GROUP_LABEL NOT IN
      ( SELECT PRECURSOR.GROUP_LABEL
        FROM PRECURSOR
        INNER JOIN TEMP ON PRECURSOR.GROUP_LABEL = TEMP.ID)"

  query2 <- "DELETE FROM PRECURSOR_PEPTIDE_MAPPING
  WHERE PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID NOT IN
    (SELECT PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID
    FROM PRECURSOR
    INNER JOIN TEMP ON PRECURSOR.GROUP_LABEL = TEMP.ID
    INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID = PRECURSOR.ID)"

  query3 <- "DELETE FROM PEPTIDE
  WHERE PEPTIDE.ID NOT IN
    (SELECT PEPTIDE.ID
    FROM PRECURSOR
    INNER JOIN TEMP ON PRECURSOR.GROUP_LABEL = TEMP.ID
    INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID = PRECURSOR.ID
    LEFT JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID)"

  query4 <- "DELETE FROM SCORE_PEPTIDE
  WHERE SCORE_PEPTIDE.PEPTIDE_ID NOT IN
    (SELECT SCORE_PEPTIDE.PEPTIDE_ID
    FROM PRECURSOR
    INNER JOIN TEMP ON PRECURSOR.GROUP_LABEL = TEMP.ID
    INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID = PRECURSOR.ID
    LEFT JOIN PEPTIDE ON PEPTIDE.ID = PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID
    INNER JOIN SCORE_PEPTIDE ON SCORE_PEPTIDE.PEPTIDE_ID = PEPTIDE.ID )"

  query5 <- "DELETE FROM FEATURE WHERE FEATURE.ID NOT IN
    (SELECT FEATURE.ID
    FROM PRECURSOR
    INNER JOIN TEMP ON PRECURSOR.GROUP_LABEL = TEMP.ID
    INNER JOIN FEATURE ON PRECURSOR.ID = FEATURE.PRECURSOR_ID )"

  query6 <- "DELETE FROM FEATURE_MS2 WHERE FEATURE_MS2.FEATURE_ID NOT IN
   (SELECT FEATURE_MS2.FEATURE_ID
    FROM PRECURSOR
    INNER JOIN TEMP ON PRECURSOR.GROUP_LABEL = TEMP.ID
    INNER JOIN FEATURE ON FEATURE.PRECURSOR_ID = PRECURSOR.ID
    LEFT JOIN FEATURE_MS2 ON FEATURE_MS2.FEATURE_ID = FEATURE.ID )"

  query7 <- "DELETE FROM SCORE_MS2 WHERE SCORE_MS2.FEATURE_ID NOT IN
    (SELECT SCORE_MS2.FEATURE_ID
    FROM PRECURSOR
    INNER JOIN TEMP ON PRECURSOR.GROUP_LABEL = TEMP.ID
    INNER JOIN FEATURE ON FEATURE.PRECURSOR_ID = PRECURSOR.ID
    LEFT JOIN FEATURE_MS2 ON FEATURE_MS2.FEATURE_ID = FEATURE.ID
    LEFT JOIN SCORE_MS2 ON SCORE_MS2.FEATURE_ID = FEATURE.ID )"

  query8 <- "DELETE FROM FEATURE_TRANSITION WHERE FEATURE_TRANSITION.FEATURE_ID NOT IN
   (SELECT FEATURE_TRANSITION.FEATURE_ID
    FROM PRECURSOR
    INNER JOIN TEMP ON PRECURSOR.GROUP_LABEL = TEMP.ID
    INNER JOIN FEATURE ON FEATURE.PRECURSOR_ID = PRECURSOR.ID
    LEFT JOIN FEATURE_TRANSITION ON FEATURE_TRANSITION.FEATURE_ID = FEATURE.ID )"

  query9 <- "DELETE FROM TRANSITION_PRECURSOR_MAPPING WHERE TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID NOT IN
    (SELECT TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID
    FROM PRECURSOR
    INNER JOIN TEMP ON PRECURSOR.GROUP_LABEL = TEMP.ID
    LEFT JOIN TRANSITION_PRECURSOR_MAPPING ON TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID = PRECURSOR.ID )"

  query10 <- "DELETE FROM TRANSITION WHERE TRANSITION.ID NOT IN
    (SELECT TRANSITION.ID
    FROM PRECURSOR
    INNER JOIN TEMP ON PRECURSOR.GROUP_LABEL = TEMP.ID
    LEFT JOIN TRANSITION_PRECURSOR_MAPPING ON TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID = PRECURSOR.ID
    LEFT JOIN TRANSITION ON TRANSITION.ID = TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID )"

  query11 <- "DELETE FROM PEPTIDE_PROTEIN_MAPPING
  WHERE PEPTIDE_PROTEIN_MAPPING.PEPTIDE_ID NOT IN
    (SELECT PEPTIDE_PROTEIN_MAPPING.PEPTIDE_ID
    FROM PRECURSOR
    INNER JOIN TEMP ON PRECURSOR.GROUP_LABEL = TEMP.ID
    INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID = PRECURSOR.ID
    LEFT JOIN PEPTIDE ON PEPTIDE.ID = PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID
    LEFT JOIN PEPTIDE_PROTEIN_MAPPING ON PEPTIDE_PROTEIN_MAPPING.PEPTIDE_ID = PEPTIDE.ID )"

  query12 <- "DELETE FROM PROTEIN
  WHERE PROTEIN.ID NOT IN
    (SELECT PROTEIN.ID
    FROM PRECURSOR
    INNER JOIN TEMP ON PRECURSOR.GROUP_LABEL = TEMP.ID
    INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID = PRECURSOR.ID
    INNER JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID
    LEFT JOIN PEPTIDE_PROTEIN_MAPPING ON PEPTIDE_PROTEIN_MAPPING.PEPTIDE_ID = PEPTIDE.ID
    LEFT JOIN PROTEIN ON PROTEIN.ID = PEPTIDE_PROTEIN_MAPPING.PROTEIN_ID )"

  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = "merged.osw")
  DBI::dbBegin(con)
  DBI::dbWriteTable(con, "TEMP", TEMP)
  DBI::dbExecute(con, query1)
  DBI::dbExecute(con, query2)
  DBI::dbExecute(con, query3)
  DBI::dbExecute(con, query4)
  DBI::dbExecute(con, query5)
  DBI::dbExecute(con, query6)
  DBI::dbExecute(con, query7)
  DBI::dbExecute(con, query8)
  DBI::dbExecute(con, query9)
  DBI::dbExecute(con, query10)
  DBI::dbExecute(con, query11)
  DBI::dbExecute(con, query12)
  DBI::dbRemoveTable(con, "TEMP")
  #DBI::dbRollback(con)
  DBI::dbCommit(con)
  DBI::dbExecute(con, "VACUUM")
  DBI::dbDisconnect(con)
}

