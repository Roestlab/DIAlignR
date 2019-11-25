context("Generate test data")

generateDIAlignRdata <- function(){
  dataPath <- "../../data/testData2/"
  oswMerged <- FALSE
  nameCutPattern <- "(.*)(/)(.*)"
  filenames <- getRunNames(dataPath, oswMerged, nameCutPattern)
  oswFiles_DIAlignR <- getOswFiles(dataPath, filenames, maxFdrQuery = 0.05, analyteFDR = 0.01,
                                   oswMerged, peptides = NULL, runType = "DIA_proteomics")
  refAnalytes <- getAnalytesName(oswFiles_DIAlignR, 1.0, commonAnalytes = TRUE)
  set.seed(1)
  refAnalytes <- sample(refAnalytes, 400)
  refAnalytes <- c(refAnalytes,"KLYAGAILEV_2")
  for (i in seq_along(oswFiles_DIAlignR)){
    keep <- which(oswFiles_DIAlignR[[i]]$transition_group_id %in% refAnalytes)
    keep <- c(keep, sample(which(!oswFiles_DIAlignR[[i]]$transition_group_id %in% refAnalytes), 100))
    oswFiles_DIAlignR[[i]] <- oswFiles_DIAlignR[[i]][keep,]
  }
  oswFiles_DIAlignR[[1]]$filename <- "DIAlignR/testData2/dia_files/191124_PM_BD-ZH12_shbh_10%_DIA_#1_400-650mz_msms41.mzML"
  oswFiles_DIAlignR[[2]]$filename <- "DIAlignR/testData2/dia_files/191124_PM_BD-ZH12_shbh_10%_DIA_#1_400-650mz_msms42.mzML"
  oswFiles_DIAlignR[[3]]$filename <- "DIAlignR/testData2/dia_files/191124_PM_BD-ZH12_shbh_10%_DIA_#1_400-650mz_msms35.mzML"
  save(oswFiles_DIAlignR, file = "oswFiles_DIAlignR.RData")

  analytes <- "KLYAGAILEV_2"
  runs <- c("run0" = "191124_PM_BD-ZH12_shbh_10%_DIA_#1_400-650mz_msms41",
            "run1" = "191124_PM_BD-ZH12_shbh_10%_DIA_#1_400-650mz_msms42",
            "run2" = "191124_PM_BD-ZH12_shbh_10%_DIA_#1_400-650mz_msms35",
            stringsAsFactors=FALSE)
  runs <- c("run0" = "170407_AM_BD-ZH12_Spleen_W_10%_DIA_#1_400-650mz_msms41",
            "run1" = "170407_AM_BD-ZH12_Spleen_W_10%_DIA_#2_400-650mz_msms42",
            "run2" = "170413_AM_BD-ZH12_BoneMarrow_W_10%_DIA_#2_400-650mz_msms35")
  outData <- getXICs4AlignObj(dataPath, runs, oswFiles, analytes,
                              SgolayFiltOrd = 4, SgolayFiltLen = 5)
  XIC_KLYAGAILEV_2 <- outData
  save(XIC_KLYAGAILEV_2, file = "XIC_KLYAGAILEV_2_DIAlignR.RData")
}
