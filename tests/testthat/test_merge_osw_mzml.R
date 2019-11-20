context("Merging OSW and mzML files.")


test_that("test_mergeOsw_ChromHeader", {
  filenames <- data.frame("filename" = c("HLA-Ligand-Atlas/BD-ZH12_Spleen_Class-1/dia_files/170407_AM_BD-ZH12_Spleen_W_10%_DIA_#1_400-650mz_msms41.mzML",
                                         "HLA-Ligand-Atlas/BD-ZH12_Spleen_Class-1/dia_files/170407_AM_BD-ZH12_Spleen_W_10%_DIA_#2_400-650mz_msms42.mzML",
                                         "HLA-Ligand-Atlas/BD-ZH12_BoneMarrow_Class-1/dia_files/170413_AM_BD-ZH12_BoneMarrow_W_10%_DIA_#2_400-650mz_msms35.mzML"),
                          "runs" = c("170407_AM_BD-ZH12_Spleen_W_10%_DIA_#1_400-650mz_msms41",
                                     "170407_AM_BD-ZH12_Spleen_W_10%_DIA_#2_400-650mz_msms42",
                                     "170413_AM_BD-ZH12_BoneMarrow_W_10%_DIA_#2_400-650mz_msms35"),
                          row.names = c("run0", "run1", "run2"),
                          stringsAsFactors=FALSE)
  oswName <- "../../data/testData2/osw/170407_AM_BD-ZH12_Spleen_W_10%_DIA_#1_400-650mz_msms41.osw"
  mzmlName <- "../../data/testData2/mzml/170407_AM_BD-ZH12_Spleen_W_10%_DIA_#1_400-650mz_msms41.chrom.mzML"
  runType <- "DIA_proteomics"
  analytesInfo <- fetchAnalytesInfo(oswName, maxFdrQuery = 0.05, oswMerged = FALSE,
                                    peptides = NULL,
                                    filename = "HLA-Ligand-Atlas/BD-ZH12_Spleen_Class-1/dia_files/170407_AM_BD-ZH12_Spleen_W_10%_DIA_#1_400-650mz_msms41.mzML",
                                    runType)
  chromHead <- readChromatogramHeader(mzmlName)
  mergeOsw_ChromHeader(chromHead, analytesInfo, peptides, runType)

})
