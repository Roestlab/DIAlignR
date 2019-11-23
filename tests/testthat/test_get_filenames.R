context("Read osw mzML files.")

test_that("test_filenamesFromOSW", {
  dataPath <- "../../data/testData"
  expOutput <- data.frame("filename" = c("HLA-Ligand-Atlas/BD-ZH12_Lung_Class-1/dia_files/170407_AM_BD-ZH12_Lung_W_10%_DIA_#1_400-650mz_msms47.mzML",
                                       "HLA-Ligand-Atlas/BD-ZH12_Lung_Class-1/dia_files/170407_AM_BD-ZH12_Lung_W_10%_DIA_#1_400-650mz_msms51.mzML",
                                       "HLA-Ligand-Atlas/BD-ZH12_BoneMarrow_Class-1/dia_files/170413_AM_BD-ZH12_BoneMarrow_W_10%_DIA_#1_400-650mz_msms67.mzML"),
                        stringsAsFactors=FALSE)
  expect_identical(filenamesFromOSW(dataPath, "*.osw"), expOutput)
  expect_identical(filenamesFromOSW(dataPath, "*.merged.osw"), expOutput)
  expect_message(filenamesFromOSW(dataPath, "*.mzML"), "Only .osw and .merged.osw files can be read.")
})

test_that("test_filenamesFromMZML", {
  dataPath <- "../../data/testData"
  expOutput <- c("170407_AM_BD-ZH12_Lung_W_10%_DIA_#1_400-650mz_msms47.chrom.mzML" = "170407_AM_BD-ZH12_Lung_W_10%_DIA_#1_400-650mz_msms47",
               "170407_AM_BD-ZH12_Lung_W_10%_DIA_#1_400-650mz_msms51.chrom.mzML" = "170407_AM_BD-ZH12_Lung_W_10%_DIA_#1_400-650mz_msms51",
               "170413_AM_BD-ZH12_BoneMarrow_W_10%_DIA_#1_400-650mz_msms67.chrom.mzML" = "170413_AM_BD-ZH12_BoneMarrow_W_10%_DIA_#1_400-650mz_msms67")
  expect_identical(filenamesFromMZML(dataPath), expOutput)
  expect_message(filenamesFromMZML("."), "0 .chrom.mzML files are found.")
  expect_identical(names(filenamesFromMZML(".")), character(0))
  expect_identical(length(filenamesFromMZML(".")), 0L)
})

test_that("test_getRunNames", {
  dataPath <- "../../data/testData"
  expOutput <- data.frame("filename" = c("HLA-Ligand-Atlas/BD-ZH12_Lung_Class-1/dia_files/170407_AM_BD-ZH12_Lung_W_10%_DIA_#1_400-650mz_msms47.mzML",
                                       "HLA-Ligand-Atlas/BD-ZH12_Lung_Class-1/dia_files/170407_AM_BD-ZH12_Lung_W_10%_DIA_#1_400-650mz_msms51.mzML",
                                       "HLA-Ligand-Atlas/BD-ZH12_BoneMarrow_Class-1/dia_files/170413_AM_BD-ZH12_BoneMarrow_W_10%_DIA_#1_400-650mz_msms67.mzML"),
                        "runs" = c("170407_AM_BD-ZH12_Lung_W_10%_DIA_#1_400-650mz_msms47",
                                   "170407_AM_BD-ZH12_Lung_W_10%_DIA_#1_400-650mz_msms51",
                                   "170413_AM_BD-ZH12_BoneMarrow_W_10%_DIA_#1_400-650mz_msms67"),
                        row.names = c("run0", "run1", "run2"),
                        stringsAsFactors=FALSE)
  expect_identical(getRunNames(dataPath, oswMerged = TRUE, nameCutPattern = "(.*)(/)(.*)"), expOutput)
  expect_identical(getRunNames(dataPath, oswMerged = FALSE, nameCutPattern = "(.*)(/)(.*)"), expOutput)
  expect_error(getRunNames(".", oswMerged = TRUE), "No .merged.osw file is found.")
  expect_error(getRunNames(".", oswMerged = FALSE), "No .osw files are found.")
})
