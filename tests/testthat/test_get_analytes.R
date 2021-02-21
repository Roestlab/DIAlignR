context("Fetching analytes names.")

test_that("test_analytesFromFeatures",{
  filenames <- data.frame("filename" = c("data/raw/hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt.mzML.gz",
                                         "data/raw/hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt.mzML.gz",
                                         "data/raw/hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.mzML.gz"),
                          "runs" = c("hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt",
                                     "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
                                     "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"),
                          row.names = c("run0", "run1", "run2"),
                          stringsAsFactors=FALSE)
  dataPath <- system.file("extdata", package = "DIAlignR")
  fileInfo <- getRunNames(dataPath = dataPath)
  oswFiles <- getFeatures(fileInfo, maxFdrQuery = 0.05, runType = "DIA_Proteomics")

  outData <- analytesFromFeatures(oswFiles, analyteFDR = 0.01, commonAnalytes = TRUE)
  expect_identical(length(outData), 138L)
  expect_identical(head(outData, 3), c(470L, 898L, 997L))
  outData <- analytesFromFeatures(oswFiles, analyteFDR = 0.01, commonAnalytes = FALSE)
  expect_identical(length(outData), 175L)
  expect_identical(tail(outData, 3), c(18807L, 945L, 19762L))
  outData <- analytesFromFeatures(oswFiles, analyteFDR = 1e-05, commonAnalytes = TRUE)
  expect_identical(outData, NULL)
})
