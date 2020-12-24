context("Read osw mzML files.")

test_that("test_filenamesFromOSW", {
  dataPath <- system.file("extdata", package = "DIAlignR")
  expOutput <- data.frame("spectraFile" = c("data/raw/hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt.mzML.gz",
                                       "data/raw/hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt.mzML.gz",
                                       "data/raw/hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.mzML.gz"),
                          "spectraFileID" = c("125704171604355508", "6752973645981403097", "2234664662238281994"),
                          "featureFile" = file.path(dataPath, "osw", "merged.osw"),
                        stringsAsFactors=FALSE)
  outData <- filenamesFromOSW(dataPath = dataPath, "*merged.osw$")
  expect_identical(outData, expOutput)
  expect_message(filenamesFromOSW(dataPath = dataPath, "*.mzML"), "Only .osw and merged.osw files can be read.")
})

test_that("test_filenamesFromMZML", {
  dataPath <- system.file("extdata", package = "DIAlignR")
  expOutput <- data.frame("runName" = c("hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt",
                                 "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
                               "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"),
                 "chromatogramFile" = c(file.path(dataPath, "mzml", "hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt.chrom.mzML"),
                                        file.path(dataPath, "mzml", "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML"),
                                        file.path(dataPath, "mzml", "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML")),
                 stringsAsFactors=FALSE)
  expect_identical(filenamesFromMZML(dataPath = dataPath), expOutput)
  expect_message(filenamesFromMZML(dataPath = "."), "0 .chrom.mzML files are found.")
  expect_identical(names(filenamesFromMZML(dataPath = ".")), c("runName", "chromatogramFile"))
  expect_identical(nrow(filenamesFromMZML(dataPath = ".")), 0L)
})

test_that("test_getRunNames", {
  dataPath <- system.file("extdata", package = "DIAlignR")
  expOutput <- data.frame("runName" = c("hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt",
                                        "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
                                        "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"),
                          "spectraFile" = c("data/raw/hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt.mzML.gz",
                                         "data/raw/hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt.mzML.gz",
                                         "data/raw/hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.mzML.gz"),
                          "spectraFileID" = c("125704171604355508", "6752973645981403097", "2234664662238281994"),
                          "featureFile" = file.path(dataPath, "osw", "merged.osw"),
                          "chromatogramFile" = c(file.path(dataPath, "mzml", "hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt.chrom.mzML"),
                                                 file.path(dataPath, "mzml", "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML"),
                                                 file.path(dataPath, "mzml", "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML")),
                          row.names = c("run0", "run1", "run2"),
                          stringsAsFactors=FALSE)
  expect_identical(getRunNames(dataPath = dataPath, oswMerged = TRUE), expOutput)
  expect_error(getRunNames(dataPath = ".", oswMerged = TRUE), "No merged.osw file is found.")
  expect_error(getRunNames(dataPath = ".", oswMerged = FALSE), "No .osw files are found.")
})

test_that("test_updateFileInfo", {

})


test_that("test_addMasterToOSW", {

})
