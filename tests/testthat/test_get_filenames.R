context("Read osw mzML files.")

test_that("test_filenamesFromOSW", {
  dataPath <- "../../data/example"
  expOutput <- data.frame("filename" = c("data/raw/hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt.mzML.gz",
                                       "data/raw/hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt.mzML.gz",
                                       "data/raw/hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.mzML.gz"),
                        stringsAsFactors=FALSE)
  expect_identical(filenamesFromOSW(dataPath, "*.osw"), expOutput)
  expect_identical(filenamesFromOSW(dataPath, "*merged.osw"), expOutput)
  expect_message(filenamesFromOSW(dataPath, "*.mzML"), "Only .osw and merged.osw files can be read.")
})

test_that("test_filenamesFromMZML", {
  dataPath <- "../../data/example"
  expOutput <- c("hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt.chrom.mzML" = "hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt",
               "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML" = "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
               "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML" = "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
  expect_identical(filenamesFromMZML(dataPath), expOutput)
  expect_message(filenamesFromMZML("."), "0 .chrom.mzML files are found.")
  expect_identical(names(filenamesFromMZML(".")), character(0))
  expect_identical(length(filenamesFromMZML(".")), 0L)
})

test_that("test_getRunNames", {
  dataPath <- "../../data/example"
  expOutput <- data.frame("filename" = c("data/raw/hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt.mzML.gz",
                                         "data/raw/hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt.mzML.gz",
                                         "data/raw/hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.mzML.gz"),
                          "runs" = c("hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt",
                                     "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
                                     "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"),
                          row.names = c("run0", "run1", "run2"),
                          stringsAsFactors=FALSE)
  expect_identical(getRunNames(dataPath, oswMerged = TRUE, nameCutPattern = "(.*)(/)(.*)"), expOutput)
  expect_error(getRunNames(".", oswMerged = TRUE), "No merged.osw file is found.")
  expect_error(getRunNames(".", oswMerged = FALSE), "No .osw files are found.")
})
