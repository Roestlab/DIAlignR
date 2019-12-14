context("Fetching analytes names.")

test_that("test_getAnalytesName",{
  filenames <- data.frame("filename" = c("data/raw/hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt.mzML.gz",
                                         "data/raw/hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt.mzML.gz",
                                         "data/raw/hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.mzML.gz"),
                          "runs" = c("hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt",
                                     "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
                                     "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"),
                          row.names = c("run0", "run1", "run2"),
                          stringsAsFactors=FALSE)
  oswFiles <- getOswFiles(dataPath = "../../inst/extdata", filenames, maxFdrQuery = 0.05, analyteFDR = 0.01,
                                   oswMerged = TRUE, analytes = NULL, runType = "DIA_proteomics")
  outData <- getAnalytesName(oswFiles, analyteFDR = 0.01, commonAnalytes = TRUE)
  expect_identical(length(outData), 137L)
  expect_identical(head(outData, 3), c("AAMIGGADATSNVR_2", "ADAFPGSLSGGQK_2", "ADKAGNVQALIGK_2"))
  outData <- getAnalytesName(oswFiles, analyteFDR = 0.01, commonAnalytes = FALSE)
  expect_identical(length(outData), 174L)
  expect_identical(tail(outData, 3), c("YVGFMNTAE_2", "MSDVQYNPTEQGIVMNKK_3", "WKEGAETLTPSLDLVGK_2"))
  outData <- getAnalytesName(oswFiles, analyteFDR = 1e-05, commonAnalytes = TRUE)
  expect_identical(outData, NULL)
})

test_that("test_getAnalytes",{
  dataPath <- "../../inst/extdata"
  runs <- c("hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
            "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
  outData <- getAnalytes(dataPath, runs, oswMerged = TRUE, maxFdrQuery = 0.01, commonAnalytes = TRUE)
  expect_identical(head(outData,3), c("AAMIGGADATSNVR_2", "ADAFPGSLSGGQK_2", "ADKAGNVQALIGK_2"))
  expect_identical(tail(outData,3), c("YGFDRNEVIMVGDQLMTDIR_3", "YLSEIISAR_2", "YTPSQVAVATFTLAVNR_3"))
  expect_identical(length(outData), 144L)
  outData <- getAnalytes(dataPath, runs, oswMerged = TRUE, maxFdrQuery = 0.01, commonAnalytes = FALSE)
  expect_identical(tail(outData,3), c("VKDIIDEMNISDMTAR_3", "YITGAFAEYDLAK_3", "YVGFMNTAE_2"))
  expect_identical(length(outData), 172L)
})
