context("Read mzml files.")

test_that("test_readChromatogramHeader",{
  mzmlName <- file.path(system.file("extdata", package = "DIAlignR"), "mzml", "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML")
  expOutput <- data.frame("chromatogramId" = c("103114", "103115", "103116"),
                          "chromatogramIndex" = c(1L,2L,3L),
                          "polarity" = c(-1L,-1L,-1L),
                          "precursorIsolationWindowTargetMZ" = c(696.037, 696.037, 696.037),
                          "precursorIsolationWindowLowerOffset" = c(0,0,0),
                          "precursorIsolationWindowUpperOffset" = c(0,0,0),
                          "precursorCollisionEnergy" = c(NA_real_, NA_real_, NA_real_),
                          "productIsolationWindowTargetMZ" = c(717.399, 558.277, 1031.560),
                          "productIsolationWindowLowerOffset" = c(0,0,0),
                          "productIsolationWindowUpperOffset" = c(0,0,0),
                          stringsAsFactors=FALSE)
  expect_identical(dim(readChromatogramHeader(mzmlName)), c(72L, 10L))
  expect_equal(readChromatogramHeader(mzmlName)[1:3,1:8], expOutput[,1:8], tolerance=1e-6)
  mzmlName <- "getError"
  expect_error(readChromatogramHeader(mzmlName))
})

test_that("test_getMZMLpointers",{
  dataPath <- system.file("extdata", package = "DIAlignR")
  runs <- c("run0" = "hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt",
           "run1" = "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
           "run2" = "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
  outData <- getMZMLpointers(dataPath = dataPath, runs = runs)
  expect_is(outData[["run0"]], "mzRpwiz")
  expect_is(outData[["run1"]], "mzRpwiz")
  expect_is(outData[["run2"]], "mzRpwiz")
})
