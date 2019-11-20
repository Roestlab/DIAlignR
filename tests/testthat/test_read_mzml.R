context("Read mzml files.")

test_that("test_readChromatogramHeader",{
  mzmlName <- "../../data/testData/mzml/170407_AM_BD-ZH12_Lung_W_10%_DIA_#1_400-650mz_msms47.chrom.mzML"
  expect_error(readChromatogramHeader(mzmlName))
})


test_that("test_readChromatogramHeader",{
  expOutput <- data.frame("chromatogramId" = c("130110", "154511"),
                          "chromatogramIndex" = c(1L,2L),
                          "polarity" = c(-1L,-1L),
                          "precursorIsolationWindowTargetMZ" = c(422.7318, 417.9083),
                          "precursorIsolationWindowLowerOffset" = c(0,0),
                          "precursorIsolationWindowUpperOffset" = c(0,0),
                          "precursorCollisionEnergy" = c(NA_real_, NA_real_),
                          "productIsolationWindowTargetMZ" = c(286.1761, 286.6341),
                          "productIsolationWindowLowerOffset" = c(0,0),
                          "productIsolationWindowUpperOffset" = c(0,0),
                          stringsAsFactors=FALSE)
  mzmlName <- "../../data/testData2/mzml/170407_AM_BD-ZH12_Spleen_W_10%_DIA_#1_400-650mz_msms41.chrom.mzML"
  expect_identical(dim(readChromatogramHeader(mzmlName)), c(175525L, 10L))
  expect_equal(readChromatogramHeader(mzmlName)[1:2,], expOutput, tolerance=1e-6)
})
