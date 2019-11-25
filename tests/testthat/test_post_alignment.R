context("Post alignment analysis")

# mzML and osw files are not required.

test_that("test_pickNearestFeature", {
  data(oswFiles_DIAlignR, package="DIAlignR")
  oswFiles <- oswFiles_DIAlignR
  outData <- pickNearestFeature(eXpRT = 4388.1, analyte = "KLYAGAILEV_2", oswFiles,
                                runname = "run0", adaptiveRT = 90.56635,
                                featureFDR = 0.05)
  expData <- list("leftWidth" = c(4369.325),
                           "rightWidth" = c(4413.365),
                           "RT" = c(4390.35),
                           "Intensity" = c(6644520),
                           "peak_group_rank" = c(1L),
                           "m_score" = c(0.002128436))
  expect_equal(outData, expData, tolerance = 1e-05)
})
