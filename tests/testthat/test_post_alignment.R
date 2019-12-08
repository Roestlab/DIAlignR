context("Post alignment analysis")

# mzML and osw files are not required.

test_that("test_pickNearestFeature", {
  data(oswFiles_DIAlignR, package="DIAlignR")
  oswFiles <- oswFiles_DIAlignR
  outData <- pickNearestFeature(eXpRT = 5237.8, analyte = "14299_QFNNTDIVLLEDFQK/3", oswFiles,
                                runname = "run2", adaptiveRT = 77.82315,
                                featureFDR = 0.05)
  expData <- list("leftWidth" = c(5217.361),
                           "rightWidth" = c(5275.395),
                           "RT" = c(5240.79),
                           "Intensity" = c(255.496),
                           "peak_group_rank" = c(1L),
                           "m_score" = c(5.692077e-05))
  expect_equal(outData, expData, tolerance = 1e-05)
})

test_that("test_mapIdxToTime", {
  timeVec <- c(1.3,5.6,7.8)
  idx <- c(NA, NA, 1L, 2L, NA, NA, 3L, NA)
  outData <- mapIdxToTime(timeVec, idx)
  expData <- c(1.3, 1.3, 1.3, 5.6, 5.6, 5.6, 7.8, 7.8)
  expect_equal(outData, expData)
})

test_that("test_mappedRTfromAlignObj", {
  AlignObj <- testAlignObj()
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR
  tVec.ref <-  XICs[["run1"]][["14299_QFNNTDIVLLEDFQK/3"]][[1]][, "time"]
  tVec.eXp <-  XICs[["run2"]][["14299_QFNNTDIVLLEDFQK/3"]][[1]][, "time"]
  expect_equal(mappedRTfromAlignObj(refRT= 5238.35, tVec.ref, tVec.eXp, AlignObj), 5237.8)
})
