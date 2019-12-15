context("Pairwise analyte peak correspondence.")

test_that("test_getAlignObj", {
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR
  data(oswFiles_DIAlignR, package="DIAlignR")
  oswFiles <- oswFiles_DIAlignR
  XICs.ref <- XICs[["run1"]][["14299_QFNNTDIVLLEDFQK/3"]]
  XICs.eXp <- XICs[["run2"]][["14299_QFNNTDIVLLEDFQK/3"]]
  Loess.fit <- getLOESSfit(oswFiles, ref = "run1", eXp = "run2", maxFdrGlobal = 0.05, spanvalue = 0.1)
  outData <- getAlignObj(XICs.ref, XICs.eXp, Loess.fit, adaptiveRT = 77.82315, samplingTime = 3.414,
                         normalization = "mean", simType = "dotProductMasked", goFactor = 0.125, geFactor = 40,
                         cosAngleThresh = 0.3, OverlapAlignment = TRUE,
                         dotProdThresh = 0.96, gapQuantile = 0.5, hardConstrain = FALSE,
                         samples4gradient = 100, objType = "light")
  expData <- testAlignObj()
  expect_equal(outData, expData, tolerance = 1e-03)
})

test_that("test_getMappedRT", {
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR
  data(oswFiles_DIAlignR, package="DIAlignR")
  oswFiles <- oswFiles_DIAlignR
  XICs.ref <- XICs[["run1"]][["14299_QFNNTDIVLLEDFQK/3"]]
  XICs.eXp <- XICs[["run2"]][["14299_QFNNTDIVLLEDFQK/3"]]
  adaptiveRT <- 77.82315 #3.5*Loess.fit$s
  Loess.fit <- getLOESSfit(oswFiles, ref = "run2", eXp = "run0", maxFdrGlobal = 0.05, spanvalue = 0.1)
  outData <- getMappedRT(refRT = 5238.35, XICs.ref, XICs.eXp, Loess.fit, alignType = "hybrid",
                         adaptiveRT = adaptiveRT, samplingTime = 3.414,
                         normalization = "mean", simMeasure = "dotProductMasked",
                         goFactor = 0.125, geFactor = 40, cosAngleThresh = 0.3,
                         OverlapAlignment = TRUE,
                         dotProdThresh = 0.96, gapQuantile = 0.5, hardConstrain = FALSE,
                         samples4gradient = 100)
  expect_equal(outData, 5237.8, tolerance = 1e-03)
})
