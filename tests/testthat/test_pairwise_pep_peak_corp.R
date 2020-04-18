context("Pairwise analyte peak correspondence.")

test_that("test_getAlignObj", {
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR
  data(oswFiles_DIAlignR, package="DIAlignR")
  oswFiles <- oswFiles_DIAlignR
  XICs.ref <- XICs[["run1"]][["14299_QFNNTDIVLLEDFQK/3"]]
  XICs.ref <- smoothXICs(XICs.ref, type = "sgolay", kernelLen = 13, polyOrd = 4)
  XICs.eXp <- XICs[["run2"]][["14299_QFNNTDIVLLEDFQK/3"]]
  XICs.eXp <- smoothXICs(XICs.eXp, type = "sgolay", kernelLen = 13, polyOrd = 4)
  globalFit <- getLOESSfit(oswFiles, ref = "run1", eXp = "run2", maxFdrGlobal = 0.05, spanvalue = 0.1)
  outData <- getAlignObj(XICs.ref, XICs.eXp, globalFit, alignType = "hybrid", adaptiveRT = 38.6594179136227,
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
                         adaptiveRT = adaptiveRT, normalization = "mean", simMeasure = "dotProductMasked",
                         goFactor = 0.125, geFactor = 40, cosAngleThresh = 0.3,
                         OverlapAlignment = TRUE,
                         dotProdThresh = 0.96, gapQuantile = 0.5, hardConstrain = FALSE,
                         samples4gradient = 100)
  expect_equal(outData, 5237.8, tolerance = 1e-03)
})

test_that("test_getAlignedIndices", {
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR
  data(oswFiles_DIAlignR, package="DIAlignR")
  oswFiles <- oswFiles_DIAlignR
  XICs.ref <- XICs[["run1"]][["14299_QFNNTDIVLLEDFQK/3"]]
  XICs.eXp <- XICs[["run2"]][["14299_QFNNTDIVLLEDFQK/3"]]
  adaptiveRT <- 77.82315 #3.5*Loess.fit$s
  Loess.fit <- getLOESSfit(oswFiles, ref = "run2", eXp = "run0", maxFdrGlobal = 0.05, spanvalue = 0.1)
  outData <- getAlignedIndices(XICs.ref, XICs.eXp, Loess.fit, alignType = "hybrid",
                         adaptiveRT = adaptiveRT, normalization = "mean", simMeasure = "dotProductMasked",
                         goFactor = 0.125, geFactor = 40, cosAngleThresh = 0.3,
                         OverlapAlignment = TRUE,
                         dotProdThresh = 0.96, gapQuantile = 0.5, hardConstrain = FALSE,
                         samples4gradient = 100, objType = "light")
  expData1 <- list(c(4978.4, 4981.8, 4985.2), c(NA_real_, NA_real_, NA_real_))
  expData2 <- list(c(5569.0, 5572.4, 5575.8), c(5572.40, 5575.80, 5582.60))
  expect_equal(outData[[1]][1:3], expData1[[1]], tolerance = 1e-03)
  expect_equal(outData[[2]][1:3], expData1[[2]], tolerance = 1e-03)
  expect_equal(outData[[1]][174:176], expData2[[1]], tolerance = 1e-03)
  expect_equal(outData[[2]][174:176], expData2[[2]], tolerance = 1e-03)
  expect_identical(sapply(outData, length), c(176L, 176L))
})
