context("Pairwise analyte peak correspondence.")

test_that("test_getAlignObj", {
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR
  data(oswFiles_DIAlignR, package="DIAlignR")
  oswFiles <- oswFiles_DIAlignR
  XICs.ref <- XICs[["hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]]
  XICs.ref <- smoothXICs(XICs.ref, type = "sgolay", kernelLen = 13, polyOrd = 4)
  XICs.eXp <- XICs[["hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]]
  XICs.eXp <- smoothXICs(XICs.eXp, type = "sgolay", kernelLen = 13, polyOrd = 4)
  RUNS_RT <- getRTdf(oswFiles, ref = "run1", eXp = "run2", maxFdrGlobal = 0.05)
  globalFit <- loess(RT.eXp ~ RT.ref, data = RUNS_RT, span = 0.1, control=loess.control(surface="direct"))
  outData <- getAlignObj(XICs.ref, XICs.eXp, globalFit, alignType = "hybrid", adaptiveRT = 38.6594179136227,
                         normalization = "mean", simType = "dotProductMasked", goFactor = 0.125, geFactor = 40,
                         cosAngleThresh = 0.3, OverlapAlignment = TRUE, dotProdThresh = 0.96,
                         gapQuantile = 0.5, kerLen = 9, hardConstrain = FALSE,
                         samples4gradient = 100, objType = "light")
  expData <- testAlignObj()
  expect_equal(outData, expData, tolerance = 1e-03)
})

test_that("test_getMappedRT", {
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR
  data(oswFiles_DIAlignR, package="DIAlignR")
  oswFiles <- oswFiles_DIAlignR
  XICs.ref <- XICs[["hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]]
  XICs.eXp <- XICs[["hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]]
  adaptiveRT <- 77.82315 #3.5*Loess.fit$s
  RUNS_RT <- getRTdf(oswFiles, ref = "run2", eXp = "run0", maxFdrGlobal = 0.05)
  Loess.fit <- loess(RT.eXp ~ RT.ref, data = RUNS_RT, span = 0.1, control=loess.control(surface="direct"))
  outData <- getMappedRT(refRT = 5238.35, XICs.ref, XICs.eXp, Loess.fit, alignType = "hybrid",
                         adaptiveRT = adaptiveRT, normalization = "mean", simMeasure = "dotProductMasked",
                         goFactor = 0.125, geFactor = 40, cosAngleThresh = 0.3,
                         OverlapAlignment = TRUE, dotProdThresh = 0.96, gapQuantile = 0.5,
                         kerLen = 9, hardConstrain = FALSE, samples4gradient = 100)
  expect_equal(outData, 5237.8, tolerance = 1e-03)
})

test_that("test_getAlignedTimes", {
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR
  data(oswFiles_DIAlignR, package="DIAlignR")
  oswFiles <- oswFiles_DIAlignR
  XICs.ref <- XICs[["hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]]
  XICs.eXp <- XICs[["hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]]
  adaptiveRT <- 77.82315 #3.5*Loess.fit$s
  RUNS_RT <- getRTdf(oswFiles, ref = "run2", eXp = "run0", maxFdrGlobal = 0.05)
  Loess.fit <- loess(RT.eXp ~ RT.ref, data = RUNS_RT, span = 0.1, control=loess.control(surface="direct"))
  outData <- getAlignedTimes(XICs.ref, XICs.eXp, Loess.fit, alignType = "hybrid",
                         adaptiveRT = adaptiveRT, normalization = "mean", simMeasure = "dotProductMasked",
                         goFactor = 0.125, geFactor = 40, cosAngleThresh = 0.3,
                         OverlapAlignment = TRUE, dotProdThresh = 0.96, gapQuantile = 0.5,
                         kerLen = 9, hardConstrain = FALSE, samples4gradient = 100, objType = "light")
  expData1 <- list(c(4978.4, 4981.8, 4985.2), c(NA_real_, NA_real_, NA_real_))
  expData2 <- list(c(5569.0, 5572.4, 5575.8), c(5572.40, 5575.80, 5582.60))
  expect_equal(outData[[1]][1:3], expData1[[1]], tolerance = 1e-03)
  expect_equal(outData[[2]][1:3], expData1[[2]], tolerance = 1e-03)
  expect_equal(outData[[1]][174:176], expData2[[1]], tolerance = 1e-03)
  expect_equal(outData[[2]][174:176], expData2[[2]], tolerance = 1e-03)
  expect_identical(sapply(outData, length), c(176L, 176L))
})

test_that("test_getAlignedIndices", {
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR
  data(oswFiles_DIAlignR, package="DIAlignR")
  XICs.ref <- XICs[["hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]]
  XICs.eXp <- XICs[["hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]]
  adaptiveRT <- 77.82315 #3.5*Loess.fit$s
  RUNS_RT <- getRTdf(oswFiles_DIAlignR, ref = "run2", eXp = "run0", maxFdrGlobal = 0.05)
  Loess.fit <- loess(RT.eXp ~ RT.ref, data = RUNS_RT, span = 0.1, control=loess.control(surface="direct"))
  outData <- getAlignedIndices(XICs.ref, XICs.eXp, Loess.fit, alignType = "hybrid",
                             adaptiveRT = adaptiveRT, normalization = "mean", simMeasure = "dotProductMasked",
                             goFactor = 0.125, geFactor = 40, cosAngleThresh = 0.3,
                             OverlapAlignment = TRUE, dotProdThresh = 0.96, gapQuantile = 0.5,
                             kerLen = 9, hardConstrain = FALSE, samples4gradient = 100, objType = "light")
  expect_equal(outData[1:3, "indexAligned.ref"], c(1L, 2L, 3L))
  expect_equal(outData[1:3, "indexAligned.eXp"], rep(NA_integer_, 3))
  expect_equal(outData[203:205, "indexAligned.ref"], c(NA_integer_, 176L, NA_integer_))
  expect_equal(outData[203:205, "indexAligned.eXp"], c(174L, 175L, 176L))
  expect_equal(outData[203:205, "score"], c(4167.9972, 4169.8598, 4169.8598), tolerance = 1e-03)
  expect_identical(dim(outData), c(205L, 3L))
})


test_that("test_getAlignedTimes", {
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  data(oswFiles_DIAlignR, package="DIAlignR")
  params <- paramsDIAlignR()
  params$globalAlignment <- "loess"
  params$kernelLen <- 0L # No filter
  run1 <- "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"
  run2 <- "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"
  XICs.ref <- lapply(XIC_QFNNTDIVLLEDFQK_3_DIAlignR[[run1]][["4618"]], as.matrix)
  XICs.eXp <- lapply(XIC_QFNNTDIVLLEDFQK_3_DIAlignR[[run2]][["4618"]], as.matrix)
  adaptiveRT <- 77.82315
  fit <- getLOESSfit(oswFiles_DIAlignR, ref = "run2", eXp = "run0", maxFdrGlobal = 0.05, spanvalue = 0.1)
  fit <- extractFit(fit, params$globalAlignment)
  outData <- getAlignedTimesFast(XICs.ref, XICs.eXp, fit, adaptiveRT, params)
  expect_equal(outData[1:3,1], c(4978.4, 4981.8, 4985.2), tolerance = 1e-03)
  expect_equal(outData[1:3,2], c(NA_real_, NA_real_, NA_real_), tolerance = 1e-03)
  expect_equal(outData[174:176,2], c(5569.0, 5572.4, 5575.8), tolerance = 1e-03)
  expect_equal(outData[174:176,2], c(5572.40, 5575.80, 5582.60), tolerance = 1e-03)
  expect_identical(dim(outData), c(176L, 2L))
})

