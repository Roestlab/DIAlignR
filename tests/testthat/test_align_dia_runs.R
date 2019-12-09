context("Align DIA runs")

test_that("test_alignTargetedruns",{
  runs <- c("hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt", "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
  analyte <- "QFNNTDIVLLEDFQK_3"
  outData <- alignTargetedruns(dataPath= "../../data/example", alignType = "hybrid",
                               analyteInGroupLabel = FALSE, oswMerged = TRUE,
                               runs = runs, analytes = analyte, nameCutPattern = "(.*)(/)(.*)",
                               maxFdrQuery = 0.05, maxFdrLoess = 0.01, analyteFDR = 0.01,
                               spanvalue = 0.1, runType = "DIA_Proteomics",
                               normalization = "mean", simMeasure = "dotProductMasked",
                               XICfilter = "sgolay", SgolayFiltOrd = 4, SgolayFiltLen = 13,
                               goFactor = 0.125, geFactor = 40,
                               cosAngleThresh = 0.3, OverlapAlignment = TRUE,
                               dotProdThresh = 0.96, gapQuantile = 0.5,
                               hardConstrain = FALSE, samples4gradient = 100,
                               samplingTime = 3.4,  RSEdistFactor = 3.5, saveFiles = FALSE)
  expData <- matrix(c(310.01, 255.496), nrow = 1, ncol = 2, byrow = TRUE,
                    dimnames = list(c("QFNNTDIVLLEDFQK_3"), c("hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
                                                           "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")))
  expect_equal(outData, expData)

  outData <- alignTargetedruns(dataPath= "../../data/example", alignType = "hybrid",
                               analyteInGroupLabel = TRUE, oswMerged = TRUE,
                               runs = runs, analytes = "14299_QFNNTDIVLLEDFQK/3", nameCutPattern = "(.*)(/)(.*)",
                               maxFdrQuery = 0.05, maxFdrLoess = 0.01, analyteFDR = 0.01,
                               spanvalue = 0.1, runType = "DIA_Proteomics",
                               normalization = "mean", simMeasure = "dotProductMasked",
                               XICfilter = "sgolay", SgolayFiltOrd = 4, SgolayFiltLen = 13,
                               goFactor = 0.125, geFactor = 40,
                               cosAngleThresh = 0.3, OverlapAlignment = TRUE,
                               dotProdThresh = 0.96, gapQuantile = 0.5,
                               hardConstrain = FALSE, samples4gradient = 100,
                               samplingTime = 3.4,  RSEdistFactor = 3.5, saveFiles = FALSE)
  expData <- matrix(c(310.01, 255.496), nrow = 1, ncol = 2, byrow = TRUE,
                    dimnames = list(c("14299_QFNNTDIVLLEDFQK/3"), c("hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
                                                              "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")))
  expect_equal(outData, expData)
})

test_that("test_getAlignObjs",{
  runs <- c("hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt", "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
  refRun <- "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"

  analyte <- "QFNNTDIVLLEDFQK_3"
  outData <- getAlignObjs(analytes = analyte, runs = runs, dataPath = "../../data/example",
               alignType = "hybrid", runType = "DIA_Proteomics", refRun = refRun,
               analyteInGroupLabel = FALSE, oswMerged = TRUE, nameCutPattern = "(.*)(/)(.*)",
               maxFdrQuery = 0.05, maxFdrLoess = 0.01, analyteFDR = 1.00, spanvalue = 0.1,
               normalization = "mean", simMeasure = "dotProductMasked",
               SgolayFiltOrd = 4, SgolayFiltLen = 13,
               goFactor = 0.125, geFactor = 40,
               cosAngleThresh = 0.3, OverlapAlignment = TRUE,
               dotProdThresh = 0.96, gapQuantile = 0.5,
               hardConstrain = FALSE, samples4gradient = 100,
               samplingTime = 3.4,  RSEdistFactor = 3.5, objType = "light")
  expData <- testAlignObj()
  expect_equal(outData[[analyte]][[1]], expData, tolerance = 1e-05)
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR
  expect_equal(outData[[analyte]][[refRun]], XICs[["run1"]][["14299_QFNNTDIVLLEDFQK/3"]], tolerance = 1e-05)
  expect_equal(outData[[analyte]][[setdiff(runs, refRun)]], XICs[["run2"]][["14299_QFNNTDIVLLEDFQK/3"]], tolerance = 1e-05)
  expData <- data.frame("leftWidth" = 5220.758, "RT" = 5238.35, "rightWidth" = 5261.723)
  expect_equal(as.data.frame(outData[[analyte]][[4]]), expData, tolerance = 1e-05)

  analyte <- "14299_QFNNTDIVLLEDFQK/3"
  outData <- getAlignObjs(analytes = analyte, runs = runs, dataPath = "../../data/example",
                          alignType = "hybrid", runType = "DIA_Proteomics", refRun = refRun,
                          analyteInGroupLabel = TRUE, oswMerged = TRUE, nameCutPattern = "(.*)(/)(.*)",
                          maxFdrQuery = 0.05, maxFdrLoess = 0.01, analyteFDR = 1.00, spanvalue = 0.1,
                          normalization = "mean", simMeasure = "dotProductMasked",
                          SgolayFiltOrd = 4, SgolayFiltLen = 13,
                          goFactor = 0.125, geFactor = 40,
                          cosAngleThresh = 0.3, OverlapAlignment = TRUE,
                          dotProdThresh = 0.96, gapQuantile = 0.5,
                          hardConstrain = FALSE, samples4gradient = 100,
                          samplingTime = 3.4,  RSEdistFactor = 3.5, objType = "light")
  expData <- testAlignObj(analyteInGroupLabel = TRUE)
  expect_equal(outData[[analyte]][[1]], expData, tolerance = 1e-05)
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR
  expect_equal(outData[[analyte]][[refRun]], XICs[["run1"]][["14299_QFNNTDIVLLEDFQK/3"]], tolerance = 1e-05)
  expect_equal(outData[[analyte]][[setdiff(runs, refRun)]], XICs[["run2"]][["14299_QFNNTDIVLLEDFQK/3"]], tolerance = 1e-05)
  expData <- data.frame("leftWidth" = 5220.758, "RT" = 5238.35, "rightWidth" = 5261.723)
  expect_equal(as.data.frame(outData[[analyte]][[4]]), expData, tolerance = 1e-05)
})
