context("get_peak_chromatograms")

test_that("test_extractXIC_group", {
  mzmlName <- "../../data/example/mzml/hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML"
  mz <- mzR::openMSfile(mzmlName, backend = "pwiz")
  chromIndices <- c(37L, 38L, 39L, 40L, 41L, 42L)
  outData <- extractXIC_group(mz, chromIndices, SgolayFiltOrd = 4, SgolayFiltLen = 13)
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["run2"]][["14299_QFNNTDIVLLEDFQK/3"]]
  expect_identical(length(outData), 6L)
  expect_equal(outData[[2]][,1], XICs[[2]][,1], tolerance = 1e-04)
  expect_equal(outData[[2]][,2], XICs[[2]][,2], tolerance = 1e-04)
})

test_that("test_getXICs4AlignObj", {
  dataPath <- "../../data/example"
  runs <- c("run1" = "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
            "run0" =  "hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt")
  filenames <- data.frame("filename" = c("data/raw/hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt.mzML.gz",
                                         "data/raw/hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt.mzML.gz",
                                         "data/raw/hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.mzML.gz"),
                          "runs" = c("hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt",
                                     "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
                                     "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"),
                          row.names = c("run0", "run1", "run2"),
                          stringsAsFactors=FALSE)
  oswFiles <- getOswFiles(dataPath, filenames, maxFdrQuery = 0.05,
              analyteFDR = 0.01, oswMerged = TRUE, analyteInGroupLabel = FALSE)
  analytes <- "QFNNTDIVLLEDFQK_3"
  outData <- getXICs4AlignObj(dataPath, runs, oswFiles, analytes,
                   SgolayFiltOrd = 4, SgolayFiltLen = 13)
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR
  expect_identical(names(outData), c("run1", "run0"))
  expect_identical(names(outData[["run1"]]), "QFNNTDIVLLEDFQK_3")
  expect_equal(outData[["run0"]][["QFNNTDIVLLEDFQK_3"]], XICs[["run0"]][["14299_QFNNTDIVLLEDFQK/3"]], tolerance = 1e-03)
  expect_equal(outData[["run1"]][["QFNNTDIVLLEDFQK_3"]], XICs[["run1"]][["14299_QFNNTDIVLLEDFQK/3"]], tolerance = 1e-03)

  oswFiles <- getOswFiles(dataPath, filenames, maxFdrQuery = 0.05,
                          analyteFDR = 0.01, oswMerged = TRUE, analyteInGroupLabel = TRUE)
  outData <- getXICs4AlignObj(dataPath, runs, oswFiles, analytes = "14299_QFNNTDIVLLEDFQK/3",
                              SgolayFiltOrd = 4, SgolayFiltLen = 13)
  expect_equal(outData[["run0"]][["14299_QFNNTDIVLLEDFQK/3"]], XICs[["run0"]][["14299_QFNNTDIVLLEDFQK/3"]], tolerance = 1e-03)
  expect_equal(outData[["run1"]][["14299_QFNNTDIVLLEDFQK/3"]], XICs[["run1"]][["14299_QFNNTDIVLLEDFQK/3"]], tolerance = 1e-03)
})

test_that("test_getAlignObj", {
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR
  data(oswFiles_DIAlignR, package="DIAlignR")
  oswFiles <- oswFiles_DIAlignR
  XICs.ref <- XICs[["run1"]][["14299_QFNNTDIVLLEDFQK/3"]]
  XICs.eXp <- XICs[["run2"]][["14299_QFNNTDIVLLEDFQK/3"]]
  Loess.fit <- getLOESSfit(oswFiles, ref = "run1", eXp = "run2", maxFdrLoess = 0.05, spanvalue = 0.1)
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
  Loess.fit <- getLOESSfit(oswFiles, ref = "run2", eXp = "run0", maxFdrLoess = 0.05, spanvalue = 0.1)
  outData <- getMappedRT(refRT = 5238.35, XICs.ref, XICs.eXp, Loess.fit, alignType = "hybrid",
                         adaptiveRT = adaptiveRT, samplingTime = 3.414,
                         normalization = "mean", simMeasure = "dotProductMasked",
                         goFactor = 0.125, geFactor = 40, cosAngleThresh = 0.3,
                         OverlapAlignment = TRUE,
                         dotProdThresh = 0.96, gapQuantile = 0.5, hardConstrain = FALSE,
                         samples4gradient = 100)
  expect_equal(outData, 5237.8, tolerance = 1e-03)
})

test_that("test_getXICs", {
  runs <- c("hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt", "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
  outData <- getXICs(analytes = "QFNNTDIVLLEDFQK_3", runs = runs, dataPath = "../../data/example",
          maxFdrQuery = 1.0, SgolayFiltOrd = 4, SgolayFiltLen = 13, runType = "DIA_proteomics",
          oswMerged = TRUE, nameCutPattern = "(.*)(/)(.*)", analyteInGroupLabel = FALSE)
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR
  expect_equal(outData[["run0"]][["QFNNTDIVLLEDFQK_3"]], XICs[["run0"]][["14299_QFNNTDIVLLEDFQK/3"]], tolerance = 1e-03)
  expect_equal(outData[["run2"]][["QFNNTDIVLLEDFQK_3"]], XICs[["run2"]][["14299_QFNNTDIVLLEDFQK/3"]], tolerance = 1e-03)
  outData <- getXICs(analytes = "14299_QFNNTDIVLLEDFQK/3", runs = runs, dataPath = "../../data/example",
                     maxFdrQuery = 1.0, SgolayFiltOrd = 4, SgolayFiltLen = 13, runType = "DIA_proteomics",
                     oswMerged = TRUE, nameCutPattern = "(.*)(/)(.*)", analyteInGroupLabel = TRUE)
  expect_equal(outData[["run0"]][["14299_QFNNTDIVLLEDFQK/3"]], XICs[["run0"]][["14299_QFNNTDIVLLEDFQK/3"]], tolerance = 1e-03)
  expect_equal(outData[["run2"]][["14299_QFNNTDIVLLEDFQK/3"]], XICs[["run2"]][["14299_QFNNTDIVLLEDFQK/3"]], tolerance = 1e-03)
})
