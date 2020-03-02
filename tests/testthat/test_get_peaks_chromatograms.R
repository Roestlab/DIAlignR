context("get_peak_chromatograms")

test_that("test_extractXIC_group", {
  mzmlName <- file.path(system.file("extdata", package = "DIAlignR"), "mzml", "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML")
  mz <- mzR::openMSfile(filename = mzmlName, backend = "pwiz")
  chromIndices <- c(37L, 38L, 39L, 40L, 41L, 42L)
  outData <- extractXIC_group(mz, chromIndices)
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["run2"]][["14299_QFNNTDIVLLEDFQK/3"]]
  expect_identical(length(outData), 6L)
  expect_equal(outData[[2]][,1], XICs[[2]][,1], tolerance = 1e-04)
  expect_equal(outData[[1]][,2], XICs[[1]][,2], tolerance = 1e-04)
})

test_that("test_getXICs4AlignObj", {
  dataPath <- system.file("extdata", package = "DIAlignR")
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
  oswFiles <- getOswFiles(dataPath = dataPath, filenames, maxFdrQuery = 0.05,
              analyteFDR = 0.01, oswMerged = TRUE, analyteInGroupLabel = FALSE)
  analytes <- "QFNNTDIVLLEDFQK_3"
  outData <- getXICs4AlignObj(dataPath = dataPath, runs, oswFiles, analytes,
                   SgolayFiltOrd = 4, SgolayFiltLen = 13)
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR
  expect_identical(names(outData), c("run1", "run0"))
  expect_identical(names(outData[["run1"]]), "QFNNTDIVLLEDFQK_3")
  expect_equal(outData[["run0"]][["QFNNTDIVLLEDFQK_3"]], XICs[["run0"]][["14299_QFNNTDIVLLEDFQK/3"]], tolerance = 1e-03)
  expect_equal(outData[["run1"]][["QFNNTDIVLLEDFQK_3"]], XICs[["run1"]][["14299_QFNNTDIVLLEDFQK/3"]], tolerance = 1e-03)

  oswFiles <- getOswFiles(dataPath = dataPath, filenames, maxFdrQuery = 0.05,
                          analyteFDR = 0.01, oswMerged = TRUE, analyteInGroupLabel = TRUE)
  outData <- getXICs4AlignObj(dataPath = dataPath, runs, oswFiles, analytes = "14299_QFNNTDIVLLEDFQK/3",
                              SgolayFiltOrd = 4, SgolayFiltLen = 13)
  expect_equal(outData[["run0"]][["14299_QFNNTDIVLLEDFQK/3"]], XICs[["run0"]][["14299_QFNNTDIVLLEDFQK/3"]], tolerance = 1e-03)
  expect_equal(outData[["run1"]][["14299_QFNNTDIVLLEDFQK/3"]], XICs[["run1"]][["14299_QFNNTDIVLLEDFQK/3"]], tolerance = 1e-03)
})

test_that("test_getXICs", {
  dataPath <- system.file("extdata", package = "DIAlignR")
  runs <- c("hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt", "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
  outData <- getXICs(analytes = "QFNNTDIVLLEDFQK_3", runs = runs, dataPath = dataPath,
          maxFdrQuery = 1.0, SgolayFiltOrd = 4, SgolayFiltLen = 13, runType = "DIA_proteomics",
          oswMerged = TRUE, nameCutPattern = "(.*)(/)(.*)", analyteInGroupLabel = FALSE)
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR
  expect_equal(outData[["hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt"]][["QFNNTDIVLLEDFQK_3"]],
               XICs[["run0"]][["14299_QFNNTDIVLLEDFQK/3"]], tolerance = 1e-03)
  expect_equal(outData[["hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"]][["QFNNTDIVLLEDFQK_3"]],
               XICs[["run2"]][["14299_QFNNTDIVLLEDFQK/3"]], tolerance = 1e-03)
  outData <- getXICs(analytes = "14299_QFNNTDIVLLEDFQK/3", runs = runs, dataPath = dataPath,
                     maxFdrQuery = 1.0, SgolayFiltOrd = 4, SgolayFiltLen = 13, runType = "DIA_proteomics",
                     oswMerged = TRUE, nameCutPattern = "(.*)(/)(.*)", analyteInGroupLabel = TRUE)
  expect_equal(outData[["hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt"]][["14299_QFNNTDIVLLEDFQK/3"]], XICs[["run0"]][["14299_QFNNTDIVLLEDFQK/3"]], tolerance = 1e-03)
  expect_equal(outData[["hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"]][["14299_QFNNTDIVLLEDFQK/3"]], XICs[["run2"]][["14299_QFNNTDIVLLEDFQK/3"]], tolerance = 1e-03)
})
