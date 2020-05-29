context("Align DIA runs")

test_that("test_alignTargetedRuns",{
  dataPath <- system.file("extdata", package = "DIAlignR")
  expect_warning(
    alignTargetedRuns(dataPath = dataPath,  outFile = "temp.csv", oswMerged = TRUE,
                               runs = NULL, runType = "DIA_Proteomics", context = "experiment-wide",
                      maxPeptideFdr = 1.00, maxFdrQuery = 0.05, XICfilter = "sgolay", polyOrd = 4, kernelLen = 9,
                    globalAlignment = "loess", globalAlignmentFdr = 0.01, globalAlignmentSpan = 0.1,
                    RSEdistFactor = 3.5, normalization = "mean", simMeasure = "dotProductMasked",
                    alignType = "hybrid", goFactor = 0.125, geFactor = 40,
                    cosAngleThresh = 0.3, OverlapAlignment = TRUE,
                    dotProdThresh = 0.96, gapQuantile = 0.5,
                    hardConstrain = FALSE, samples4gradient = 100,
                    analyteFDR = 1.0,
                    unalignedFDR = 0.01, alignedFDR = 0.05,
                    baselineType = "base_to_base", integrationType = "intensity_sum",
                    fitEMG = FALSE, recalIntensity = FALSE, fillMissing = TRUE, smoothPeakArea = FALSE)
  )
  outData <- read.table("temp.csv", stringsAsFactors = FALSE, sep = ",", header = TRUE)
  expData <- read.table("test.csv", stringsAsFactors = FALSE, sep = ",", header = TRUE)
  expect_identical(dim(outData), dim(expData))
  expect_identical(colnames(outData), colnames(expData))
  expect_identical(outData[["peptide"]], expData[["peptide"]])
  expect_identical(outData[["run"]], expData[["run"]])
  for(i in 3:13){
    expect_equal(outData[[i]], expData[[i]], tolerance = 1e-04)
  }
  file.remove("temp.csv")

  runs <- c("hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt",
            "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt")
  expect_warning(
    outData <- alignTargetedRuns(dataPath = dataPath,  outFile = "temp.csv", oswMerged = TRUE,
                               runs = runs, runType = "DIA_Proteomics", context = "experiment-wide",
                               maxPeptideFdr = 1.00, maxFdrQuery = 0.05, XICfilter = "sgolay", polyOrd = 4, kernelLen = 9,
                               globalAlignment = "loess", globalAlignmentFdr = 0.01, globalAlignmentSpan = 0.1,
                               RSEdistFactor = 3.5, normalization = "mean", simMeasure = "dotProductMasked",
                               alignType = "hybrid", goFactor = 0.125, geFactor = 40,
                               cosAngleThresh = 0.3, OverlapAlignment = TRUE,
                               dotProdThresh = 0.96, gapQuantile = 0.5,
                               hardConstrain = FALSE, samples4gradient = 100,
                               analyteFDR = 1.00,
                               unalignedFDR = 0.01, alignedFDR = 0.05,
                               baselineType = "base_to_base", integrationType = "intensity_sum",
                               fitEMG = FALSE, recalIntensity = FALSE, fillMissing = TRUE, smoothPeakArea = FALSE)
  )
  outData <- read.table("temp.csv", stringsAsFactors = FALSE, sep = ",", header = TRUE)
  expData <- read.table("test2.csv", stringsAsFactors = FALSE, sep = ",", header = TRUE)
  expect_identical(dim(outData), dim(expData))
  expect_identical(colnames(outData), colnames(expData))
  expect_identical(outData[["peptide"]], expData[["peptide"]])
  expect_identical(outData[["run"]], expData[["run"]])
  for(i in 3:13){
    expect_equal(outData[[i]], expData[[i]], tolerance = 1e-04)
  }
  file.remove("temp.csv")
})

test_that("test_getAlignObjs",{
  runs <- c("hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt", "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
  refRun <- "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"
  dataPath <- system.file("extdata", package = "DIAlignR")
  analytes <- c(32L, 898L, 4618L)
  expect_warning(
    outData <- getAlignObjs(analytes, runs, dataPath = dataPath, refRun = refRun,
               oswMerged = TRUE, runType = "DIA_Proteomics", maxFdrQuery = 0.05,
               analyteFDR = 0.01, XICfilter = "sgolay", polyOrd = 4,
               kernelLen = 13, globalAlignment = "loess",
               globalAlignmentFdr = 0.01, globalAlignmentSpan = 0.1,
               RSEdistFactor = 3.5, normalization = "mean",
               simMeasure = "dotProductMasked", alignType = "hybrid",
               goFactor = 0.125, geFactor = 40, cosAngleThresh = 0.3,
               OverlapAlignment = TRUE, dotProdThresh = 0.96, gapQuantile = 0.5,
               hardConstrain = FALSE, samples4gradient = 100, objType = "light")
    )
  expData <- testAlignObj()
  expect_equal(outData[[2]][["4618"]][["run1_run2"]][[1]], expData, tolerance = 1e-05)
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR
  expect_equal(outData[[2]][["4618"]][["run1_run2"]][["ref"]], XICs[["run1"]][["14299_QFNNTDIVLLEDFQK/3"]], tolerance = 1e-05)
  expect_equal(outData[[2]][["4618"]][["run1_run2"]][["eXp"]], XICs[["run2"]][["14299_QFNNTDIVLLEDFQK/3"]], tolerance = 1e-05)
  expData <- data.frame("leftWidth" = 5220.758, "RT" = 5238.35, "rightWidth" = 5261.723)
  expect_equal(as.data.frame(outData[[2]][["4618"]][["run1_run2"]][["peak"]]), expData, tolerance = 1e-05)
  expect_identical(outData[[2]][["32"]], NULL)
  expect_identical(outData[[2]][["898"]], NULL)
})
