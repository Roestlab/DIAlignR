context("get_peak_chromatograms")

testAlignObj <- function(){
  AlignObj <- new("AffineAlignObjLight",
                indexA_aligned = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,0,52,0,53,0,54,0,55,0,56,0,57,0,58,0,59,0,60,0,61,0,62,0,63,0,64,0,65,0,66,0,67,0,68,0,69,0,70,0,71,0,72,0,73,0,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,0,162,0,163,164,165,166,0,167,0,168,169,170,171,172,173,174,175,176),
                indexB_aligned = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,2,3,4,5,6,7,8,9,0,10,0,11,0,12,0,13,0,14,0,15,0,16,0,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,0,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176),
                score = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8.258636,21.70066,37.77944,56.28491,77.49291,102.7556,131.1216,159.8879,186.0392,185.7486,199.3454,199.0548,208.9874,208.6968,224.2224,223.9318,235.3893,235.0987,239.616,239.3254,239.8545,239.5639,239.7825,239.4919,264.0039,345.6801,518.224,760.9532,1002.997,1171.843,1255.645,1283.27,1289.328,1289.037,1289.66,1289.37,1289.694,1289.404,1289.829,1289.538,1290.128,1289.838,1290.421,1290.13,1290.703,1290.413,1291.675,1291.385,1292.812,1292.522,1293.994,1293.704,1295.055,1294.764,1295.694,1295.403,1296.342,1296.051,1297.042,1296.752,1297.599,1297.308,1298.15,1297.859,1298.766,1298.475,1299.462,1299.172,1300.227,1299.936,1300.543,1300.252,1300.388,1300.097,1300.442,1300.151,1301.528,1301.238,1350.073,1535.732,1925.536,2471.722,2995.898,3355.553,3518.7,3568.132,3579.685,3585.838,3592.471,3599.173,3606.964,3614.774,3621.711,3627.712,3632.502,3635.622,3637.255,3638.218,3638.903,3639.589,3639.299,3640.431,3642.042,3643.971,3645.901,3648.116,3650.643,3653.25,3655.58,3657.816,3660.221,3662.642,3664.741,3666.717,3668.866,3670.707,3672.353,3674.09,3676.172,3678.516,3680.837,3683.553,3687.36,3692.948,3700.384,3710.218,3721.758,3735.07,3749.023,3761.307,3771.312,3779.496,3786.653,3793.472,3800.401,3810.434,3826.845,3850.323,3881.263,3915.979,3945.088,3962.835,3970.613,3974.154,3976.556,3979.059,3981.622,3983.965,3986.197,3987.975,3989.738,3990.982,3992.252,3993.539,3994.95,3996.329,3997.578,3998.949,4000.323,4001.774,4003.358,4005.047,4006.826,4008.501,4010.002,4011.12,4010.829,4011.757,4011.466,4012.571,4013.726,4014.719,4015.612,4015.322,4016.153,4015.862,4016.773,4017.729,4018.441,4019.06,4019.624,4020.201,4020.895,4021.846,4023.249))
  AlignObj
}

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
