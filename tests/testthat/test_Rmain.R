context("Rmain.cpp")

test_that("test_getSeqSimMatCpp",{
  seq1 <- "GCAT"
  seq2 <- "CAGTG"
  outData <- getSeqSimMatCpp(seq1, seq2, match = 10, misMatch = -2)
  expData <- matrix(c(-2, 10, -2, -2, -2, -2, 10, -2, 10, -2, -2,
                      -2, -2, -2, -2, 10, 10, -2, -2, -2), 4, 5, byrow = FALSE)
  expect_equal(outData, expData)
})

test_that("test_getChromSimMatCpp",{
  r1 <- list(c(1.0,3.0,2.0,4.0), c(0.0,0.0,0.0,1.0), c(4.0,4.0,4.0,5.0))
  r2 <- list(c(1.4,2.0,1.5,4.0), c(0.0,0.5,0.0,0.0), c(2.0,3.0,4.0,0.9))
  outData <- round(getChromSimMatCpp(r1, r2, "L2", "dotProductMasked"), 3)
  expData <- matrix(c(0.125, 0.162, 0.144, 0.208, 0.186, 0.240,0.213, 0.313, 0.233,
           0.273, 0.253, 0.346, 0.101, 0.208, 0.154, 0.273), 4, 4, byrow = FALSE)
  expect_equal(outData, expData)

  outData <- round(getChromSimMatCpp(r1, r2, "L2", "dotProduct"), 3)
  expData <- matrix(c(0.125, 0.162, 0.144, 0.208, 0.186,0.240, 0.213, 0.313,
                      0.233, 0.273, 0.253, 0.346, 0.101, 0.208, 0.154, 0.273), 4, 4, byrow = FALSE)
  expect_equal(outData, expData)

  outData <- round(getChromSimMatCpp(r1, r2, "L2", "cosineAngle"), 3)
  expData <- matrix(c(0.934, 0.999, 0.989, 0.986, 0.933, 0.989, 0.983, 0.996,
                      0.994, 0.960, 0.995, 0.939, 0.450, 0.761, 0.633, 0.772), 4, 4, byrow = FALSE)
  expect_equal(outData, expData)

  outData <- round(getChromSimMatCpp(r1, r2, "L2", "cosine2Angle"), 3)
  expData <- matrix(c(0.744, 0.998, 0.957, 0.944, 0.740, 0.956, 0.932, 0.985,
                      0.974, 0.842, 0.978, 0.764, -0.596, 0.158, -0.200, 0.190), 4, 4, byrow = FALSE)

  outData <- round(getChromSimMatCpp(r1, r2, "mean", "euclideanDist"), 3)
  expData <- matrix(c(0.608, 0.614, 0.680, 0.434, 0.530, 0.742, 0.659, 0.641,
                      0.520, 0.541, 0.563, 0.511, 0.298,0.375, 0.334, 0.355), 4, 4, byrow = FALSE)
  expect_equal(outData, expData)

  outData <- round(getChromSimMatCpp(r1, r2, "L2", "covariance"), 3)
  expData <- matrix(c(0.025, 0.028, 0.027, 0.028, 0.032, 0.034, 0.033, 0.034,
                      0.055, 0.051, 0.053, 0.051, -0.004, 0.028, 0.012, 0.028), 4, 4, byrow = FALSE)
  expect_equal(outData, expData)

  outData <- round(getChromSimMatCpp(r1, r2, "L2", "correlation"), 3)
  expData <- matrix(c(0.874, 0.999, 0.974, 0.999, 0.923, 0.986, 0.993, 0.986, 0.991, 0.911,
         0.990, 0.911, -0.065, 0.477, 0.214, 0.477), 4, 4, byrow = FALSE)
  expect_equal(outData, expData)
})

test_that("test_getGlobalAlignMaskCpp",{
  tA <- c(3353.2, 3356.6, 3360.0, 3363.5)
  tB <- c(3325.9, 3329.3, 3332.7, 3336.1)
  B1p <- 3325.751; B2p <- 3336.119
  noBeef <- 1
  mask <- getGlobalAlignMaskCpp(tA, tB, B1p, B2p, noBeef, FALSE)
  outData <- round(mask, 3)
  expData <-matrix(c(0.000, 0.000, 0.707, 1.414, 0.000, 0.000, 0.000, 0.707, 0.707, 0.000,
             0.000, 0.000, 1.414, 0.707, 0.000, 0.000), 4, 4, byrow = FALSE)
  expect_equal(outData, expData)
})

test_that("test_constrainSimCpp",{
  sim <- matrix(c(-2, 10, -2, -2, -2, -2, 10, -2, 10, -2, -2, -2, -2,
                  -2, -2, 10, 10, -2,-2, -2), 4, 5, byrow = FALSE)
  MASK <- matrix(c(0.000, 0.000, 0.707, 1.414, 0.000, 0.000, 0.000, 0.707, 0.707,
                   0.000, 0.000, 0.000, 1.414, 0.707, 0, 0, 2.121, 1.414, 0, 0), 4, 5, byrow = FALSE)
  outData <- constrainSimCpp(sim, MASK, 10)
  expData <-matrix(c(-2, 10, -3.414, -4.828, -2, -2, 10, -3.414, 8.586, -2, -2, -2, -4.828,
            -3.414, -2, 10, 5.758, -4.828, -2, -2), 4, 5, byrow = FALSE)
  expect_equal(outData, expData)
})

test_that("test_getBaseGapPenaltyCpp",{
  sim <- matrix(c(-12, 1.0, 12, -2.3, -2, -2, 1.07, -2, 1.80,
                  2, 22, 42, -2, -1.5, -2, 10), 4, 4, byrow = FALSE)
  expect_equal(getBaseGapPenaltyCpp(sim, "dotProductMasked", 0.5), -0.25)
})

test_that("test_areaIntegrator",{
  time <- c( 2.23095,2.239716667,2.248866667,2.25765,2.266416667,
             2.275566667,2.2847,2.293833333,2.304066667,2.315033333,2.325983333,2.336566667,
             2.3468,2.357016667,2.367283333,2.377183333,2.387083333,2.39735,2.40725,2.4175,
             2.4274,2.4373,2.44755,2.45745,2.4677,2.477966667,2.488216667,2.498516667,2.5084,
             2.5183,2.5282,2.538466667,2.548366667,2.558266667,2.568516667,2.578783333,
             2.588683333,2.59895,2.6092,2.619466667,2.630066667,2.64065,2.65125,2.662116667,
             2.672716667,2.6833,2.6939,2.7045,2.715083333,2.725683333,2.736266667,2.746866667,
             2.757833333,2.768416667,2.779016667,2.789616667,2.8002,2.810116667,2.820033333,
             2.830316667,2.840216667,2.849766667,2.859316667,2.868866667,2.878783333,2.888683333,
             2.898233333,2.907783333,2.916033333,2.924266667,2.93215,2.940383333,2.947933333,
             2.955816667,2.964066667,2.97195,2.979833333,2.987716667,2.995616667,3.003516667,
             3.011416667,3.01895,3.026833333,3.034366667,3.042266667,3.0498,3.05735,3.065233333,
             3.073133333,3.080666667,3.0882,3.095733333,3.103633333,3.111533333,3.119066667,
             3.126966667,3.134866667,3.14275,3.15065,3.15855,3.166433333,3.174333333,3.182233333,
             3.190133333,3.198016667,3.205916667,3.213166667
  )

  intensity <- c(
    1447,2139,1699,755,1258,1070,944,1258,1573,1636,
    1762,1447,1133,1321,1762,1133,1447,2391,692,1636,2957,1321,1573,1196,1258,881,
    1384,2076,1133,1699,1384,692,1636,1133,1573,1825,1510,2391,4342,10382,17618,
    51093,153970,368094,632114,869730,962547,966489,845055,558746,417676,270942,
    184865,101619,59776,44863,31587,24036,20450,20324,11074,9879,10508,7928,7110,
    6733,6481,5726,6921,6670,5537,4971,4719,4782,5097,5789,4279,5411,4530,3524,
    2139,3335,3083,4342,4279,3083,3649,4216,4216,3964,2957,2202,2391,2643,3524,
    2328,2202,3649,2706,3020,3335,2580,2328,2894,3146,2769,2517
  )
  left <- 2.472833334
  right <- 3.022891666
  outData <- areaIntegrator(list(time), list(intensity), left, right, integrationType = "intensity_sum",
                                baselineType = "base_to_base", fitEMG = FALSE)
  expect_equal(outData, 6645331.33866)

  outData <- areaIntegrator(list(time, time), list(intensity, intensity), left, right, integrationType = "trapezoid",
                                baselineType = "vertical_division_min", fitEMG = FALSE)
  expect_equal(outData, 2*71063.59368, tolerance = 0.01)
})

test_that("test_alignChromatogramsCpp",{
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR
  data(oswFiles_DIAlignR, package="DIAlignR")
  oswFiles <- oswFiles_DIAlignR
  Loess.fit <- getLOESSfit(oswFiles, ref = "run1", eXp = "run2", maxFdrGlobal = 0.05, spanvalue = 0.1)
  XICs.ref <- XICs[["run1"]][["14299_QFNNTDIVLLEDFQK/3"]]
  XICs.ref <- smoothXICs(XICs.ref, type = "sgolay", kernelLen = 13, polyOrd = 4)
  XICs.eXp <- XICs[["run2"]][["14299_QFNNTDIVLLEDFQK/3"]]
  XICs.eXp <- smoothXICs(XICs.eXp, type = "sgolay", kernelLen = 13, polyOrd = 4)
  tVec.ref <- XICs.ref[[1]][["time"]] # Extracting time component
  tVec.eXp <- XICs.eXp[[1]][["time"]] # Extracting time component
  B1p <- predict(Loess.fit, tVec.ref[1])
  B2p <- predict(Loess.fit, tVec.ref[length(tVec.ref)])
  noBeef <- 38.6594179136227/3.414
  l1 <- lapply(XICs.ref, `[[`, 2)
  l2 <- lapply(XICs.eXp, `[[`, 2)
  outData <- alignChromatogramsCpp(l1, l2, alignType = "hybrid",
                                   tA = tVec.ref, tB = tVec.eXp, normalization = "mean", simType = "dotProductMasked",
                                   B1p = B1p, B2p = B2p, noBeef = noBeef,
                         goFactor = 0.125, geFactor = 40,
                         cosAngleThresh = 0.3, OverlapAlignment = TRUE,
                         dotProdThresh = 0.96, gapQuantile = 0.5, hardConstrain = FALSE,
                         samples4gradient = 100, objType = "light")
  expData <- testAlignObj()
  expect_equal(outData, expData, tolerance = 1e-03)
})

test_that("test_doAlignmentCpp",{
  s <- getSeqSimMatCpp(seq1 = "GCAT", seq2 = "CAGTG", match = 10, misMatch = -2)
  outData <- doAlignmentCpp(s, 22, FALSE)
  expData <- c(-2, -4, -6, 4, -18)
  expect_equal(outData@score, expData)
  outData <- doAlignmentCpp(s, 22, TRUE)
  expData <- c(0, 10, 20, 18, 18, 18)
  expect_equal(outData@score, expData)

  s <- getSeqSimMatCpp(seq1 = "TTTC", seq2 = "TGC", match = 1, misMatch = -1)
  outData <- doAlignmentCpp(s, 2, FALSE)
  expData <- matrix(data = c(1,1,1,1,1,1,1,1,1,2,1,2,1,
                             3,3,1,1,3,6,3), nrow = 5, ncol =4, byrow = TRUE)
  expect_equal(outData@optionalPaths, expData)
  expData <- matrix(data = c(0,-2,-4,-6,-2,-7,-22,-45,-4,-20,-72,-184,-6,-41,-178,-547,-8,-72,-366,-1274), nrow = 5, ncol =4, byrow = TRUE)
  expect_equal(outData@M_forw, expData)
})

test_that("test_doAffineAlignmentCpp",{
  Match <- 10
  MisMatch <- -2
  s <- getSeqSimMatCpp(seq1 = "GCAT", seq2=  "CAGTG", Match, MisMatch)
  outData <- doAffineAlignmentCpp(s, 22, 7, FALSE)
  expData <- c(-2, -4, -6, 4, -18)
  expect_equal(outData@score, expData)
  outData <- doAffineAlignmentCpp(s, 22, 7, TRUE)
  expData <- c(0, 10, 20, 18, 18, 18)
  expect_equal(outData@score, expData)

  s <- getSeqSimMatCpp(seq1 = "CAT", seq2 = "CAGTG", Match, MisMatch)
  outData <- doAffineAlignmentCpp(s, 22, 7, FALSE)
  expData <- c(10, 20, -2, -9, -11)
  expect_equal(outData@score, expData)
  outData <- doAffineAlignmentCpp(s, 22, 7, TRUE)
  expData <- c(10, 20, 18, 18, 18)
  expect_equal(outData@score, expData)
})
