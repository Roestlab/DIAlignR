context("Smoothing chromatograms")

test_that("test_smoothSingleXIC", {
  time <- seq(from = 3003.4, to = 3048, by = 3.4)
  y <- c(0.2050595, 0.8850070, 2.2068768, 3.7212677, 5.1652605, 5.8288915, 5.5446804,
             4.5671360, 3.3213154, 1.9485889, 0.9520709, 0.3294218, 0.2009581, 0.1420923)
  chrom <- cbind(time, y)

  ## Savitzky- Golay smoothing
  outData1 <- smoothSingleXIC(chrom, type = "sgolay", kernelLen = 9, polyOrd = 5)
  expOutput <- cbind(time, "intensity" = c(0.20636662, 0.88411087, 2.19019973, 3.76695006,
                          5.12687085, 5.77230554, 5.56200672, 4.5968725 , 3.2886408 , 1.97239146,
                          0.93076564, 0.34700936, 0.19229358, 0.14383756))
  expect_equal(outData1, expOutput, tolerance = 1e-03)

  ## Boxcar smoothing
  # kernelLen = 4 gives different output across 32-bit and 64-bit systems (Numerical instability).
  outData2 <- smoothSingleXIC(chrom, type = "boxcar", samplingTime = 3.4, kernelLen = 3.9)
  expOutput <- cbind(time, "intensity" = c(0.5450333, 1.0989811, 2.2710505, 3.6978017, 4.9051399,
                                                5.5129441, 5.3135693, 4.4777106, 3.2790134, 2.0739917,
                                                1.0766939, 0.4941503, 0.2241574, 0.1715252))
  expect_equal(outData2, expOutput, tolerance = 1e-03)

  ## Gaussian smoothing
  outData3 <- smoothSingleXIC(chrom, type = "gaussian", samplingTime = 3.4, kernelLen = 4)
  expOutput <- cbind(time, "intensity" = c(1.032401, 1.630418, 2.516603, 3.572851, 4.497432,
                                                4.975738, 4.860700, 4.215738, 3.249920, 2.219575,
                                                1.345037, 0.746991, 0.416589, 0.263005))
  expect_equal(outData3, expOutput, tolerance = 1e-03)

  ## Loess smoothing
  # Python code
  # from statsmodels.nonparametric.smoothers_lowess import lowess
  # import numpy as np
  # x = np.arange(3003.4, 3048, 3.4)
  # y = np.array([0.2050595, 0.8850070, 2.2068768, 3.7212677, 5.1652605, 5.8288915, 5.5446804,
  #               4.5671360, 3.3213154, 1.9485889, 0.9520709, 0.3294218, 0.2009581, 0.1420923])
  # z = lowess(y, x, frac= 3.9/len(x), it = 0)
  # z[:,1]
  # kernelLen = 4 gives different output across 32-bit and 64-bit systems (Numerical instability).
  outData4 <- smoothSingleXIC(chrom, type = "loess", kernelLen = 3.9, polyOrd = 1)
  expOutput <- cbind(time, "intensity" = c(0.2050595, 0.8850070, 2.2068768, 3.7212677, 5.1652605,
                                                5.8288915, 5.5446804, 4.5671360, 3.3213154, 1.9485889,
                                                0.9520709, 0.3294218, 0.2009581, 0.1420923))
  expect_equal(outData4, expOutput, tolerance = 1e-03)

  ## None
  outData5 <- smoothSingleXIC(chrom, type = "none")
  expOutput <- cbind(time, "intensity" = y)
  expect_equal(outData5, expOutput, tolerance = 1e-05)
})

test_that("test_smoothXICs", {
  time <- seq(from = 3003.4, to = 3048, by = 3.4)
  y <- c(0.2050595, 0.8850070, 2.2068768, 3.7212677, 5.1652605, 5.8288915, 5.5446804,
         4.5671360, 3.3213154, 1.9485889, 0.9520709, 0.3294218, 0.2009581, 0.1420923)
  chrom <- cbind(time, y)
  XICs <- list(chrom, chrom)

  ## Savitzky- Golay smoothing
  outData <- smoothXICs(XICs, type = "sgolay", kernelLen = 9, polyOrd = 5)
  expOutput <- cbind(time, "intensity1" = c(0.20636662, 0.88411087, 2.19019973, 3.76695006,
                                                5.12687085, 5.77230554, 5.56200672, 4.5968725 , 3.2886408 , 1.97239146,
                                                0.93076564, 0.34700936, 0.19229358, 0.14383756))
  expOutput <- list(expOutput, expOutput)
  colnames(expOutput[[2]]) <- c("time", "intensity2")
  expect_equal(outData, expOutput, tolerance = 1e-05)

  ## Boxcar smoothing
  # kernelLen = 4 gives different output across 32-bit and 64-bit systems (Numerical instability).
  outData <- smoothXICs(XICs, type = "boxcar", samplingTime = 3.4, kernelLen = 3.9)
  expOutput <- cbind(time, "intensity1" = c(0.5450333, 1.0989811, 2.2710505, 3.6978017, 4.9051399,
                                                 5.5129441, 5.3135693, 4.4777106, 3.2790134, 2.0739917,
                                                 1.0766939, 0.4941503, 0.2241574, 0.1715252))
  expOutput <- list(expOutput, expOutput)
  colnames(expOutput[[2]]) <- c("time", "intensity2")
  expect_equal(outData, expOutput, tolerance = 1e-05)

  ## Gaussian smoothing
  outData <- smoothXICs(XICs, type = "gaussian", samplingTime = 3.4, kernelLen = 4)
  expOutput <- cbind(time, "intensity1" = c(1.032401, 1.630418, 2.516603, 3.572851, 4.497432,
                                                4.975738, 4.860700, 4.215738, 3.249920, 2.219575,
                                                1.345037, 0.746991, 0.416589, 0.263005))
  expOutput <- list(expOutput, expOutput)
  colnames(expOutput[[2]]) <- c("time", "intensity2")
  expect_equal(outData, expOutput, tolerance = 1e-05)

  ## Loess smoothing
  # Python code
  # from statsmodels.nonparametric.smoothers_lowess import lowess
  # import numpy as np
  # x = np.arange(3003.4, 3048, 3.4)
  # y = np.array([0.2050595, 0.8850070, 2.2068768, 3.7212677, 5.1652605, 5.8288915, 5.5446804,
  #               4.5671360, 3.3213154, 1.9485889, 0.9520709, 0.3294218, 0.2009581, 0.1420923])
  # z = lowess(y, x, frac= 4/len(x), it = 0)
  # z[:,1]
  # kernelLen = 4 gives different output across 32-bit and 64-bit systems (Numerical instability).
  outData <- smoothXICs(XICs, type = "loess", kernelLen = 3.9, polyOrd = 1)
  expOutput <- cbind(time, "intensity1" = c(0.2050595, 0.8850070, 2.2068768, 3.7212677, 5.1652605,
                                                 5.8288915, 5.5446804, 4.5671360, 3.3213154, 1.9485889,
                                                 0.9520709, 0.3294218, 0.2009581, 0.1420923))
  expOutput <- list(expOutput, expOutput)
  colnames(expOutput[[2]]) <- c("time", "intensity2")
  expect_equal(outData, expOutput, tolerance = 1e-05)

  ## None
  outData <- smoothXICs(XICs, type = "none")
  expOutput <- cbind(time, "intensity1" = y)
  expOutput <- list(expOutput, expOutput)
  colnames(expOutput[[2]]) <- c("time", "intensity2")
  expect_equal(outData, expOutput, tolerance = 1e-05)
  })


test_that("test_smoothSingleXIC", {
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt"]][["4618"]]
  outData <- trimXICs(XICs, 0.04)[[2]]
  time <- c(5248.1, 5251.5, 5254.9, 5258.3, 5261.7, 5265.1, 5268.5, 5272.0, 5275.4)
  X27707 <- c(0.5315190, 0.3346700, 1.4568230, 0.5905989, 0.5905989, 0.5905989, 0.7086883,
              0.1378156, 0)
  expOutput <- data.frame(time, X27707)
  expect_equal(outData, expOutput, tolerance = 1e-05)
  })
