context("Impute chromatogram")

test_that("test_sgolayFill", {
  time <- seq(from = 3003.4, to = 3048, by = 3.4)
  y <- c(0.2050595, 0.8850070, 2.2068768, 3.7212677, 5.1652605, 5.8288915, 5.5446804,
         4.5671360, 3.3213154, 1.9485889, 0.9520709, 0.3294218, 0.2009581, 0.1420923)
  chrom <- data.frame(time, y)
  chrom$y[6] <- NA_real_
  outData1 <- sgolayFill(chrom, polyOrd = 3, kernelLen = 9)
  expOutput <- data.frame(time, "intensity" = c(0.2050595, 0.8850070, 2.2068768, 3.7212677, 5.1652605,
                                                5.3457775, 5.5446804, 4.5671360, 3.3213154, 1.9485889,
                                                0.9520709, 0.3294218, 0.2009581, 0.1420923))
  expect_equal(outData1, expOutput, tolerance = 1e-03)

  chrom$y[1] <- NA_real_
  outData2 <- sgolayFill(chrom, polyOrd = 3, kernelLen = 9)
  expOutput <- data.frame(time, "intensity" = c(NA_real_, 0.8850070, 2.2068768, 3.7212677, 5.1652605,
                                                5.3457775, 5.5446804, 4.5671360, 3.3213154, 1.9485889,
                                                0.9520709, 0.3294218, 0.2009581, 0.1420923))
  expect_equal(outData2, expOutput, tolerance = 1e-03)

  chrom$y[6] <-5.8288915
  outData3 <- sgolayFill(chrom, polyOrd = 4, kernelLen = 9)
  expOutput <- data.frame(time, "intensity" = c(1.0760958, 0.8850070, 2.2068768, 3.7212677, 5.1652605,
                                                5.8288915, 5.5446804, 4.5671360, 3.3213154, 1.9485889,
                                                0.9520709, 0.3294218, 0.2009581, 0.1420923))
  expect_equal(outData3, expOutput, tolerance = 1e-03)
})

test_that("test_splineFill", {
  time <- seq(from = 3003.4, to = 3048, by = 3.4)
  y <- c(0.2050595, 0.8850070, 2.2068768, 3.7212677, 5.1652605, 5.8288915, 5.5446804,
         4.5671360, 3.3213154, 1.9485889, 0.9520709, 0.3294218, 0.2009581, 0.1420923)
  chrom <- data.frame(time, y)
  chrom$y[6] <- NA_real_
  outData1 <- splineFill(chrom, method = "fmm")
  expOutput <- data.frame(time, "intensity" = c(0.2050595, 0.8850070, 2.2068768, 3.7212677, 5.1652605,
                                                5.8089716, 5.5446804, 4.5671360, 3.3213154, 1.9485889,
                                                0.9520709, 0.3294218, 0.2009581, 0.1420923))
  expect_equal(outData1, expOutput, tolerance = 1e-03)

  chrom$y[1] <- NA_real_
  outData2 <- splineFill(chrom, method = "natural")
  expOutput <- data.frame(time, "intensity" = c(-0.3918884, 0.8850070, 2.2068768, 3.7212677, 5.1652605,
                                                5.8075496, 5.5446804, 4.5671360, 3.3213154, 1.9485889,
                                                0.9520709, 0.3294218, 0.2009581, 0.1420923))
  expect_equal(outData2, expOutput, tolerance = 1e-03)
})

test_that("test_approxFill", {
  time <- seq(from = 3003.4, to = 3048, by = 3.4)
  y <- c(0.2050595, 0.8850070, 2.2068768, 3.7212677, 5.1652605, 5.8288915, 5.5446804,
         4.5671360, 3.3213154, 1.9485889, 0.9520709, 0.3294218, 0.2009581, 0.1420923)
  chrom <- data.frame(time, y)
  chrom$y[6] <- NA_real_
  chrom$y[1] <- NA_real_
  outData <- approxFill(chrom)
  expOutput <- data.frame(time, "intensity" = c(0.8850070, 0.8850070, 2.2068768, 3.7212677, 5.1652605,
                                                5.3549704, 5.5446804, 4.5671360, 3.3213154, 1.9485889,
                                                0.9520709, 0.3294218, 0.2009581, 0.1420923))
  expect_equal(outData, expOutput, tolerance = 1e-03)
})

test_that("test_imputeChromatogram", {
  time <- seq(from = 3003.4, to = 3048, by = 3.4)
  y <- c(0.2050595, 0.8850070, 2.2068768, 3.7212677, 5.1652605, 5.8288915, 5.5446804,
         4.5671360, 3.3213154, 1.9485889, 0.9520709, 0.3294218, 0.2009581, 0.1420923)
  chrom <- data.frame(time, y)
  chrom$y[6] <- NA_real_
  chrom$y[1] <- NA_real_
  outData1 <- imputeChromatogram(chrom, method = "sgolay", polyOrd = 3, kernelLen = 9,
                                splineMethod = "fmm")
  expOutput <- data.frame(time, "intensity" = c(-0.0842669, 0.8850070, 2.2068768, 3.7212677, 5.1652605,
                                                5.3457775, 5.5446804, 4.5671360, 3.3213154, 1.9485889,
                                                0.9520709, 0.3294218, 0.2009581, 0.1420923))
  expect_equal(outData1, expOutput, tolerance = 1e-03)

  outData2 <- imputeChromatogram(chrom, method = "spline", splineMethod = "natural")
  expOutput <- data.frame(time, "intensity" = c(-0.3918884, 0.8850070, 2.2068768, 3.7212677, 5.1652605,
                                                5.8075496, 5.5446804, 4.5671360, 3.3213154, 1.9485889,
                                                0.9520709, 0.3294218, 0.2009581, 0.1420923))
  expect_equal(outData2, expOutput, tolerance = 1e-03)


  outData3 <- imputeChromatogram(chrom, method = "linear")
  expOutput <- data.frame(time, "intensity" = c(0.8850070, 0.8850070, 2.2068768, 3.7212677, 5.1652605,
                                                5.3549704, 5.5446804, 4.5671360, 3.3213154, 1.9485889,
                                                0.9520709, 0.3294218, 0.2009581, 0.1420923))
  expect_equal(outData3, expOutput, tolerance = 1e-03)
})
