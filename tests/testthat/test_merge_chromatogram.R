context("Merge two chromatograms")

test_that("test_childXICs",{
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  data(alignObj_DIAlignR, package="DIAlignR")
  XICs.ref <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]]
  XICs.eXp <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]]
  alignedIndices <- cbind(alignObj_DIAlignR@indexA_aligned, alignObj_DIAlignR@indexB_aligned)
  colnames(alignedIndices) <- c("indexAligned.ref", "indexAligned.eXp")
  alignedIndices[, 1:2][alignedIndices[, 1:2] == 0] <- NA_integer_
  outData <- childXICs(XICs.ref, XICs.eXp, alignedIndices)

  expect_identical(length(outData), 2L)
  expect_identical(dim(outData[[1]][[6]]), c(177L, 2L))

  data(masterXICs_DIAlignR, package="DIAlignR")
  expData <- masterXICs_DIAlignR
  expect_equal(outData, expData)
})

test_that("test_childXIC",{
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  data(alignObj_DIAlignR, package="DIAlignR")
  XIC.ref <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]][[1]]
  XIC.eXp <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]][[1]]
  alignedIndices <- cbind(alignObj_DIAlignR@indexA_aligned, alignObj_DIAlignR@indexB_aligned)
  colnames(alignedIndices) <- c("indexAligned.ref", "indexAligned.eXp")
  alignedIndices[, 1:2][alignedIndices[, 1:2] == 0] <- NA_integer_

  outData <- childXIC(XIC.ref, XIC.eXp, alignedIndices, keepFlanks = TRUE)
  expect_identical(length(outData), 2L)
  # 177 = 176 + 1. 176 indices from Reference chromatogram.
  # +1 from the flanking signal at the tail of experiment chromatogram.
  expect_identical(dim(outData[[1]]), c(177L, 2L))
  expect_identical(dim(outData[[2]]), c(205L, 5L))
  expData <- as.matrix(data.frame(indexAligned.ref = c(NA_integer_, 176L, NA_integer_),
             indexAligned.eXp = c(174L, 175L, 176L),
             tAligned.ref = c(NA_real_, 5575.8, NA_real_),
             tAligned.eXp = c(5579.2, 5582.6, 5586.1),
             alignedChildTime = c(NA_real_, 5579.2, 5582.70)))
  expect_equal(outData[[2]][203:205,], expData)

  outData2 <- childXIC(XIC.ref, XIC.eXp, alignedIndices, keepFlanks = FALSE)
  # 158 = 176 -19 +1. 176 indices from Reference chromatogram. 1:18 indices are flanking in the reference chromatogram.
  expect_identical(dim(outData2[[1]]), c(158L, 2L))
  expData <- as.matrix(data.frame(indexAligned.ref = c(NA_integer_, 176L, NA_integer_),
                                  indexAligned.eXp = c(174L, 175L, 176L),
                                  tAligned.ref = c(NA, 5575.8, NA),
                                  tAligned.eXp = c(5579.2, 5582.6, 5586.1),
                                  alignedChildTime = c(NA, 5579.2, NA)))
  expect_equal(outData2[[2]][203:205,], expData)
})

test_that("test_addFlankToLeft",{
  time <- seq(from = 3003.4, to = 3048, by = 3.4)
  y <- c(0.2050595, 0.8850070, 2.2068768, 3.7212677, 5.1652605, 5.8288915, 5.5446804,
         4.5671360, 3.3213154, 1.9485889, 0.9520709, 0.3294218, 0.2009581, 0.1420923)
  chrom <- data.frame(time, y)
  chrom2 <- data.frame(time = c(3013.4, 3016, 3020), intensity = c(1.2, 3.4, 5.6))
  flankSeq <- as.logical(c(1,1,0,0,0,0,0,0,0,0,0,0,1,1))
  outData <- addFlankToLeft(flankSeq, chrom, chrom2)
  expData <-  data.frame(time = c(3006.6, 3010, 3013.4, 3016, 3020),
                         intensity = c(0.2050595, 0.885007, 1.2, 3.4, 5.6))
  expect_equal(outData, expData, tolerance = 1e-03)

  flankSeq <- as.logical(c(0,0,0,0,0,0,0,0,0,0,0,0,1,1))
  expect_error(addFlankToLeft(flankSeq, chrom, chrom2))
})

test_that("test_addFlankToRight",{
  time <- seq(from = 3003.4, to = 3048, by = 3.4)
  y <- c(0.2050595, 0.8850070, 2.2068768, 3.7212677, 5.1652605, 5.8288915, 5.5446804,
         4.5671360, 3.3213154, 1.9485889, 0.9520709, 0.3294218, 0.2009581, 0.1420923)
  chrom <- data.frame(time, y)
  chrom2 <- data.frame(time = c(3013.4, 3016, 3020), intensity = c(1.2, 3.4, 5.6))
  flankSeq <- as.logical(c(1,1,0,0,0,0,0,0,0,0,0,0,1,1))
  outData <- addFlankToRight(flankSeq, chrom, chrom2)
  expData <-  data.frame(time = c(3013.4, 3016, 3020, 3023.4, 3026.8),
                         intensity = c(1.2, 3.4, 5.6, 0.2009581, 0.1420923))
  expect_equal(outData, expData, tolerance = 1e-03)

  flankSeq <- as.logical(c(1,1,0,0,0,0,0,0,0,0,0,0,0,0))
  expect_error(addFlankToRight(flankSeq, chrom, chrom2))
})

test_that("test_mergeXIC",{
  t1 <- c(seq(from = 3003.4, by = 3.4, length.out = 10), 3035,
          seq(from = 3037.4, by = 3.4, length.out = 4))
  y1 <- c(0.2050595, 0.8850070, 2.2068768, 3.7212677, 5.1652605, 5.8288915, 5.5446804,
         4.5671360, 3.3213154, 1.9485889, 1.5, 0.9520709, 0.3294218, 0.2009581, 0.1420923)
  chrom1 <- data.frame(time = t1, y = y1)

  t2 <-  c(seq(from = 3053.4, by = 3.4, length.out = 10), 3084.5,
           seq(from = 3087.4, by = 3.4, length.out = 4))
  set.seed(1)
  y2 <- y1 + rnorm(length(y1), mean = 1)
  chrom2 <- data.frame(time = t2, y = y2)

  intensity <- rowMeans(cbind(y1, y2))
  outData1 <- mergeXIC(chrom1, chrom2, w.ref = 0.5, mergeStrategy = "ref")
  expect_equal(outData1, data.frame(time = t1, intensity))
  outData2 <- mergeXIC(chrom1, chrom2, w.ref = 0.5, mergeStrategy = "avg")
  expect_equal(outData2, data.frame(time = rowMeans(cbind(t1, t2)), intensity))
  outData3 <- mergeXIC(chrom1, chrom2, w.ref = 0.5, mergeStrategy = "refStart")
  expect_equal(outData3, data.frame(time = seq(from = 3028.4, by = 3.4, length.out = 15),
                                    intensity), tolerance = 1e-03)
  outData4 <- mergeXIC(chrom1, chrom2, w.ref = 0.5, mergeStrategy = "refEnd")
  time <- rev(seq(from = 3072.6, by = -3.4, length.out = 15))
  expect_equal(outData4, data.frame(time, intensity), tolerance = 1e-03)

  outData5 <- mergeXIC(chrom1, chrom2, w.ref = 0, mergeStrategy = "ref")
  expect_equal(outData5, data.frame(time = t1, intensity = y2))

  expect_error(mergeXIC(chrom1, chrom2, w.ref = 0.5, mergeStrategy = "none"),
               "mergeStrategy is not correct. Must be selected from 'ref', 'avg', 'refStart' and 'refEnd'.")
})

test_that("test_alignedXIC",{
  time <- seq(from = 3003.4, to = 3048, by = 3.4)
  y <- c(0.2050595, 0.8850070, 2.2068768, 3.7212677, 5.1652605, 5.8288915, 5.5446804,
    4.5671360, 3.3213154, 1.9485889, 0.9520709, 0.3294218, 0.2009581, 0.1420923)
  chrom <- data.frame(time, y)
  indices <- as.integer(c(NA, 1,2,3,4, NA, 5,6,7, NA, NA, 8,9,10,11,12,13, NA, 14, NA, NA))
  outData <- alignedXIC(chrom, indices, method = "spline", polyOrd = 4, kernelLen = 9, splineMethod = "fmm")

  time <- c(NA, 3003.4, 3006.8, 3010.2, 3013.6, 3015.3, 3017, 3020.4, 3023.8,
            3024.93, 3026.07, 3027.2, 3030.6, 3034, 3037.4, 3040.8, 3044.2, 3045.9,
            3047.6, NA, NA)
  y <- c(NA, 0.2051, 0.885, 2.2069, 3.7213, 4.4958, 5.1653, 5.8289, 5.5447, 5.2724, 4.9389,
         4.5671, 3.3213, 1.9486, 0.9521, 0.3294, 0.201, 0.1938, 0.1421, NA, NA)
  expect_equal(outData, data.frame(time, y), tolerance = 1e-03)

})
