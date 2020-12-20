context("Post alignment analysis")

# mzML and osw files are not required.

test_that("test_pickNearestFeature", {
  data(oswFiles_DIAlignR, package="DIAlignR")
  oswFiles <- oswFiles_DIAlignR
  outData <- pickNearestFeature(eXpRT = 5237.8, analyte = 4618L, oswFiles,
                                runname = "run2", adaptiveRT = 77.82315,
                                featureFDR = 0.05)
  expData <- list("leftWidth" = c(5217.361),
                           "rightWidth" = c(5275.395),
                           "RT" = c(5240.79),
                           "intensity" = c(255.496),
                           "peak_group_rank" = c(1L),
                           "m_score" = c(5.692077e-05))
  expect_equal(outData, expData, tolerance = 1e-05)
})

test_that("test_mapIdxToTime", {
  timeVec <- c(1.3,5.6,7.8)
  idx <- c(NA, NA, 1L, 2L, NA, NA, 3L, NA)
  outData <- mapIdxToTime(timeVec, idx)
  expData <- c(NA, NA, 1.3, 5.6, 6.333333, 7.066667, 7.8, NA)
  expect_equal(outData, expData, tolerance = 1e-04)
})

test_that("test_mappedRTfromAlignObj", {
  AlignObj <- testAlignObj()
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR
  tVec.ref <-  XICs[["hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]][[1]][, "time"]
  tVec.eXp <-  XICs[["hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]][[1]][, "time"]
  expect_equal(mappedRTfromAlignObj(refRT= 5238.35, tVec.ref, tVec.eXp, AlignObj), 5241.3)
})

test_that("test_setAlignmentRank", {
  data(multipeptide_DIAlignR, package="DIAlignR")
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  params <- paramsDIAlignR()
  adaptiveRT <- 38.66
  df <- setDT(multipeptide_DIAlignR[["14383"]])
  df[2, alignment_rank := 1L]
  XICs.ref <- XICs.eXp <- list()
  XICs.ref[["4618"]] <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]]
  XICs.eXp[["4618"]] <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]]
  alignObj <- testAlignObj()
  tAligned <- alignedTimes2(alignObj, XICs.ref[["4618"]], XICs.eXp[["4618"]])
  outData <- setAlignmentRank(df, ref = "run1", eXp = "run2", tAligned, XICs.eXp, params, adaptiveRT)
  df[3, alignment_rank := 1L]
  expect_equal(df[3,], outData)

  # 2nd case
  data(multipeptide_DIAlignR, package="DIAlignR")
  df <- setDT(multipeptide_DIAlignR[["14383"]])
  df$alignment_rank[1] <- 1L; df$m_score[3] <- 0.03
  outData <- setAlignmentRank(df, ref = "run1", eXp = "run2", tAligned, XICs.eXp, params, adaptiveRT)
  expect_equal(df[3,], outData)

  # case 3
  outData <- setAlignmentRank(df, ref = "run0", eXp = "run2", tAligned, XICs.eXp, params, adaptiveRT)
  df[3, alignment_rank := 1L]
  expect_equal(df[3,], outData)

  # case 4
  data(multipeptide_DIAlignR, package="DIAlignR")
  df <- setDT(multipeptide_DIAlignR[["14383"]])
  df$alignment_rank[2] <- 1L; df$m_score[3] <- 0.06
  outData <- setAlignmentRank(df, ref = "run1", eXp = "run2", tAligned, XICs.eXp, params, adaptiveRT)
  expect_equal(df[3,], outData[1,])
  expect_equal(outData[2,], data.table("transition_group_id" = 4618L, feature_id = bit64::NA_integer64_,
                                       RT = 5241.30, intensity = 189.304,  leftWidth = 5224.2, rightWidth = 5265.2, peak_group_rank = NA_integer_,
                                       m_score = NA_real_, run = "run2", alignment_rank = 1),   tolerance = 1e-06)

  # case 5
  data(multipeptide_DIAlignR, package="DIAlignR")
  df <- setDT(multipeptide_DIAlignR[["14383"]])
  df$m_score[3] <- 0.06
  outData <- setAlignmentRank(df, ref = "run0", eXp = "run2", tAligned, XICs.eXp, params, adaptiveRT)
  expect_equal(df[3,], outData)

  # case 6
  data(multipeptide_DIAlignR, package="DIAlignR")
  df <- setDT(multipeptide_DIAlignR[["14383"]])
  df$alignment_rank[2] <- 1L; df$m_score[3] <- NA_real_
  outData <- setAlignmentRank(df, ref = "run1", eXp = "run2", tAligned, XICs.eXp, params, adaptiveRT)
  expect_equal(df[3,], outData[1,])
  expect_equal(outData[2,], data.table("transition_group_id" = 4618L, feature_id = bit64::NA_integer64_,
                                       RT = 5241.30, intensity = 189.304,  leftWidth = 5224.2, rightWidth = 5265.2, peak_group_rank = NA_integer_,
                                       m_score = NA_real_, run = "run2", alignment_rank = 1),
               tolerance = 1e-06)

  # case 7
  data(multipeptide_DIAlignR, package="DIAlignR")
  df <- setDT(multipeptide_DIAlignR[["14383"]])
  df$alignment_rank[2] <- 1L; df$m_score[2] <- NA_real_; df$m_score[3] <- 0.03
  outData <- setAlignmentRank(df, ref = "run1", eXp = "run2", tAligned, XICs.eXp, params, adaptiveRT)
  df$alignment_rank[3] <- 1L
  expect_equal(df[3,], outData)

  # case 8
  data(multipeptide_DIAlignR, package="DIAlignR") #
  df <- setDT(multipeptide_DIAlignR[["14383"]])
  df$alignment_rank[2] <- 1L; df$m_score[1:3] <- NA_real_
  outData <- setAlignmentRank(df, ref = "run1", eXp = "run2", tAligned, XICs.eXp, params, adaptiveRT)
  expect_equal(df[3,], outData[1,])
  expect_equal(outData[2,], data.table("transition_group_id" = 4618L, feature_id = bit64::NA_integer64_,
                                       RT = 5241.30, intensity = 189.304,  leftWidth = 5224.2, rightWidth = 5265.2, peak_group_rank = NA_integer_,
                                       m_score = NA_real_, run = "run2", alignment_rank = 1),
               tolerance = 1e-06)
  # case 9
  data(multipeptide_DIAlignR, package="DIAlignR") # data.table get modified in-place
  df <- setDT(multipeptide_DIAlignR[["14383"]])
  df$m_score[1:3] <- NA_real_
  outData <- setAlignmentRank(df, ref = "run1", eXp = "run2", tAligned, XICs.eXp, params, adaptiveRT)
  expect_equal(outData, df[3,])
  # case 10
  # bit64 does not return NA, instead returns 9218868437227407266 https://stackoverflow.com/a/27283100/6484844
  outData <- setAlignmentRank(df[4,], ref = "run1", eXp = "run2", tAligned, XICs.eXp, params, adaptiveRT)
  expData <- data.frame("transition_group_id" = NA_integer_, feature_id = bit64::as.integer64("9218868437227407266"),
                        RT = NA_real_, intensity = NA_real_,  leftWidth = NA_real_, rightWidth = NA_real_, peak_group_rank = NA_integer_,
                        m_score = NA_real_, run = NA_character_, alignment_rank = NA_integer_)
  # expect_equal(outData, expData)
})

test_that("test_setOtherPrecursors", {
  data(multipeptide_DIAlignR, package="DIAlignR")
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  params <- paramsDIAlignR()

  df <- setDT(multipeptide_DIAlignR[["14383"]])
  XICs.eXp <- list()
  XICs.eXp[["4618"]] <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]]

  outData <- setOtherPrecursors(df[3,], XICs.eXp, analytes = 4618L, params)
  data.table::setindex(outData, NULL)
  expect_equal(outData, df[3,])

  dataPath <- system.file("extdata", package = "DIAlignR")
  mz <- mzR::openMSfile(file.path(dataPath, "mzml","hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML"))
  df <- setDT(multipeptide_DIAlignR[["9861"]])
  chromIndices <- list(c(43, 44, 45, 46, 47, 48), c(49, 50, 51, 52, 53, 54))
  XICs.eXp <- lapply(chromIndices, function(i) extractXIC_group(mz, chromIndices = i))
  names(XICs.eXp) <- c("9719", "9720")
  df[4,"alignment_rank"] <- 1L
  outData <- setOtherPrecursors(df[4,], XICs.eXp, analytes = c(9719L, 9720L), params)
  data.table::setindex(outData, NULL)
  expect_equal(df[4,], outData[1,])
  expect_equal(outData[2,], data.table("transition_group_id" = 9719L, feature_id = bit64::NA_integer64_,
                                       RT = 2607.05, intensity = 11.80541,  leftWidth = 2591.431, rightWidth = 2625.569, peak_group_rank = NA_integer_,
                                       m_score = NA_real_, run = "run2", alignment_rank = 1),
               tolerance = 1e-06)

  data(multipeptide_DIAlignR, package="DIAlignR")
  df <- setDT(multipeptide_DIAlignR[["9861"]])
  outData <- setOtherPrecursors(df[4,], XICs.eXp, analytes = c(9719L, 9720L), params)
  data.table::setindex(outData, NULL)
  expect_equal(df[4,], outData)

  df[2:3,"alignment_rank"] <- 1L
  outData <- setOtherPrecursors(df[2:3,], XICs.eXp, analytes = c(9719L, 9720L), params)
  data.table::setindex(outData, NULL)
  expect_equal(df[2:3,], outData)
})

test_that("test_reIntensity", {
  data(multipeptide_DIAlignR, package="DIAlignR")
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  params <- paramsDIAlignR()
  XICs.eXp <- list()
  df <- setDT(multipeptide_DIAlignR[["14383"]])
  XICs.eXp[["4618"]] <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]]

  outData <- reIntensity(df[3,], XICs.eXp, params)
  data.table::setindex(outData, NULL)
  expect_equal(outData, df[3,])

  df[3, alignment_rank:= 1L]
  outData <- reIntensity(df[3,], XICs.eXp, params)
  expect_equal(outData[1, intensity], 211.3709, tolerance = 1e-05)

  dataPath <- system.file("extdata", package = "DIAlignR")
  mz <- mzR::openMSfile(file.path(dataPath, "mzml","hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML"))
  df <- setDT(multipeptide_DIAlignR[["9861"]])
  chromIndices <- list(c(43, 44, 45, 46, 47, 48), c(49, 50, 51, 52, 53, 54))
  XICs.eXp <- lapply(chromIndices, function(i) extractXIC_group(mz, chromIndices = i))
  names(XICs.eXp) <- c("9719", "9720")

  df.eXp <- df[4,]
  outData <- reIntensity(df.eXp, XICs.eXp, params)
  expect_equal(outData, df.eXp)

  df.eXp$alignment_rank <- 1L
  outData <- reIntensity(df.eXp, XICs.eXp, params)
  expect_equal(outData[1, intensity], 24.50702, tolerance = 1e-05)

  mz <- mzR::openMSfile(file.path(dataPath, "mzml","hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML"))
  XICs.ref <- lapply(chromIndices, function(i) extractXIC_group(mz, chromIndices = i))
  names(XICs.ref) <- c("9719", "9720")

  df.ref <- df[2:3,]
  df.ref$alignment_rank[2] <- 1L
  outData <- reIntensity(df.ref, XICs.ref, params)
  df.ref$intensity[2] <- 20.11727
  expect_equal(outData, df.ref)
})














