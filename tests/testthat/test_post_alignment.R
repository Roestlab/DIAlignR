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
  df <- data.table::data.table(multipeptide_DIAlignR[["14383"]])
  setkeyv(df, "run")
  df[3, alignment_rank := 1L]
  XICs.ref <- XICs.eXp <- list()
  XICs.ref[["4618"]] <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]]
  XICs.eXp[["4618"]] <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]]
  alignObj <- testAlignObj()
  tAligned <- alignedTimes2(alignObj, XICs.ref[["4618"]], XICs.eXp[["4618"]])
  setAlignmentRank(df, refIdx = 3L, eXp = "run2", tAligned, XICs.eXp, params, adaptiveRT)
  expect_equal(c(NA_integer_, NA_integer_, 1L, NA_integer_, 1L, NA_integer_), df[,alignment_rank])

  # 2nd case
  df <- data.table::data.table(multipeptide_DIAlignR[["14383"]])
  df$alignment_rank[1] <- 1L; df$m_score[5] <- 0.03
  setAlignmentRank(df, refIdx = 1L, eXp = "run2", tAligned, XICs.eXp, params, adaptiveRT)
  expect_equal(c(1L, NA_integer_, NA_integer_, NA_integer_, 1L, NA_integer_), df[,alignment_rank])

  # case 3
  setAlignmentRank(df, refIdx = 1L, eXp = "run1", tAligned, XICs.eXp, params, adaptiveRT)
  expect_equal(c(1L, NA_integer_, 1L, NA_integer_, 1L, NA_integer_), df[,alignment_rank])

  # case 4
  df <- data.table::data.table(multipeptide_DIAlignR[["14383"]])
  df$alignment_rank[3] <- 1L; df$m_score[5] <- 0.06
  setAlignmentRank(df, refIdx = 3L, eXp = "run2", tAligned, XICs.eXp, params, adaptiveRT)
  expect_equal(df[,alignment_rank], c(NA_integer_, NA_integer_, 1L, NA_integer_, NA_integer_, 1L))
  expect_equal(df[6,], data.table("transition_group_id" = 4618L, feature_id = bit64::NA_integer64_,
                                       RT = 5241.30, intensity = 189.304,  leftWidth = 5224.2, rightWidth = 5265.2, peak_group_rank = NA_integer_,
                                       m_score = NA_real_, run = "run2", alignment_rank = 1, key = "run"), tolerance = 1e-06)

  # case 5
  df <- data.table::data.table(multipeptide_DIAlignR[["14383"]])
  df$m_score[5] <- 0.06
  setAlignmentRank(df, refIdx = 1L, eXp = "run2", tAligned, XICs.eXp, params, adaptiveRT)
  expect_equal(df[6,], data.table("transition_group_id" = 4618L, feature_id = bit64::NA_integer64_,
                                  RT = 5224.20, intensity = 99.77859,  leftWidth = 5203.7, rightWidth = 5251.5, peak_group_rank = NA_integer_,
                                  m_score = NA_real_, run = "run2", alignment_rank = 1, key = "run"), tolerance = 1e-06)


  # case 6
  df <- data.table::data.table(multipeptide_DIAlignR[["14383"]])
  df$alignment_rank[3] <- 1L; df$m_score[5] <- NA_real_
  setAlignmentRank(df, refIdx = 3L, eXp = "run2", tAligned, XICs.eXp, params, adaptiveRT)
  expect_equal(df[6,], data.table("transition_group_id" = 4618L, feature_id = bit64::NA_integer64_,
                                       RT = 5241.30, intensity = 189.304,  leftWidth = 5224.2, rightWidth = 5265.2, peak_group_rank = NA_integer_,
                                       m_score = NA_real_, run = "run2", alignment_rank = 1, key = "run"),
               tolerance = 1e-06)

  # case 7
  df <- data.table::data.table(multipeptide_DIAlignR[["14383"]])
  df$alignment_rank[3] <- 1L; df$m_score[3] <- NA_real_; df$m_score[5] <- 0.03
  setAlignmentRank(df, refIdx = 3L, eXp = "run2", tAligned, XICs.eXp, params, adaptiveRT)
  expect_equal(c(NA_integer_, NA_integer_, 1L, NA_integer_, 1L, NA_integer_), df[,alignment_rank])

  # case 8
  df <- data.table::data.table(multipeptide_DIAlignR[["14383"]])
  df$alignment_rank[3] <- 1L; df$m_score[1:6] <- NA_real_
  setAlignmentRank(df, refIdx = 3L, eXp = "run2", tAligned, XICs.eXp, params, adaptiveRT)
  expect_equal(df[6,], data.table("transition_group_id" = 4618L, feature_id = bit64::NA_integer64_,
                                       RT = 5241.30, intensity = 189.304,  leftWidth = 5224.2, rightWidth = 5265.2, peak_group_rank = NA_integer_,
                                       m_score = NA_real_, run = "run2", alignment_rank = 1, key = "run"),
              tolerance = 1e-06)

  # case 9
  df <- data.table::data.table(multipeptide_DIAlignR[["14383"]])
  df$m_score[1:6] <- NA_real_
  params$fillMissing <- FALSE
  setAlignmentRank(df, refIdx = 3L, eXp = "run2", tAligned, XICs.eXp, params, adaptiveRT)
  expect_equal(df$alignment_rank, rep(NA_integer_, 6))

  # case 10
  # bit64 does not return NA, instead returns 9218868437227407266 https://stackoverflow.com/a/27283100/6484844
  expect_error(setAlignmentRank(df, refIdx = integer(0), eXp = "run2", tAligned, XICs.eXp, params, adaptiveRT))
})

test_that("test_setOtherPrecursors", {
  data(multipeptide_DIAlignR, package="DIAlignR")
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  params <- paramsDIAlignR()

  df <- data.table::data.table(multipeptide_DIAlignR[["14383"]])
  XICs.eXp <- list()
  XICs.eXp[["4618"]] <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]]

  setOtherPrecursors(df, 5L, XICs.eXp, analytes = 4618L, params)
  expect_equal(df[,alignment_rank], rep(NA_integer_, 6))

  dataPath <- system.file("extdata", package = "DIAlignR")
  mz <- mzR::openMSfile(file.path(dataPath, "mzml","hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML"))
  df <- data.table::data.table(multipeptide_DIAlignR[["9861"]])
  chromIndices <- list(c(43, 44, 45, 46, 47, 48), c(49, 50, 51, 52, 53, 54))
  XICs.eXp <- lapply(chromIndices, function(i) extractXIC_group(mz, chromIndices = i))
  names(XICs.eXp) <- c("9719", "9720")
  df[10L, alignment_rank := 1L]
  setOtherPrecursors(df, 10L, XICs.eXp, analytes = c(9719L, 9720L), params)
  expect_equal(df[9L,], data.table("transition_group_id" = 9719L, feature_id = bit64::NA_integer64_,
       RT = 2607.05, intensity = 11.80541,  leftWidth = 2591.431, rightWidth = 2625.569, peak_group_rank = NA_integer_,
       m_score = NA_real_, run = "run2", alignment_rank = 1, key = "run"),
   tolerance = 1e-06)

  df <- data.table::data.table(multipeptide_DIAlignR[["9861"]])
  setOtherPrecursors(df, 10L, XICs.eXp, analytes = c(9719L, 9720L), params)
  expect_equal(df[,alignment_rank], c(rep(NA_integer_, 8), 1L, rep(NA_integer_, 3)))

  df[6L, alignment_rank := 1L]
  setOtherPrecursors(df, 6L, XICs.eXp, analytes = c(9719L, 9720L), params)
  expect_equal(df[,alignment_rank], c(rep(NA_integer_, 4), 1L, 1L, NA_integer_, NA_integer_, 1L, rep(NA_integer_, 3)))
})

test_that("test_reIntensity", {
  data(multipeptide_DIAlignR, package="DIAlignR")
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  params <- paramsDIAlignR()
  XICs.eXp <- list()
  df <- data.table::data.table(multipeptide_DIAlignR[["14383"]])
  XICs.eXp[["4618"]] <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]]

  reIntensity(df, "run2", XICs.eXp, params)
  expect_equal(df[5,intensity], 255.496)

  df[5, alignment_rank:= 1L]
  reIntensity(df, "run2", XICs.eXp, params)
  expect_equal(df[5L, intensity], 211.3709, tolerance = 1e-05)

  dataPath <- system.file("extdata", package = "DIAlignR")
  mz <- mzR::openMSfile(file.path(dataPath, "mzml","hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML"))
  df <- data.table::data.table(multipeptide_DIAlignR[["9861"]])
  chromIndices <- list(c(43, 44, 45, 46, 47, 48), c(49, 50, 51, 52, 53, 54))
  XICs.eXp <- lapply(chromIndices, function(i) extractXIC_group(mz, chromIndices = i))
  names(XICs.eXp) <- c("9719", "9720")

  df$alignment_rank[c(6L,10L)] <- 1L
  reIntensity(df, "run2", XICs.eXp, params)
  expect_equal(df[6L, intensity], 52.95950)
  expect_equal(df[10L, intensity], 24.50702, tolerance = 1e-05)

  mz <- mzR::openMSfile(file.path(dataPath, "mzml","hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML"))
  XICs.ref <- lapply(chromIndices, function(i) extractXIC_group(mz, chromIndices = i))
  names(XICs.ref) <- c("9719", "9720")

  reIntensity(df, "run1", XICs.ref, params)
  expect_equal(df[6L, intensity], 20.11727)
})
