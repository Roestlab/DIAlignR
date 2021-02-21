context("Utility functions")

test_that("test_getRefRun", {
  dataPath <- system.file("extdata", package = "DIAlignR")
  fileInfo <- getRunNames(dataPath, oswMerged = TRUE)
  precursors <- getPrecursors(fileInfo, oswMerged = TRUE, context = "experiment-wide", maxPeptideFdr = 1.00)
  peptideIDs <- unique(precursors$peptide_id)
  peptidesInfo <- getPeptideScores(fileInfo, peptideIDs, oswMerged = TRUE, runType = "DIA_Proteomics", context = "experiment-wide")
  peptidesInfo <- lapply(peptideIDs, function(pep) dplyr::filter(peptidesInfo, .data$peptide_id == pep))
  names(peptidesInfo) <- as.character(peptideIDs)
  outData <- getRefRun(peptidesInfo)

  expData <- data.table("peptide_id" = c(15L, 105L),
                        "run" = c("run1","run0"), key = "peptide_id")
  expect_identical(dim(outData), c(297L, 2L))
  expect_identical(outData[2:3,], expData)
})

test_that("test_getMultipeptide", {
  dataPath <- system.file("extdata", package = "DIAlignR")
  fileInfo <- getRunNames(dataPath, oswMerged = TRUE)
  precursors <- getPrecursors(fileInfo, oswMerged = TRUE, context = "experiment-wide", maxPeptideFdr = 0.05)
  features <- getFeatures(fileInfo, maxFdrQuery = 0.05)
  outData <- getMultipeptide(precursors, features)

  expData <- data.table("transition_group_id" = rep(c(9719L, 9720L),3),
                        "feature_id" = bit64::as.integer64(c(NA, "7742456255764097691", NA, NA, "6462000664077079508", "5135268764240690321", NA, NA, NA, "298844719207353347", NA, NA)),
                        "RT" = c(NA_real_, 2541.83, NA_real_, NA_real_, 2586.12, 2585.61, NA_real_, NA_real_, NA_real_, 2607.05, NA_real_, NA_real_),
                        "intensity" = c(NA_real_, 34.7208, NA_real_, NA_real_, 26.2182, 52.9595, NA_real_, NA_real_, NA_real_, 33.5961, NA_real_, NA_real_),
                        "leftWidth" = c(NA_real_, 2526.555, NA_real_, NA_real_, 2571.738, 2564.094, NA_real_, NA_real_, NA_real_, 2591.431, NA_real_, NA_real_),
                        "rightWidth" = c(NA_real_, 2560.693, NA_real_, NA_real_, 2609.288, 2605.060, NA_real_, NA_real_, NA_real_, 2625.569, NA_real_, NA_real_),
                        "peak_group_rank" = c(NA_integer_, 1L, NA_integer_, NA_integer_, 1L, 1L, NA_integer_, NA_integer_, NA_integer_, 1L, NA_integer_, NA_integer_),
                        "m_score" = c(NA_real_, 5.692077e-05, NA_real_, NA_real_, 1.041916e-03, 5.692077e-05, NA_real_, NA_real_, NA_real_, 2.005418e-04, NA_real_, NA_real_),
                        "run" = rep(c("run0", "run1", "run2"), each = 4),
                        "alignment_rank" = c(rep(NA_integer_,6)),
                        key = "run")
  expect_identical(length(outData), 229L)
  expect_equal(outData[[110]], expData, tolerance = 1e-04)
  expect_equal(outData[["9861"]], expData, tolerance = 1e-03)
})

test_that("test_writeTables", {
  # To check if the two files are same, we can use tools::md5sum(). However, this will
  # not tell us which value is different.
  dataPath <- system.file("extdata", package = "DIAlignR")
  fileInfo <- getRunNames(dataPath, oswMerged = TRUE)
  # data(multipeptide_DIAlignR, package="DIAlignR")
  precursors <- getPrecursors(fileInfo, oswMerged = TRUE, context = "experiment-wide")
  features <- getFeatures(fileInfo, maxFdrQuery = 0.05)
  allIDs <- unique(unlist(lapply(features, function(df) df[m_score <= 0.05, transition_group_id]),
                          recursive = FALSE, use.names = FALSE))
  precursors <- precursors[data.table::CJ(unique(peptide_id), allIDs) , nomatch=0L]
  multipeptide <- getMultipeptide(precursors, features)
  outData <- writeTables(fileInfo, multipeptide, precursors)
  expect_identical(nrow(outData), 500L)

  multipeptide[["7040"]][c(1,3), alignment_rank:= 1L]
  multipeptide[["3200"]][c(1,3, 5), alignment_rank:= 1L]
  outData <- writeTables(fileInfo, multipeptide, precursors)

  expect_identical(dim(outData), c(502L, 14L))
  expData <- data.table(peptide_id = c(3200L, 7040L), precursor = c(523L, 32L),
                        run = fileInfo$runName[1:2], RT = c(1529.54, NA_real_),
                        intensity = c(13.5686, NA_real_), leftWidth = c(1512.653,NA_real_), rightWidth= c(1546.791,NA_real_),
                        peak_group_rank = c(1L, NA_integer_), m_score = c(0.01440752, NA_real_), alignment_rank = 1L,
                        feature_id = c("2382660248384924660", NA_character_), sequence = c("DVDQYPR", "GNNSVYMNNFLNLILQNER"),
                        charge = c(2L, 3L), group_label = c("10483_DVDQYPR/2", "10030_GNNSVYMNNFLNLILQNER/3"),
                        key = NULL)
  expect_equal(outData[c(82,167),], expData, tolerance = 1e-05)
})

test_that("test_checkParams", {
})

test_that("test_paramsDIAlignR", {
})

test_that("test_alignmentStats", {
})


test_that("test_checkOverlap",{
  expect_true(checkOverlap(c(1.1, 3.1), c(2.1, 3.1)))
  expect_true(checkOverlap(c(1.1, 3.1), c(3.1, 4.1)))
  expect_true(checkOverlap(c(1.1, 3.1), c(2.1, 2.5)))
  expect_true(checkOverlap(c(1.1, 3.1), c(0.1, 9.1)))
  expect_false(checkOverlap(c(1.1, 3.1), c(3.2, 7.1)))

})
