context("Utility functions")

test_that("test_getRefRun", {
  dataPath <- system.file("extdata", package = "DIAlignR")
  fileInfo <- getRunNames(dataPath, oswMerged = TRUE)
  precursors <- getPrecursors(fileInfo, oswMerged = TRUE, context = "experiment-wide", maxPeptideFdr = 1.00)
  peptideIDs <- unique(precursors$peptide_id)
  peptidesInfo <- getPeptideScores(fileInfo, peptideIDs, oswMerged = TRUE, runType = "DIA_proteomics", context = "experiment-wide")
  peptidesInfo <- lapply(peptideIDs, function(pep) dplyr::filter(peptidesInfo, .data$peptide_id == pep))
  names(peptidesInfo) <- as.character(peptideIDs)
  outData <- getRefRun(peptidesInfo)

  expData <- data.frame("peptide_id" = c(7040L, 4170L),
                        "run" = c("run0","run2"),
                        stringsAsFactors = FALSE)
  expect_identical(dim(outData), c(297L, 2L))
  expect_identical(outData[1:2,], expData)
})

test_that("test_getMultipeptide", {
  dataPath <- system.file("extdata", package = "DIAlignR")
  fileInfo <- getRunNames(dataPath, oswMerged = TRUE)
  precursors <- getPrecursors(fileInfo, oswMerged = TRUE, context = "experiment-wide", maxPeptideFdr = 0.05)
  features <- getFeatures(fileInfo, maxFdrQuery = 0.05)
  outData <- getMultipeptide(precursors, features)

  expData <- data.frame("transition_group_id" = c(9720L, 9719L, 9720L, 9720L),
                        "feature_id" = bit64::as.integer64(c(7742456255764097691, 6462000664077079508, 5135268764240690321, 298844719207353347)),
                        "RT" = c(2541.83, 2586.12, 2585.61, 2607.05),
                        "intensity" = c(34.7208, 26.2182, 52.9595, 33.5961),
                        "leftWidth" = c(2526.555, 2571.738, 2564.094, 2591.431),
                        "rightWidth" = c(2560.693, 2609.288, 2605.060, 2625.569),
                        "peak_group_rank" = c(1L),
                        "m_score" = c(5.692077e-05, 1.041916e-03, 5.692077e-05, 2.005418e-04),
                        "run" = c("run0", "run1", "run1", "run2"),
                        "alignment_rank" = c(rep(NA_integer_,4)),
                        stringsAsFactors = FALSE)
  expect_identical(length(outData), 229L)
  expect_equal(outData[[108]], expData, tolerance = 1e-04)
  expect_equal(outData[["9861"]], expData, tolerance = 1e-04)
})

test_that("test_writeTables", {
  # To check if the two files are same, we can use tools::md5sum(). However, this will
  # not tell us which value is different.

  dataPath <- system.file("extdata", package = "DIAlignR")
  fileInfo <- getRunNames(dataPath, oswMerged = TRUE)
  # data(multipeptide_DIAlignR, package="DIAlignR")
  precursors <- getPrecursors(fileInfo, oswMerged = TRUE, context = "experiment-wide")
  features <- getFeatures(fileInfo, maxFdrQuery = 0.05)
  allIDs <- unique(unlist(lapply(features, function(df) df[df[["m_score"]] <= 0.05,
                          "transition_group_id"]),
                          recursive = FALSE, use.names = FALSE))
  precursors <- precursors[precursors[["transition_group_id"]] %in% allIDs, ]
  multipeptide <- getMultipeptide(precursors, features)
  outData <- writeTables(fileInfo, multipeptide, precursors)
  expect_identical(nrow(outData), 500L)

  multipeptide[["7040"]]$alignment_rank[1:2] <- 1L
  multipeptide[["3200"]]$alignment_rank <- 1L
  outData <- writeTables(fileInfo, multipeptide, precursors)

  expect_identical(dim(outData), c(501L, 14L))
  expData <- data.frame(peptide_id = c(3200L, 7040L), precursor = c(523L, 32L),
                        run = fileInfo$runName[1:2], RT = c(1529.54, NA_real_),
                        intensity = c(13.5686, NA_real_), leftWidth = c(1512.653,NA_real_), rightWidth= c(1546.791,NA_real_),
                        peak_group_rank = c(1L, NA_integer_), m_score = c(0.01440752, NA_real_), alignment_rank = 1L,
                        feature_id = c("2382660248384924660", NA_character_), sequence = c("DVDQYPR", "GNNSVYMNNFLNLILQNER"),
                        charge = c(2L, 3L), group_label = c("10483_DVDQYPR/2", "10030_GNNSVYMNNFLNLILQNER/3"), row.names = c(82L, 167L))
  expect_equal(outData[c(82,167),], expData, tolerance = 1e-05)
})

test_that("test_checkParams", {
})

test_that("test_paramsDIAlignR", {
})

test_that("test_alignmentStats", {
})

