context("Utility functions")

test_that("test_getRefRun", {
  dataPath <- system.file("extdata", package = "DIAlignR")
  fileInfo <- getRunNames(dataPath, oswMerged = TRUE)
  precursors <- getPrecursors(fileInfo, oswMerged = TRUE, context = "experiment-wide", maxPeptideFdr = 1.00)
  features <- getFeatures(fileInfo, maxFdrQuery = 0.05)
  allIDs <- unique(unlist(lapply(features, `[[`, "transition_group_id"),
                          recursive = FALSE, use.names = FALSE))
  precursors <- precursors[precursors[["transition_group_id"]] %in% allIDs, ]
  multipeptide <- getMultipeptide(precursors, features)
  outData <- getRefRun(multipeptide)

  expData <- data.frame("transition_group_id" = c(32L, 470L),
                        "run" = c("run0","run0"),
                        stringsAsFactors = FALSE)
  expect_identical(dim(outData), c(199L, 2L))
  expect_identical(outData[1:2,], expData)
})

test_that("test_getMultipeptide", {
  dataPath <- system.file("extdata", package = "DIAlignR")
  fileInfo <- getRunNames(dataPath, oswMerged = TRUE)
  precursors <- getPrecursors(fileInfo, oswMerged = TRUE, context = "experiment-wide", maxPeptideFdr = 1.00)
  features <- getFeatures(fileInfo, maxFdrQuery = 0.05)
  outData <- getMultipeptide(precursors, features)

  expData <- data.frame("transition_group_id" = 9723L,
                        "feature_id" = c(bit64::NA_integer64_, bit64::as.integer64(5930589188108275441), bit64::NA_integer64_),
                        "run" = c("run0", "run1", "run2"),
                        "RT" = c(NA_real_, 4057.14, NA_real_),
                        "intensity" = c(NA_real_, 36.1802, NA_real_),
                        "leftWidth" = c(NA_real_, 4049.51, NA_real_),
                        "rightWidth" = c(NA_real_, 4083.649, NA_real_),
                        "peak_group_rank" = c(NA_integer_, 1L, NA_integer_),
                        "m_score" = c(NA_real_, 0.03737512, NA_real_),
                        "alignment_rank" = c(rep(NA_integer_,3)),
                        stringsAsFactors = FALSE)
  expect_identical(length(outData), 312L)
  expect_equal(outData[[145]], expData, tolerance = 1e-04)
  expect_equal(outData[["9723"]], expData, tolerance = 1e-04)
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

  expData <- read.table("test.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE)
  expect_identical(dim(outData), dim(expData))
  expect_identical(colnames(outData), colnames(expData))
  expect_identical(outData[["peptide"]], expData[["peptide"]])
  expect_identical(outData[["run"]], expData[["run"]])
  # Not checking intensity and other scores as alignment is not yet done.
  for(i in 10:13){
    expect_equal(outData[[i]], expData[[i]], tolerance = 1e-04)
  }
})

test_that("test_checkParams", {
})

test_that("test_paramsDIAlignR", {
})

test_that("test_alignmentStats", {
})

