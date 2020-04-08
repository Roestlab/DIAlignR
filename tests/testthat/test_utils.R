context("Utility functions")

test_that("test_getRefRun", {
  dataPath <- system.file("extdata", package = "DIAlignR")
  fileInfo <- getRunNames(dataPath, oswMerged = TRUE)
  precursors <- getPrecursors(fileInfo, oswMerged = TRUE)
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
  precursors <- getPrecursors(fileInfo, oswMerged = TRUE)
  features <- getFeatures(fileInfo, maxFdrQuery = 0.05)
  outData <- getMultipeptide(precursors, features)

  expData <- data.frame("transition_group_id" = 9723L,
                        "run" = c("run0", "run1", "run2"),
                        "RT" = c(NA_real_, 4057.14, NA_real_),
                        "intensity" = c(NA_real_, 36.1802, NA_real_),
                        "leftWidth" = c(NA_real_, 4049.51, NA_real_),
                        "rightWidth" = c(NA_real_, 4083.649, NA_real_),
                        "peak_group_rank" = c(NA_integer_, 1L, NA_integer_),
                        "m_score" = c(NA_real_, 0.03737512, NA_real_),
                        stringsAsFactors = TRUE)
  expect_identical(length(outData), 322L)
  expect_equal(outData[[150]], expData, tolerance = 1e-04)
  expect_equal(outData[["9723"]], expData, tolerance = 1e-04)
})


test_that("test_selectChromIndices", {
  data(oswFiles_DIAlignR, package="DIAlignR")
  oswFiles <- oswFiles_DIAlignR
  outData <- selectChromIndices(oswFiles, runname = "run0", analyte = "14299_QFNNTDIVLLEDFQK/3")
  expData <- c(37L,38L,39L,40L,41L,42L)
  expect_identical(outData, expData)
  outData <- selectChromIndices(oswFiles, runname = "run2", analyte = "AQPVPAHVY_2")
  expect_identical(outData, NULL)
})


test_that("test_writeTables", {

  # To check if the two files are same, we can use tools::md5sum(). However, this will
  # not tell us which value is different.
  refAnalytes <- c("14299_QFNNTDIVLLEDFQK/3", "13597_VVAGGELFKESVVVNDK/3")
  runs <- c("run0" = "hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt",
           "run1" = "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
           "run2" = "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
  filename <- "temp.csv"
  intensityTbl <-  list(c(157.864, 310.010, 255.496), c(NA, 30.1907, 34.6354))
  rtTbl <-  list(c(5222.12, 5238.35, 5240.79), c(NA, 3658.39, 3701.29))
  lwTbl <-  list(c(5203.6821289, 5220.7578125, 5217.3608398), c(NA, 3657.1279297, 3674.2170410))
  rwTbl <-  list(c(5248.0629883, 5261.7231445, 5275.3950195), c(NA, 3701.5061035, 3715.1831055))
  writeTables(refAnalytes, runs, filename, intensityTbl, rtTbl, lwTbl, rwTbl)

  outData <- read.table("temp.csv", stringsAsFactors = FALSE, sep = ",", header = TRUE)
  expData <- read.table("test2.csv", stringsAsFactors = FALSE, sep = ",", header = TRUE)
  expect_identical(dim(outData), dim(expData))
  expect_identical(colnames(outData), colnames(expData))
  expect_identical(outData[["peptide"]], expData[["peptide"]])
  expect_identical(outData[["run"]], expData[["run"]])
  for(i in 3:6){
    expect_equal(outData[[i]], expData[[i]], tolerance = 1e-04)
  }
  file.remove("temp.csv")
})
