context("Utility functions")

test_that("test_getRefRun", {
  data(oswFiles_DIAlignR, package="DIAlignR")
  oswFiles <- oswFiles_DIAlignR
  expect_identical(getRefRun(oswFiles, analyte = "14299_QFNNTDIVLLEDFQK/3"), 1L)
  expect_identical(getRefRun(oswFiles, analyte = "AQPPVSTEY_2"), NULL)
  expect_identical(getRefRun(oswFiles, analyte = "19051_KLIVTSEGC[160]FK/2"), 2L)
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
