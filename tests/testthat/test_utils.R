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
