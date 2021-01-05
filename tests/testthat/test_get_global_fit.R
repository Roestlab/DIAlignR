context("Global Fit")

test_that("test_getLOESSfit", {
    data(oswFiles_DIAlignR, package="DIAlignR")
    oswFiles <- oswFiles_DIAlignR
    Loess.fit <- getLOESSfit(oswFiles, ref = "run1", eXp = "run2", maxFdrGlobal = 0.05, spanvalue = 0.1)
    # Testing for Residual standard error
    expect_equal(getRSE(Loess.fit, "loess"), 22.00103, tolerance = 1e-05)
    expect_identical(names(Loess.fit), c("x", "y", "RT.ref", "RT.eXp"))
    # Add predict function as well.
    lfun <- stats::approxfun(Loess.fit)
    expect_equal(lfun(4978.4), 4966.134, tolerance = 1e-05)
    expect_equal(lfun(5575.8), 5566.859, tolerance = 1e-05)
  })

test_that("test_dialignrLoess", {
  df <- data.frame("transition_group_id" = 1:10, "RT.eXp" = 2:11, "RT.ref" = 10:19)
  # Testing for loess
  outData <- dialignrLoess(df, 0.1)
  expect_equal(outData$pars$span, 0.8, tolerance = 1e-05)
  expect_equal(predict(outData, newdata = data.frame("RT.ref"= 13.5))[[1]], 5.5, tolerance = 1e-05)
  # Testing for linear
  outData <- dialignrLoess(df[1:4,], 0.1)
  expect_equal(outData$coefficients[1][[1]], -8, tolerance = 1e-02)
})


test_that("test_getLinearfit", {
  data(oswFiles_DIAlignR, package="DIAlignR")
  oswFiles <- oswFiles_DIAlignR
  Linear.fit <- getLinearfit(oswFiles, ref = "run1", eXp = "run2", maxFdrGlobal = 0.05)
  # Testing for Residual standard error
  expect_equal(getRSE(Linear.fit, "linear"), 30.12705, tolerance = 1e-05)
  outData <- sum(coef(Linear.fit)*c(1,4978.4))
  expect_equal(outData, 4990.682, tolerance = 1e-05)
  outData <- sum(coef(Linear.fit)*c(1,5575.8))
  expect_equal(outData, 5577.561, tolerance = 1e-05)
})

test_that("test_getGlobalAlignment", {
  data(oswFiles_DIAlignR, package="DIAlignR")
  oswFiles <- oswFiles_DIAlignR
  globalFit <- getGlobalAlignment(oswFiles, ref = "run1", eXp = "run2", maxFdrGlobal = 0.05, spanvalue = 0.1, fitType = "loess")
  # Testing for Residual standard error
  expect_equal(getRSE(globalFit, "loess"), 22.00103, tolerance = 1e-05)
  # Add predict function as well.
  lfun <- stats::approxfun(globalFit)
  expect_equal(lfun(4978.4), 4966.134, tolerance = 1e-05)
  expect_equal(lfun(5575.8), 5566.859, tolerance = 1e-05)

  globalFit <- getGlobalAlignment(oswFiles, ref = "run1", eXp = "run2", maxFdrGlobal = 0.05, fitType = "linear")
  # Testing for Residual standard error
  expect_equal(getRSE(globalFit, "linear"), 30.12705, tolerance = 1e-05)
  outData <- sum(coef(globalFit)*c(1,4978.4))
  expect_equal(outData, 4990.682, tolerance = 1e-05)
  outData <- sum(coef(globalFit)*c(1,5575.8))
  expect_equal(outData, 5577.561, tolerance = 1e-05)
})


test_that("test_getRSE", {
  data(oswFiles_DIAlignR, package="DIAlignR")
  oswFiles <- oswFiles_DIAlignR
  globalFit <- getGlobalAlignment(oswFiles, ref = "run1", eXp = "run2", maxFdrGlobal = 0.05, spanvalue = 0.1, fitType = "loess")
  expect_equal(getRSE(globalFit, "loess"), 22.00103, tolerance = 1e-05)

  globalFit <- getGlobalAlignment(oswFiles, ref = "run1", eXp = "run2", maxFdrGlobal = 0.05, fitType = "linear")
  expect_equal(getRSE(globalFit, "linear"), 30.12705, tolerance = 1e-05)
})


test_that("test_getGlobalFits", {
  dataPath <- system.file("extdata", package = "DIAlignR")
  fileInfo <- getRunNames(dataPath, oswMerged = TRUE)
  features <- getFeatures(fileInfo, maxFdrQuery = 0.05)
  refRun <- data.frame("transition_group_id" = c(32L, 18342L),
                       "run" = c("run0", "run1"),stringsAsFactors = FALSE)
  globalFits <- getGlobalFits(refRun, features, fileInfo, "loess", 0.05, 0.1)
  globalFit <- globalFits[["run1_run2"]]
  expect_equal(getRSE(globalFit, "loess"), 22.00103, tolerance = 1e-05)

  globalFits <- getGlobalFits(refRun, features, fileInfo, "linear", 0.05, 0.1)
  globalFit <- globalFits[["run1_run2"]]
  expect_equal(getRSE(globalFit, "linear"), 30.12705, tolerance = 1e-05)
})
