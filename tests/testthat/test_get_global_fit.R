context("Global Fit")

test_that("test_getLOESSfit", {
  data(oswFiles_DIAlignR, package="DIAlignR")
  oswFiles <- oswFiles_DIAlignR
  Loess.fit <- getLOESSfit(oswFiles, ref = "run1", eXp = "run2", maxFdrLoess = 0.05, spanvalue = 0.1)
  # Testing for Residual standard error
  expect_equal(Loess.fit$s, 22.23519, tolerance = 1e-05)
  # Add predict function as well.
  expect_equal(predict(Loess.fit, newdata = data.frame("RT.ref"= 4978.4))[[1]], 4964.752, tolerance = 1e-05)
  expect_equal(predict(Loess.fit, newdata = data.frame("RT.ref"= 5575.8))[[1]], 5565.462, tolerance = 1e-05)
})

test_that("test_getLinearfit", {
  data(oswFiles_DIAlignR, package="DIAlignR")
  oswFiles <- oswFiles_DIAlignR
  Linear.fit <- getLinearfit(oswFiles, ref = "run1", eXp = "run2", maxFdrLoess = 0.05)
  # Testing for Residual standard error
  expect_equal(summary(Linear.fit)[["sigma"]], 30.12705, tolerance = 1e-05)
  outData <- predict(Linear.fit, newdata = data.frame("RT.ref"=4978.4))[[1]]
  expect_equal(outData, 4990.682, tolerance = 1e-05)
  outData <- predict(Linear.fit, newdata = data.frame("RT.ref"=5575.8))[[1]]
  expect_equal(outData, 5577.561, tolerance = 1e-05)
})
