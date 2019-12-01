context("Global Fit")

test_that("test_getLOESSfit", {
  data(oswFiles_DIAlignR, package="DIAlignR")
  oswFiles <- oswFiles_DIAlignR
  Loess.fit <- getLOESSfit(oswFiles, ref = "run2", eXp = "run0", maxFdrLoess = 0.05, spanvalue = 0.1)
  # Testing for Residual standard error
  expect_equal(Loess.fit$s, 25.8761, tolerance = 1e-05)
  # Add predict function as well.
  expect_equal(predict(Loess.fit, 4090.9), 4025.508, tolerance = 1e-05)
  expect_equal(predict(Loess.fit, 4688.3), 4641.511, tolerance = 1e-05)
})

test_that("test_getLinearfit", {
  data(oswFiles_DIAlignR, package="DIAlignR")
  oswFiles <- oswFiles_DIAlignR
  Linear.fit <- getLinearfit(oswFiles, ref = "run2", eXp = "run0", maxFdrLoess = 0.05)
  # Testing for Residual standard error
  expect_equal(summary(Linear.fit)[["sigma"]], 26.58681, tolerance = 1e-05)
  outData <- predict(Linear.fit, newdata = data.frame("RT.ref"=c(4090.9)))[[1]]
  expect_equal(outData, 4027.519, tolerance = 1e-05)
  outData <- predict(Linear.fit, newdata = data.frame("RT.ref"=c(4688.3)))[[1]]
  expect_equal(outData, 4616.269, tolerance = 1e-05)
})
