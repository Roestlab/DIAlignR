context("Global Fit")

test_that("test_getLOESSfit", {
  data(oswFiles_DIAlignR, package="DIAlignR")
  oswFiles <- oswFiles_DIAlignR
  Loess.fit <- getLOESSfit(oswFiles, ref = "run2", eXp = "run0", maxFdrLoess = 0.05, spanvalue = 0.1)
  # Testing for Residual standard error
  expect_equal(Loess.fit$s, 25.8761, tolerance = 1e-05)
  # Add predict function as well.
})

test_that("test_getLinearfit", {
  data(oswFiles_DIAlignR, package="DIAlignR")
  oswFiles <- oswFiles_DIAlignR
  Linear.fit <- getLinearfit(oswFiles, ref = "run2", eXp = "run0", maxFdrLoess = 0.05)
  # Testing for Residual standard error
  expect_equal(summary(Linear.fit)[["sigma"]], 26.58681, tolerance = 1e-05)
  # Add predict function as well.
})
