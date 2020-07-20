context("Progressive alignment")

test_that("test_progAlignRuns", {
  skip_if_no_pyopenms()
  dataPath <- system.file("extdata", package = "DIAlignR")
  params <- paramsDIAlignR()
  ropenms <- get_ropenms(condaEnv =  envName)
})
