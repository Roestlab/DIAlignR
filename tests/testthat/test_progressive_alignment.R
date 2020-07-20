context("Progressive alignment")

test_that("test_progAlignRuns", {
  skip_if_no_pyopenms()
  dataPath <- system.file("extdata", package = "DIAlignR")
  params <- paramsDIAlignR()
  ropenms <- get_ropenms(condaEnv =  envName)
  expect_warning(progAlignRuns(dataPath, params = params, outFile = "temp.tsv", ropenms = ropenms))
  outData <- read.table("temp.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE)
  expData <- read.table("test3.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE)
  expect_identical(dim(outData), dim(expData))
  expect_identical(colnames(outData), colnames(expData))
  expect_identical(outData[["peptide"]], expData[["peptide"]])
  expect_identical(outData[["run"]], expData[["run"]])
  for(i in 3:13){
    expect_equal(outData[[i]], expData[[i]], tolerance = 1e-04)
  }
  file.remove("temp.tsv")
  file.remove(list.files(dataPath, pattern = "*_av.rds", full.names = TRUE))
  file.remove(list.files(file.path(dataPath, "mzml"), pattern = "^master[0-9]+\\.chrom\\.mzML$", full.names = TRUE))
})
