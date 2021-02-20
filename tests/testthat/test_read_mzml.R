context("Read xic files.")

test_that("test_readMzMLHeader",{
  mzmlName <- file.path(system.file("extdata", package = "DIAlignR"), "xics", "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML")
  expOutput <- data.frame("chromatogramId" = c("103114", "103115", "103116"),
                          "chromatogramIndex" = c(1L,2L,3L),
                          "polarity" = c(-1L,-1L,-1L),
                          "precursorIsolationWindowTargetMZ" = c(696.037, 696.037, 696.037),
                          "precursorIsolationWindowLowerOffset" = c(0,0,0),
                          "precursorIsolationWindowUpperOffset" = c(0,0,0),
                          "precursorCollisionEnergy" = c(NA_real_, NA_real_, NA_real_),
                          "productIsolationWindowTargetMZ" = c(717.399, 558.277, 1031.560),
                          "productIsolationWindowLowerOffset" = c(0,0,0),
                          "productIsolationWindowUpperOffset" = c(0,0,0),
                          stringsAsFactors=FALSE)
  expect_identical(dim(readMzMLHeader(mzmlName)), c(72L, 10L))
  expect_equal(readMzMLHeader(mzmlName)[1:3,1:3], expOutput[,1:3], tolerance=1e-6)

  mzmlName <- "getError"
  expect_error(readMzMLHeader(mzmlName))
})

test_that("test_readSqMassHeader",{
  sqName <- file.path(system.file("extdata", package = "DIAlignR"), "xics", "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt.chrom.sqMass")
  expOutput <- data.frame("chromatogramId" = c("103114", "103115", "103116"),
                          "chromatogramIndex" = c(0L,1L,2L),
                          stringsAsFactors=FALSE)
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = sqName)
  expect_identical(dim(readSqMassHeader(con)), c(72L, 2L))
  DBI::dbDisconnect(con)
})

test_that("test_getMZMLpointers",{
  dataPath <- system.file("extdata", package = "DIAlignR")
  fileInfo <- data.frame("chromatogramFile" = c(file.path(dataPath, "xics", "hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt.chrom.mzML"),
                                    file.path(dataPath, "xics", "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML"),
                                    file.path(dataPath, "xics", "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML")),
             row.names = c("run0", "run1", "run2"),
             stringsAsFactors=FALSE)
  outData <- getMZMLpointers(fileInfo)
  expect_is(outData[["run0"]], "mzRpwiz")
  expect_is(outData[["run1"]], "mzRpwiz")
  expect_is(outData[["run2"]], "mzRpwiz")

  fileInfo <- data.frame("chromatogramFile" = c(file.path(dataPath, "xics", "hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt.chrom.sqMass"),
                                                file.path(dataPath, "xics", "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt.chrom.sqMass"),
                                                file.path(dataPath, "xics", "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.chrom.sqMass")),
              row.names = c("run0", "run1", "run2"),
              stringsAsFactors=FALSE)
  outData <- getMZMLpointers(fileInfo)
  expect_is(outData[["run0"]], "SQLiteConnection")
  expect_is(outData[["run1"]], "SQLiteConnection")
  expect_is(outData[["run2"]], "SQLiteConnection")

  for(mz in outData) DBI::dbDisconnect(mz)

})
