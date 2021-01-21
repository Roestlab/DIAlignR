context("Merging OSW and mzML files.")

test_that("test_chromatogramIdAsInteger",{
  chromHead <- data.frame("chromatogramId" = c("130110", "154511"),
                          "chromatogramIndex" = c(1L,2L),
                          "polarity" = c(-1L,-1L),
                          "precursorIsolationWindowTargetMZ" = c(422.7318, 417.9083),
                          "precursorIsolationWindowLowerOffset" = c(0,0),
                          "precursorIsolationWindowUpperOffset" = c(0,0),
                          "precursorCollisionEnergy" = c(NA_real_, NA_real_),
                          "productIsolationWindowTargetMZ" = c(286.1761, 286.6341),
                          "productIsolationWindowLowerOffset" = c(0,0),
                          "productIsolationWindowUpperOffset" = c(0,0),
                          stringsAsFactors=FALSE)
  chromatogramIdAsInteger(chromHead)
  expData <- data.frame("chromatogramId" = c(130110L, 154511L),
                        "chromatogramIndex" = c(1L,2L),
                        stringsAsFactors=FALSE)
  expect_identical(expData, chromHead)
})

test_that("test_mapPrecursorToChromIndices",{
  chromHead <- data.frame("chromatogramId" = c(130110L, 154511L, 102750L, 131399L, 110509L, 153463L, 45085L, 45089L, 45095L, 45098L, 45103L),
                          "chromatogramIndex" = c(1L,2L,3L,4L,5L,6L,100743L,104255L,107437L,109555L,114846L),
                          stringsAsFactors=FALSE)
  prec2transition <- data.frame("transition_group_id" = c(rep(32L,5), rep(192L,5), 396L),
                                "transition_ids" = c(154511L, 102750L, 131399L, 2130110L, 2130120L, 45085L, 45089L, 45095L, 45098L, 45103L, 45104L))
  outData <- mapPrecursorToChromIndices(prec2transition, chromHead)
  expData <- data.frame("transition_group_id" = c(32L, 192L, 396L))
  expData[1, "chromatogramIndex"][[1]] <- list(c(2L, 3L, 4L, NA_integer_,NA_integer_))
  expData[2, "chromatogramIndex"][[1]] <- list(c(100743L, 104255L, 107437L, 109555L, 114846L))
  expData[3, "chromatogramIndex"][[1]] <- list(c(NA_integer_))
  expect_identical(expData, outData)
})

test_that("test_mergeOswAnalytes_ChromHeader", {
  oswAnalytes <- data.frame("transition_group_id" = rep("KLYAGAILEV_2", 10),
                          "filename" = rep("HLA-Ligand-Atlas/BD-ZH12_Spleen_Class-1/dia_files/170407_AM_BD-ZH12_Spleen_W_10%_DIA_#1_400-650mz_msms41.mzML", 10),
                          "RT" = c(rep(4390.35, 5), rep(4433.38, 5)),
                          "delta_rt" = c(rep(66.94013, 5), rep(109.96858, 5)),
                          "assay_RT" = rep(4326.335, 10),
                          "Intensity" = c(rep(6644520, 5), rep(914820, 5)),
                          "leftWidth" = c(rep(4369.325, 5), rep(4418.935, 5)),
                          "rightWidth" = c(rep(4413.365, 5), rep(4439.139, 5)),
                          "peak_group_rank" = c(rep(1L, 5), rep(2L, 5)),
                          "m_score" = c(rep(0.002128436, 5), rep(0.029838394, 5)),
                          "transition_id" = c(45085, 45089, 45095, 45098, 45103, 45085, 45089, 45095, 45098, 45103),
                          stringsAsFactors=FALSE)
  chromHead <- data.frame("chromatogramId" = c(130110L, 154511L, 102750L, 131399L, 110509L, 153463L, 45085L, 45089L, 45095L, 45098L, 45103L),
                          "chromatogramIndex" = c(1L,2L,3L,4L,5L,6L,100743L,104255L,107437L,109555L,114846L),
                          stringsAsFactors=FALSE)
  outData <- data.frame("transition_group_id" = rep("KLYAGAILEV_2", 2),
                        "filename" = rep("HLA-Ligand-Atlas/BD-ZH12_Spleen_Class-1/dia_files/170407_AM_BD-ZH12_Spleen_W_10%_DIA_#1_400-650mz_msms41.mzML", 2),
                        "RT" = c(4390.35, 4433.38),
                        "delta_rt" = c(66.94013, 109.96858),
                        "assay_RT" = rep(4326.335, 2),
                        "Intensity" = c(6644520, 914820),
                        "leftWidth" = c(4369.325, 4418.935),
                        "rightWidth" = c(4413.365, 4439.139),
                        "peak_group_rank" = c(1L, 2L),
                        "m_score" = c(0.002128436, 0.029838394),
                        "chromatogramIndex" = c("100743,104255,107437,109555,114846", "100743,104255,107437,109555,114846"),
                        "transition_ids" = c("45085,45089,45095,45098,45103", "45085,45089,45095,45098,45103"),
                        stringsAsFactors=FALSE)
  outData <- dplyr::as_tibble(outData)
  mergeOswAnalytes_ChromHeader(oswAnalytes, chromHead, analyteFDR = 0.01, runType = "DIA_proteomics")
  expect_identical(outData, oswAnalytes)
})

test_that("test_getOswFiles", {
  dataPath <- system.file("extdata", package = "DIAlignR")
  params <- paramsDIAlignR()
  params[["chromFile"]] <- "mzML"
  fileInfo <- getRunNames(dataPath, oswMerged = TRUE, params)
  mzPntrs <- getMZMLpointers(fileInfo)
  maxFdrQuery <- 0.05
  analyteFDR <- 0.01
  oswMerged <- TRUE
  outData <- getOswFiles(fileInfo, mzPntrs, maxFdrQuery, analyteFDR, oswMerged, analyteInGroupLabel = FALSE)
  rm(mzPntrs)

  expect_identical(length(outData), 3L)
  expect_identical(dim(outData[["run0"]]), c(211L, 11L))
  expect_identical(dim(outData[["run1"]]), c(227L, 11L))
  expect_identical(dim(outData[["run2"]]), c(212L, 11L))
  expData <- data.frame("transition_group_id" = "AAAEMGIDLGQVPGTGPK_3",
                        "RT" = 3884.56,
                        "delta_rt" = 17.0629,
                        "assay_RT" = 53,
                        "Intensity" = 42.0064,
                        "leftWidth" = 3865.266,
                        "rightWidth" = 3902.816,
                        "peak_group_rank" = 1L,
                        "m_score" = 0.0241832,
                        "chromatogramIndex" = "NA,NA,NA,NA,NA,NA",
                        "transition_ids" = "106468,106469,106470,106471,106472,106473",
                        stringsAsFactors=FALSE)
  expect_equal(as.data.frame(outData[["run1"]][1,]), expData, tolerance=1e-5)

  mzPntrs <- getMZMLpointers(fileInfo)
  outData <- getOswFiles(fileInfo, mzPntrs, maxFdrQuery, analyteFDR, oswMerged, analyteInGroupLabel = TRUE)
  rm(mzPntrs)
  expData <- data.frame("transition_group_id" = "10434_LIPNEAADVYVK/2",
                        "RT" = 3239.48,
                        "delta_rt" = 9.10682,
                        "assay_RT" = 34.5,
                        "Intensity" = 2087.8300,
                        "leftWidth" = 3217.045,
                        "rightWidth" = 3264.839,
                        "peak_group_rank" = 1L,
                        "m_score" = 5.692077e-05,
                        "chromatogramIndex" = "NA,NA,NA,NA,NA,NA",
                        "transition_ids" = "2820,2821,2822,2823,2824,2825",
                        stringsAsFactors=FALSE)
  expect_equal(as.data.frame(outData[["run1"]][1,]), expData, tolerance=1e-5)
})

test_that("test_getChromatogramIndices",{
  dataPath <- system.file("extdata", package = "DIAlignR")
  params <- paramsDIAlignR()
  params[["chromFile"]] <- "mzML"
  fileInfo <- getRunNames(dataPath, oswMerged = TRUE, params)
  precursors <- getPrecursors(fileInfo, oswMerged = TRUE,
                              context = "experiment-wide", maxPeptideFdr = 1.00)
  mzPntrs <- getMZMLpointers(fileInfo)
  outData <- getChromatogramIndices(fileInfo, precursors, mzPntrs)
  outData2 <- getChromatogramIndices(fileInfo, precursors[c(3,25,1),], mzPntrs)
  rm(mzPntrs)

  expData <- data.frame("transition_group_id" = c(9720L, 9723L),
                        row.names = c(144L, 145L))
  expData[1, "chromatogramIndex"][[1]] <- list(c(49L, 50L, 51L, 52L, 53L, 54L))
  expData[2, "chromatogramIndex"][[1]] <- list(rep(NA_integer_, 6))
  expect_identical(expData, outData[["run2"]][144:145,])

  expData <- data.frame("transition_group_id" = c(470L, 1967L, 32L))
  expData[1, "chromatogramIndex"][[1]] <- list(rep(NA_integer_, 6))
  expData[2, "chromatogramIndex"][[1]] <- list(c(13:18))
  expData[3, "chromatogramIndex"][[1]] <- list(rep(NA_integer_, 6))
  expect_identical(expData, outData2[["run2"]])

  dataPath <- system.file("extdata", package = "DIAlignR")
  params <- paramsDIAlignR()
  params[["chromFile"]] <- "sqMass"
  fileInfo <- getRunNames(dataPath = dataPath, oswMerged = TRUE, params = params)
  precursors <- getPrecursors(fileInfo, oswMerged = TRUE,
                              context = "experiment-wide", maxPeptideFdr = 1.00)
  mzPntrs <- getMZMLpointers(fileInfo)
  outData <- getChromatogramIndices(fileInfo, precursors, mzPntrs)
  outData2 <- getChromatogramIndices(fileInfo, precursors[c(3,25,1),], mzPntrs)
  for(mz in mzPntrs) DBI::dbDisconnect(mz)

  expData <- data.frame("transition_group_id" = c(9720L, 9723L),
                        row.names = c(144L, 145L))
  expData[1, "chromatogramIndex"][[1]] <- list(c(48L, 49L, 50L, 51L, 52L, 53L))
  expData[2, "chromatogramIndex"][[1]] <- list(rep(NA_integer_, 6))
  expect_identical(expData, outData[["run2"]][144:145,])

  expData <- data.frame("transition_group_id" = c(470L, 1967L, 32L))
  expData[1, "chromatogramIndex"][[1]] <- list(rep(NA_integer_, 6))
  expData[2, "chromatogramIndex"][[1]] <- list(c(12:17))
  expData[3, "chromatogramIndex"][[1]] <- list(rep(NA_integer_, 6))
  expect_identical(expData, outData2[["run2"]])

})
