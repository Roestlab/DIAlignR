context("Read osw files.")

test_that("test_fetchAnalytesInfo",{
  filenames <- data.frame("filename" = c("data/raw/hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt.mzML.gz",
                                         "data/raw/hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt.mzML.gz",
                                         "data/raw/hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.mzML.gz"),
                          "runs" = c("hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt",
                                     "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
                                     "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"),
                          row.names = c("run0", "run1", "run2"),
                          stringsAsFactors=FALSE)
  oswName <- file.path(system.file("extdata", package = "DIAlignR"), "osw", "merged.osw")
  expOutput <- data.frame("transition_group_id" = rep("19051_KLIVTSEGC[160]FK/2", 6),
                          "filename" = rep("data/raw/hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt.mzML.gz", 6),
                          "RT" = rep(2586.12, 6),
                          "delta_rt" = rep(78.9663, 6),
                          "assay_RT" = rep(13.5, 6),
                          "Intensity" = rep(26.2182, 6),
                          "leftWidth" = rep(2571.738, 6),
                          "rightWidth" = rep(2609.288, 6),
                          "peak_group_rank" = rep(1, 6),
                          "m_score" = rep(0.001041916, 6),
                          "transition_id" = c(58312, 58313, 58314, 58315, 58316, 58317),
                          stringsAsFactors=FALSE)
  outData <- fetchAnalytesInfo(oswName, maxFdrQuery = 0.05, oswMerged = TRUE,
                               analytes = c("19051_KLIVTSEGC[160]FK/2"), filename = filenames$filename[2],
                               runType = "DIA_proteomics", analyteInGroupLabel = TRUE)
  expect_equal(outData, expOutput, tolerance=1e-6)
  outData <- fetchAnalytesInfo(oswName, maxFdrQuery = 0.5, oswMerged = TRUE,
                               analytes = c("IHFLSPVRPFTLTPGDEEESFIQLITPVR_3"), filename = filenames$filename[3],
                               runType = "DIA_proteomics", analyteInGroupLabel = FALSE)
  expOutput <- data.frame("transition_group_id" = rep("IHFLSPVRPFTLTPGDEEESFIQLITPVR_3", 12),
                          "filename" = rep("data/raw/hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.mzML.gz", 12),
                          "RT" = c(rep(6483.50, 6), rep(6597.54, 6)),
                          "delta_rt" = c(rep(78.8163, 6), rep(192.8560, 6)),
                          "assay_RT" = rep(126.7, 12),
                          "Intensity" = c(rep(61.0299, 6), rep(16.7115, 6)),
                          "leftWidth" = c(rep(6468.855, 6), rep(6574.684, 6)),
                          "rightWidth" =  c(rep(6499.579, 6), rep(6615.649, 6)),
                          "peak_group_rank" = c(rep(1, 6), rep(2, 6)),
                          "m_score" = c(rep(5.692077e-05, 6), rep(3.690986e-01,6)),
                          "transition_id" = rep(c(14843, 14844, 14845, 14846, 14847, 14848), 2),
                          stringsAsFactors=FALSE)
  expect_equal(outData, expOutput, tolerance=1e-6)
})

test_that("test_getOswAnalytes",{
  dataPath <- system.file("extdata", package = "DIAlignR")
  fileInfo <- getRunNames(dataPath, oswMerged = TRUE)
  outData <- getOswAnalytes(fileInfo, oswMerged = TRUE,
                            maxFdrQuery = 0.01, runType  = "DIA_proteomics")
  expData <- data.frame("transition_group_id" = rep("AAMIGGADATSNVR_2", 2),
                        "filename" = rep("data/raw/hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.mzML.gz", 2),
                        "peak_group_rank" = c(1L, 1L),
                        "m_score" = rep(5.692077e-05, 2),
                        "transition_id" = c(81958L, 81959L),
                        stringsAsFactors=FALSE)
  expect_identical(dim(outData[["run0"]]), c(1026L, 5L))
  expect_identical(dim(outData[["run1"]]), c(1152L, 5L))
  expect_identical(dim(outData[["run2"]]), c(1086L, 5L))
  expect_equal(outData[["run2"]][1:2,], expData, tolerance=1e-6)
})

test_that("test_fetchPrecursorsInfo",{
  dataPath <- system.file("extdata", package = "DIAlignR")
  filename <- paste0(dataPath,"/osw/merged.osw")
  outData <- fetchPrecursorsInfo(filename, runType = "DIA_proteomics", NULL,
                                 context = "experiment-wide", maxPeptideFdr = 1.0)
  expData <- data.table("transition_group_id" = 32L,
             "peptide_id" = 7040L,
             "sequence" = "GNNSVYMNNFLNLILQNER",
             "charge" = 3L,
             "group_label" = "10030_GNNSVYMNNFLNLILQNER/3",
             "transition_ids" = list(c(192L, 193L, 194L, 195L, 196L, 197L)))
  expect_identical(outData[108,], expData)
  expect_identical(dim(outData), c(302L, 6L))
})

test_that("test_getPrecursors",{
  dataPath <- system.file("extdata", package = "DIAlignR")
  fileInfo <- data.frame("featureFile" = rep(file.path(dataPath, "osw", "merged.osw"),3),
                         row.names = c("run0", "run1", "run2"),
                         stringsAsFactors=FALSE)
  outData <- getPrecursors(fileInfo, oswMerged = TRUE, runType = "DIA_proteomics",
                           context = "experiment-wide", maxPeptideFdr = 0.05)
  expData <- data.table("transition_group_id" = 32L,
                        "peptide_id" = 7040L,
                        "sequence" = "GNNSVYMNNFLNLILQNER",
                        "charge" = 3L,
                        "group_label" = "10030_GNNSVYMNNFLNLILQNER/3",
                        "transition_ids" = list(c(192L, 193L, 194L, 195L, 196L, 197L)),
                        key = c("peptide_id", "transition_group_id"))
  expect_identical(outData[79,], expData)
  expect_identical(dim(outData), c(234L, 6L))
})

test_that("test_getPrecursorByID",{
  dataPath <- system.file("extdata", package = "DIAlignR")
  fileInfo <- data.frame("featureFile" = rep(file.path(dataPath, "osw", "merged.osw"),3),
                         row.names = c("run0", "run1", "run2"),
                         stringsAsFactors=FALSE)
  outData <- getPrecursorByID(c(32L, 2474L), fileInfo, oswMerged = TRUE, runType = "DIA_proteomics")
  expData <- data.table("transition_group_id" = c(32L, 2474L),
                        "peptide_id" = c(7040L, 8496L),
                        "sequence" = c("GNNSVYMNNFLNLILQNER", "IHFLSPVRPFTLTPGDEEESFIQLITPVR"),
                        "charge" = c(3L, 3L),
                        "group_label" = c("10030_GNNSVYMNNFLNLILQNER/3", "12300_IHFLSPVRPFTLTPGDEEESFIQLITPVR/3"),
                        "transition_ids" =  list(c(192L, 193L, 194L, 195L, 196L, 197L),
                                                 c(14843L, 14844L, 14845L, 14846L, 14847L, 14848L)),
                        key = c("peptide_id", "transition_group_id"))
  expect_identical(outData, expData)
})

test_that("test_fetchFeaturesFromRun",{
  dataPath <- system.file("extdata", package = "DIAlignR")
  fileInfo <- data.frame("featureFile" = rep(file.path(dataPath, "osw", "merged.osw"),3),
                         "spectraFileID" = c("125704171604355508", "6752973645981403097", "2234664662238281994"),
                         row.names = c("run0", "run1", "run2"),
                         stringsAsFactors=FALSE)
  fileInfo$featureFile <- as.factor(fileInfo$featureFile)
  outData <- fetchFeaturesFromRun(fileInfo$featureFile[1], runID = "125704171604355508", maxFdrQuery = 0.05, runType = "DIA_proteomics")
  expData <- data.table("transition_group_id" = 32L, "feature_id" = bit64::as.integer64(484069199212214166),
                        "RT" = 6528.23, "intensity" = 26.7603,
                        "leftWidth" = 6518.602, "rightWidth" = 6535.67,
                        "peak_group_rank" = 1L, "m_score" = 0.0264475,
                        key = "transition_group_id")
  expect_equal(outData[1,], expData, tolerance = 1e-04)
  expect_identical(dim(outData), c(211L, 8L))

  outData <- fetchFeaturesFromRun(fileInfo$featureFile[2], runID = "6752973645981403097", maxFdrQuery = 0.01, runType = "DIA_proteomics")
  expData <- data.table("transition_group_id" = 19954L, "feature_id" = bit64::as.integer64(3189052421957813097),
                        "RT" = 5226.47, "intensity" = 104.944,
                        "leftWidth" = 5215.051, "rightWidth" = 5228.706,
                        "peak_group_rank" = 3L, "m_score" = 0.0009634075,
                        key = "transition_group_id")
  expect_equal(outData[192,], expData, tolerance = 1e-04)
  expect_identical(dim(outData), c(192L, 8L))

  outData <- fetchFeaturesFromRun(fileInfo$featureFile[3], runID = "2234664662238281994", maxFdrQuery = 1.00, runType = "DIA_proteomics")
  expData <- data.table("transition_group_id" = 10918L, "feature_id" = bit64::as.integer64(4248772434739795103),
                        "RT" = 6019.18, "intensity" = 78.4294,
                        "leftWidth" = 6006.667, "rightWidth" = 6044.217,
                        "peak_group_rank" = 3L, "m_score" = 0.3225775,
                        key = "transition_group_id")
  expect_equal(outData[500,], expData, tolerance = 1e-04)
  expect_identical(dim(outData), c(926L, 8L))
})

test_that("test_getFeatures",{
  dataPath <- system.file("extdata", package = "DIAlignR")
  fileInfo <- data.frame("featureFile" = rep(file.path(dataPath, "osw", "merged.osw"),3),
                         "spectraFileID" = c("125704171604355508", "6752973645981403097", "2234664662238281994"),
                         row.names = c("run0", "run1", "run2"),
                         stringsAsFactors=FALSE)
  fileInfo$featureFile <- as.factor(fileInfo$featureFile)
  outData <- getFeatures(fileInfo, maxFdrQuery = 0.05, runType = "DIA_proteomics")
  expect_identical(length(outData), 3L)
  expect_identical(dim(outData[["run1"]]), c(227L, 8L))
})

test_that("test_fetchPeptidesInfo", {
  dataPath <- system.file("extdata", package = "DIAlignR")
  filename <- paste0(dataPath,"/osw/merged.osw")
  outData <- fetchPeptidesInfo(oswName = filename, runType = "DIA_proteomics", context = "experiment-wide")
  expData <- data.table("peptide_id" = c(19046L),
                        "run" = bit64::as.integer64(c(6752973645981403097, 2234664662238281994, 125704171604355508)),
                        "score" = c(7.182150, 7.664316, 7.588328),
                        "pvalue" = 5.603183e-05,
                        "qvalue" = 5.204949e-05)
  expect_identical(dim(outData), c(896L, 5L))
  expect_equal(tail(outData,3), expData, tolerance = 1e-06)

  outData2 <- fetchPeptidesInfo(oswName = filename, runType = "DIA_proteomics", context = "global")
  expData <- data.table("peptide_id" = as.integer(), "run" = as.numeric(),
                        "score" = as.numeric(), "pvalue" = as.numeric(), "qvalue" = as.numeric(),
                        stringsAsFactors = FALSE)
  expect_equal(outData2, expData)
})

test_that("test_getPeptideScores", {
  dataPath <- system.file("extdata", package = "DIAlignR")
  fileInfo <- data.frame("featureFile" = rep(file.path(dataPath, "osw", "merged.osw"),3),
                         "spectraFileID" = c("125704171604355508", "6752973645981403097", "2234664662238281994"),
                         row.names = c("run0", "run1", "run2"),
                         stringsAsFactors=FALSE)
  peptides <- c(7260L, 3L, 4L)
  outData <- getPeptideScores(fileInfo, peptides, oswMerged = TRUE, runType = "DIA_proteomics", context = "experiment-wide")
  expData <- data.table("peptide_id" = c(3L, 4L, rep(7260L, 3)),
                        "run" =c(NA_character_, NA_character_, "run0", "run1", "run2"),
                        "score" = c(NA_real_, NA_real_, 7.779751, 7.404515, 7.324655),
                        "pvalue" = c(NA_real_, NA_real_, rep(5.603183e-05, 3)),
                        "qvalue" = c(NA_real_, NA_real_, rep(5.204949e-05, 3)), key = c("peptide_id"))
  expect_equal(outData, expData, tolerance = 1e-06)

  outData2 <- getPeptideScores(fileInfo, peptides = 7260L, oswMerged = TRUE, runType = "DIA_proteomics", context = "run-specific")
  expect_equal(outData2, data.table("peptide_id" = 7260L,
                                   "run" =NA_character_,
                                   "score" = NA_real_, "pvalue" = NA_real_, "qvalue" = NA_real_, key = "peptide_id"))
})

test_that("test_fetchTransitionsFromRun",{
  dataPath <- system.file("extdata", package = "DIAlignR")
  fileInfo <- data.frame("featureFile" = rep(file.path(dataPath, "osw", "merged.osw"),3),
                         "spectraFileID" = c("125704171604355508", "6752973645981403097", "2234664662238281994"),
                         row.names = c("run0", "run1", "run2"),
                         stringsAsFactors=FALSE)
  outData <- fetchTransitionsFromRun(fileInfo$featureFile[1], runID = "125704171604355508", maxFdrQuery = 0.05, runType = "DIA_proteomics")
  expData <- data.table("transition_group_id" = 32L,
                        "feature_id" = bit64::as.integer64(484069199212214166),
                        "RT" = 6528.23, "intensity" = list(c(10.232500, 0.133768, 9.743950, 0.987916, 4.298210, 1.363980)),
                        "leftWidth" = 6518.602, "rightWidth" = 6535.67, "peak_group_rank" = 1L,
                        "m_score" = 0.0264475, key = c("transition_group_id", "peak_group_rank"))
  expect_equal(outData[1,], expData, tolerance = 1e-04)
  expect_identical(dim(outData), c(211L, 8L))

  outData <- fetchTransitionsFromRun(fileInfo$featureFile[2], runID = "6752973645981403097", maxFdrQuery = 0.01, runType = "DIA_proteomics")
  expData <- data.table("transition_group_id" = 19954L, "feature_id" = bit64::as.integer64(3189052421957813097),
                        "RT" = 5226.47, "intensity" = list(c(41.11890, 19.45290, 12.51970, 11.41050, 8.10003, 12.34190)),
                        "leftWidth" = 5215.051, "rightWidth" = 5228.706, "peak_group_rank" = 3L,
                        "m_score" = 0.0009634075, key = c("transition_group_id", "peak_group_rank"))
  expect_equal(outData[192,], expData, tolerance = 1e-04)
  expect_identical(dim(outData), c(192L, 8L))
})

test_that("test_getTransitions",{
  dataPath <- system.file("extdata", package = "DIAlignR")
  fileInfo <- data.frame("featureFile" = rep(file.path(dataPath, "osw", "merged.osw"),3),
                         "spectraFileID" = c("125704171604355508", "6752973645981403097", "2234664662238281994"),
                         row.names = c("run0", "run1", "run2"),
                         stringsAsFactors=FALSE)
  outData <- getTransitions(fileInfo, maxFdrQuery = 0.05, runType = "DIA_proteomics")
  expect_identical(length(outData), 3L)
  expect_identical(dim(outData[["run1"]]), c(227L, 8L))
})
