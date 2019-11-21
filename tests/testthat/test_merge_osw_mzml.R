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
  chromHead <- chromatogramIdAsInteger(chromHead)
  expData <- data.frame("chromatogramId" = c(130110L, 154511L),
                        "chromatogramIndex" = c(1L,2L),
                        stringsAsFactors=FALSE)
  expect_identical(expData, chromHead)
})

test_that("test_mergeOswAnalytes_ChromHeader", {

  #oswName <- "../../data/testData2/osw/170407_AM_BD-ZH12_Spleen_W_10%_DIA_#1_400-650mz_msms41.osw"
  #mzmlName <- "../../data/testData2/mzml/170407_AM_BD-ZH12_Spleen_W_10%_DIA_#1_400-650mz_msms41.chrom.mzML"
  peptides <- c()
  oswAnalytes <- data.frame("transition_group_id" = rep("KLYAGAILEV_2", 10),
                          "filename" = rep("HLA-Ligand-Atlas/BD-ZH12_Spleen_Class-1/dia_files/170407_AM_BD-ZH12_Spleen_W_10%_DIA_#1_400-650mz_msms41.mzML", 10),
                          "RT" = c(rep(4390.35, 5), rep(4433.38, 5)),
                          "delta_rt" = c(rep(66.94013, 5), rep(109.96858, 5)),
                          "assay_RT" = rep(4326.335, 10),
                          "Intensity" = c(rep(6644520, 5), rep(914820, 5)),
                          "leftWidth" = c(rep(4369.325, 5), rep(4418.935, 5)),
                          "rightWidth" = c(rep(4413.365, 5), rep(4439.139, 5)),
                          "peak_group_rank" = c(rep(1, 5), rep(2, 5)),
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
                        "peak_group_rank" = c(1, 2),
                        "m_score" = c(0.002128436, 0.029838394),
                        "chromatogramIndex" = c("100743,104255,107437,109555,114846", "100743,104255,107437,109555,114846"),
                        "transition_ids" = c("45085,45089,45095,45098,45103", "45085,45089,45095,45098,45103"),
                        stringsAsFactors=FALSE)
  outData <- dplyr::as_tibble(outData)
  oswAnalytes <- mergeOswAnalytes_ChromHeader(oswAnalytes, chromHead, analyteFDR = 0.01, runType = "DIA_proteomics")
  expect_identical(outData, oswAnalytes)
  expect_identical(peptides, NULL)
})

test_that("test_getOswFiles", {
  dataPath <- "../../data/testData2/"
  filenames <- data.frame("filename" = c("HLA-Ligand-Atlas/BD-ZH12_Spleen_Class-1/dia_files/170407_AM_BD-ZH12_Spleen_W_10%_DIA_#1_400-650mz_msms41.mzML",
                                         "HLA-Ligand-Atlas/BD-ZH12_Spleen_Class-1/dia_files/170407_AM_BD-ZH12_Spleen_W_10%_DIA_#2_400-650mz_msms42.mzML",
                                         "HLA-Ligand-Atlas/BD-ZH12_BoneMarrow_Class-1/dia_files/170413_AM_BD-ZH12_BoneMarrow_W_10%_DIA_#2_400-650mz_msms35.mzML"),
                          "runs" = c("170407_AM_BD-ZH12_Spleen_W_10%_DIA_#1_400-650mz_msms41",
                                     "170407_AM_BD-ZH12_Spleen_W_10%_DIA_#2_400-650mz_msms42",
                                     "170413_AM_BD-ZH12_BoneMarrow_W_10%_DIA_#2_400-650mz_msms35"),
                          row.names = c("run0", "run1", "run2"),
                          stringsAsFactors=FALSE)
  maxFdrQuery <- 0.05
  analyteFDR <- 0.01
  oswMerged <- FALSE
  outData <- getOswFiles(dataPath, filenames, maxFdrQuery, analyteFDR, oswMerged)
  expect_identical(length(outData), 3L)
  expect_identical(dim(outData[[1]]), c(4997L, 12L))
  expect_identical(dim(outData[[2]]), c(5090L, 12L))
  expect_identical(dim(outData[[3]]), c(3338L, 12L))
  expData <- data.frame("transition_group_id" = "AAAAAAQSVY_2",
                        "filename" = "HLA-Ligand-Atlas/BD-ZH12_Spleen_Class-1/dia_files/170407_AM_BD-ZH12_Spleen_W_10%_DIA_#2_400-650mz_msms42.mzML",
                        "RT" = 2215.56,
                        "delta_rt" = 24.45302,
                        "assay_RT" = 2171.316,
                        "Intensity" = 62024300,
                        "leftWidth" = 2197.681,
                        "rightWidth" = 2231.702,
                        "peak_group_rank" = 1L,
                        "m_score" = 0.0004164053,
                        "chromatogramIndex" = "16487,17740,19178,20625,22078",
                        "transition_ids" = "6419,6425,6427,6429,6434",
                        stringsAsFactors=FALSE)
  expect_equal(as.data.frame(outData[[2]][1,]), expData, tolerance=1e-5)
})
