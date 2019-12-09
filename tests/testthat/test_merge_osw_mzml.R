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
  dataPath <- "../../data/example/"
  filenames <- data.frame("filename" = c("data/raw/hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt.mzML.gz",
                                         "data/raw/hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt.mzML.gz",
                                         "data/raw/hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.mzML.gz"),
                          "runs" = c("hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt",
                                     "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
                                     "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"),
                          row.names = c("run0", "run1", "run2"),
                          stringsAsFactors=FALSE)
  maxFdrQuery <- 0.05
  analyteFDR <- 0.01
  oswMerged <- TRUE
  outData <- getOswFiles(dataPath, filenames, maxFdrQuery, analyteFDR, oswMerged, analyteInGroupLabel = FALSE)
  expect_identical(length(outData), 3L)
  expect_identical(dim(outData[["run0"]]), c(211L, 12L))
  expect_identical(dim(outData[["run1"]]), c(227L, 12L))
  expect_identical(dim(outData[["run2"]]), c(212L, 12L))
  expData <- data.frame("transition_group_id" = "AAAEMGIDLGQVPGTGPK_3",
                        "filename" = "data/raw/hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt.mzML.gz",
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
  outData <- getOswFiles(dataPath, filenames, maxFdrQuery, analyteFDR, oswMerged, analyteInGroupLabel = TRUE)
  expData <- data.frame("transition_group_id" = "10434_LIPNEAADVYVK/2",
                        "filename" = "data/raw/hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt.mzML.gz",
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
