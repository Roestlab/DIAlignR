context("Fetching analytes names.")

test_that("test_getAnalytesName",{
  df1 <- data.frame("transition_group_id" = rep("KLYAGAILEV_2", 2),
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
  df1 <- dplyr::as_tibble(df1)
  df2 <- data.frame("transition_group_id" = c("AQPPVSTEY_2", "AQPVPAHVY_2"),
                    "filename" = rep("HLA-Ligand-Atlas/BD-ZH12_Spleen_Class-1/dia_files/170407_AM_BD-ZH12_Spleen_W_10%_DIA_#2_400-650mz_msms42.mzML", 2),
                    "RT" = c(2646.83, 2288.47),
                    "delta_rt" = c(39.65983, -21.74987),
                    "assay_RT" = c(2591.269, 2291.749),
                    "Intensity" = c(3547210, 3160430),
                    "leftWidth" = c(2635.629, 2272.900),
                    "rightWidth" = c(2658.890, 2297.921),
                    "peak_group_rank" = c(1L, 1L),
                    "m_score" = c(0.003287926, 0.003287926),
                    "chromatogramIndex" = c("38090,38804,43646,46675,49639", "37861,38083,43271,49367,49907"),
                    "transition_ids" = c("19613,19614,19619,19623,19627", "17295,17296,17299,17302,17304"),
                    stringsAsFactors=FALSE)
  df2 <- dplyr::as_tibble(df2)
  oswFiles <- list("run0" = df1, "run1" = df2)
  outData <- getAnalytesName(oswFiles, analyteFDR = 0.01, commonAnalytes = TRUE)
  expect_identical(outData, NULL)
  outData <- getAnalytesName(oswFiles, analyteFDR = 0.01, commonAnalytes = FALSE)
  expData <- c("AQPPVSTEY_2", "AQPVPAHVY_2", "KLYAGAILEV_2")
  expect_identical(outData, expData)
  outData <- getAnalytesName(oswFiles, analyteFDR = 0.001, commonAnalytes = FALSE)
  expect_identical(outData, NULL)
})


test_that("test_getAnalytes",{
  dataPath <- "../../data/testData2"
  runs <- c("run1" = "170407_AM_BD-ZH12_Spleen_W_10%_DIA_#2_400-650mz_msms42",
            "run2" = "170413_AM_BD-ZH12_BoneMarrow_W_10%_DIA_#2_400-650mz_msms35")
  outData <- getAnalytes(dataPath, runs, oswMerged = FALSE, maxFdrQuery = 0.001076, commonAnalytes = TRUE)
  expect_identical(tail(outData,3), c("YTAQPTQGY_2", "YTLDQTYAK_2", "YVHDDGRVSY_2"))
  expect_identical(head(outData,3), c("AAAAAAQSVY_2", "AAAGPVADLF_2", "AAFPGASLY_2"))
  expect_identical(length(outData), 776L)
})
