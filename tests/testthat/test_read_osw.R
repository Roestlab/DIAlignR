context("Read osw files.")

test_that("test_fetchAnalytesInfo",{
  filenames <- data.frame("filename" = c("HLA-Ligand-Atlas/BD-ZH12_Spleen_Class-1/dia_files/170407_AM_BD-ZH12_Spleen_W_10%_DIA_#1_400-650mz_msms41.mzML",
                                         "HLA-Ligand-Atlas/BD-ZH12_Spleen_Class-1/dia_files/170407_AM_BD-ZH12_Spleen_W_10%_DIA_#2_400-650mz_msms42.mzML",
                                         "HLA-Ligand-Atlas/BD-ZH12_BoneMarrow_Class-1/dia_files/170413_AM_BD-ZH12_BoneMarrow_W_10%_DIA_#2_400-650mz_msms35.mzML"),
                          "runs" = c("170407_AM_BD-ZH12_Spleen_W_10%_DIA_#1_400-650mz_msms41",
                                     "170407_AM_BD-ZH12_Spleen_W_10%_DIA_#2_400-650mz_msms42",
                                     "170413_AM_BD-ZH12_BoneMarrow_W_10%_DIA_#2_400-650mz_msms35"),
                          row.names = c("run0", "run1", "run2"),
                          stringsAsFactors=FALSE)
  expOutput <- data.frame("transition_group_id" = rep("KLYAGAILEV_2", 10),
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
  oswName <- "../../data/testData2/osw/170407_AM_BD-ZH12_Spleen_W_10%_DIA_#1_400-650mz_msms41.osw"
  outData <- fetchAnalytesInfo(oswName, maxFdrQuery = 0.05, oswMerged = FALSE,
                               analytes = c("KLYAGAILEV_2"), filename = filenames$filename[1],
                               runType = "DIA_proteomics")
  #skip('skip')
  expect_equal(outData, expOutput, tolerance=1e-6)
})
