
test_that("test_calculateIntensity", {
  time <- c( 2.23095,2.239716667,2.248866667,2.25765,2.266416667,
             2.275566667,2.2847,2.293833333,2.304066667,2.315033333,2.325983333,2.336566667,
             2.3468,2.357016667,2.367283333,2.377183333,2.387083333,2.39735,2.40725,2.4175,
             2.4274,2.4373,2.44755,2.45745,2.4677,2.477966667,2.488216667,2.498516667,2.5084,
             2.5183,2.5282,2.538466667,2.548366667,2.558266667,2.568516667,2.578783333,
             2.588683333,2.59895,2.6092,2.619466667,2.630066667,2.64065,2.65125,2.662116667,
             2.672716667,2.6833,2.6939,2.7045,2.715083333,2.725683333,2.736266667,2.746866667,
             2.757833333,2.768416667,2.779016667,2.789616667,2.8002,2.810116667,2.820033333,
             2.830316667,2.840216667,2.849766667,2.859316667,2.868866667,2.878783333,2.888683333,
             2.898233333,2.907783333,2.916033333,2.924266667,2.93215,2.940383333,2.947933333,
             2.955816667,2.964066667,2.97195,2.979833333,2.987716667,2.995616667,3.003516667,
             3.011416667,3.01895,3.026833333,3.034366667,3.042266667,3.0498,3.05735,3.065233333,
             3.073133333,3.080666667,3.0882,3.095733333,3.103633333,3.111533333,3.119066667,
             3.126966667,3.134866667,3.14275,3.15065,3.15855,3.166433333,3.174333333,3.182233333,
             3.190133333,3.198016667,3.205916667,3.213166667
  )

  intensity <- c(
    1447,2139,1699,755,1258,1070,944,1258,1573,1636,
    1762,1447,1133,1321,1762,1133,1447,2391,692,1636,2957,1321,1573,1196,1258,881,
    1384,2076,1133,1699,1384,692,1636,1133,1573,1825,1510,2391,4342,10382,17618,
    51093,153970,368094,632114,869730,962547,966489,845055,558746,417676,270942,
    184865,101619,59776,44863,31587,24036,20450,20324,11074,9879,10508,7928,7110,
    6733,6481,5726,6921,6670,5537,4971,4719,4782,5097,5789,4279,5411,4530,3524,
    2139,3335,3083,4342,4279,3083,3649,4216,4216,3964,2957,2202,2391,2643,3524,
    2328,2202,3649,2706,3020,3335,2580,2328,2894,3146,2769,2517
  )
  df <- data.frame(time, intensity)
  left <- 2.472833334
  right <- 3.022891666
  outData <- calculateIntensity(list(df), left, right, integrationType = "intensity_sum",
                                baselineType = "base_to_base", fitEMG = FALSE)
  expect_equal(outData, 6645331.33866)

  df <- cbind(time, intensity)
  outData <- calculateIntensity(list(df), left, right, integrationType = "intensity_sum",
                              baselineType = "base_to_base", fitEMG = FALSE)
  expect_equal(outData, 6645331.33866)

  outData <- calculateIntensity(list(df, df), left, right, integrationType = "trapezoid",
                                baselineType = "vertical_division_min", fitEMG = FALSE)
  expect_equal(outData, 2*71063.59368, tolerance = 0.01)

  outData <- calculateIntensity(list(df, df), left, right, integrationType = "trapezoid",
                                baselineType = "vertical_division_min", fitEMG = FALSE, transitionIntensity = TRUE)
  expect_equal(outData, rep(71063.59368, 2), tolerance = 0.01)
})


test_that("test_recalculateIntensity", {
  peakTable <- data.frame(precursor = c(1967L, 1967L, 2474L, 2474L),
                          run = rep(c("hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt",
                      "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"), 2),
                      intensity = c(186.166, 579.832, 47.9525, 3.7413),
                      leftWidth = c(5001.76, 5025.66, 6441.51, 6516.6),
                      rightWidth = c(5076.86, 5121.25, 6475.65, 6554.2), stringsAsFactors = FALSE)
  dataPath <- system.file("extdata", package = "DIAlignR")
  outData <- recalculateIntensity(peakTable, dataPath)
  expOutput <- data.frame(precursor = c(1967L, 1967L, 2474L, 2474L),
                          run = rep(c("hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt",
                                      "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"), 2),
                          intensity = c(136.509, 526.932, 36.099, 3.741),
                          stringsAsFactors=FALSE)
  expect_equal(outData[, "intensity"], expOutput[, "intensity"], tolerance = 0.001)
  expect_identical(outData[,c(1,2)], expOutput[,c(1,2)])
})
