context("Align DIA runs")

test_that("test_alignTargetedRuns",{
  dataPath <- system.file("extdata", package = "DIAlignR")
  params <- paramsDIAlignR()
  params[["maxPeptideFdr"]] <- 0.05
  params[["XICfilter"]] <- "none"
  params[["globalAlignment"]] <- "loess"
  params[["context"]] <- "experiment-wide"
  params[["chromFile"]] <- "sqMass"
  expect_message(
    alignTargetedRuns(dataPath = dataPath,  outFile = "temp", params = params, oswMerged = TRUE,
                      runs = NULL, applyFun = lapply)
  )
  outData <- read.table("temp.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE)
  expData <- read.table("test.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE)
  expect_identical(dim(outData), dim(expData))
  expect_identical(colnames(outData), colnames(expData))
  expect_identical(outData[["peptide_id"]], expData[["peptide_id"]])
  expect_identical(outData[["precursor"]], expData[["precursor"]])
  expect_identical(outData[["run"]], expData[["run"]])
  for(i in 4:10){
    expect_equal(outData[[i]], expData[[i]], tolerance = 1e-04)
  }
  file.remove("temp.tsv")

  runs <- c("hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt",
            "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt")
  BiocParallel::register(BiocParallel::MulticoreParam())
  params <- paramsDIAlignR()
  params[["maxPeptideFdr"]] <- 0.05
  params[["batchSize"]] <- 10L
  params[["globalAlignment"]] <- "loess"
  params[["context"]] <- "experiment-wide"
  params[["chromFile"]] <- "mzML"
  alignTargetedRuns(dataPath = dataPath,  outFile = "temp", params = params, oswMerged = TRUE,
                      runs = runs, applyFun = lapply)
  outData <- read.table("temp.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE)
  expData <- read.table("test2.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE)
  expect_identical(dim(outData), dim(expData))
  expect_identical(colnames(outData), colnames(expData))
  expect_identical(outData[["peptide_id"]], expData[["peptide_id"]])
  expect_identical(outData[["precursor"]], expData[["precursor"]])
  expect_identical(outData[["run"]], expData[["run"]])
  for(i in 4:14){
    expect_equal(outData[[i]], expData[[i]], tolerance = 1e-04)
  }
  file.remove("temp.tsv")

  dataPath <- system.file("extdata", package = "DIAlignR")
  params <- paramsDIAlignR()
  params[["maxPeptideFdr"]] <- 0.05
  params[["XICfilter"]] <- "none"
  params[["context"]] <- "experiment-wide"
  params[["transitionIntensity"]] <- TRUE
  params[["chromFile"]] <- "mzML"
  params[["globalAlignment"]] <- "linear"
  expect_message(
    alignTargetedRuns(dataPath = dataPath,  outFile = "temp", params = params, oswMerged = TRUE,
                      runs = NULL, applyFun = lapply)
  )
  outData <- read.table("temp.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE)
  expData <- read.table("test.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE)
  expect_identical(dim(outData), dim(expData))
  expect_identical(colnames(outData), colnames(expData))
  expect_identical(outData[["peptide_id"]], expData[["peptide_id"]])
  expect_identical(outData[["precursor"]], expData[["precursor"]])
  expect_identical(outData[["run"]], expData[["run"]])
  x <- sapply(outData[["intensity"]], function(a) sum(as.numeric(strsplit(a, split = ",")[[1]])), USE.NAMES = FALSE)
  expect_equal(x, expData[["intensity"]], tolerance = 1e-04)
  for(i in 6:10){
    expect_equal(outData[[i]], expData[[i]], tolerance = 1e-04)
  }
  file.remove("temp.tsv")

})

test_that("test_getAlignObjs",{
  runs <- c("hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt", "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
  refRun <- "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"
  dataPath <- system.file("extdata", package = "DIAlignR")
  analytes <- c(32L, 898L, 4618L)
  params <- paramsDIAlignR()
  params[["maxPeptideFdr"]] <- 0.05
  params[["globalAlignment"]] <- "loess"
  params[["kernelLen"]] <- 13L
  params[["polyOrd"]] <- 4L
  params[["context"]] <- "experiment-wide"
  params[["chromFile"]] <- "mzML"
  expect_warning(
    outData <- getAlignObjs(analytes, runs, dataPath = dataPath, refRun = refRun,
               oswMerged = TRUE, params = params, objType = "light")
    )
  expData <- testAlignObj()
  expect_equal(outData[[2]][["4618"]][["run1_run2"]][[1]], expData, tolerance = 1e-05)
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR
  expect_equal(outData[[2]][["4618"]][["run1_run2"]][["ref"]],
               lapply(XICs[["hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]], as.matrix), tolerance = 1e-05)
  expect_equal(outData[[2]][["4618"]][["run1_run2"]][["eXp"]],
               lapply(XICs[["hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]], as.matrix), tolerance = 1e-05)
  expData <- data.frame("leftWidth" = 5220.758, "RT" = 5238.35, "rightWidth" = 5261.723)
  expect_equal(as.data.frame(outData[[2]][["4618"]][["run1_run2"]][["peak"]]), expData, tolerance = 1e-05)
  expect_identical(outData[[2]][["32"]], NULL)
  expect_identical(outData[[2]][["898"]], NULL)
})
