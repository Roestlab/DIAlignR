context("Merge two runs")

childFeatures <- function(){
  df <- data.frame(transition_group_id = 4618L,
                   feature_id = bit64::as.integer64(c("7675762503084486466", "8133647188324775335", "436075290410024368",
                                                      "1709365655702457288", "7174847956052480089", "8455234228721846034")),
                   RT = c(5237.8, 5282.2, 5393.2, 5514.4, 5155.9, 4997.2),
                   intensity = c(229.707813, 0.9945442, 29.670797, 2.622945, 4.2806256, 24.580870),
                   leftWidth = c(5217.35, 5278.8, 5365.9, 5497.3, 5137.975, 4980.1),
                   rightWidth = c(5261.7, 5302.7, 5420.5, 5521.2, 5171.25, 5019.4),
                   peak_group_rank = c(1L, 2L, 3L, 4L, 5L, 6L),
                   m_score = c(5.692e-05, 5.039e-02, 0.201, 0.33, 0.36201, 0.36669))
  df
}

test_that("test_getNodeRun",{
  skip_if_no_pyopenms()
  dir.create("mzml")
  ropenms <- get_ropenms(condaEnv = envName, useConda=TRUE)
  dataPath <- system.file("extdata", package = "DIAlignR")
  params <- paramsDIAlignR()
  params[["maxPeptideFdr"]] <- 0.05
  params[["keepFlanks"]] <- TRUE
  params[["XICfilter"]] <- "none";  params[["kernelLen"]] <- 0L
  params[["globalAlignment"]] <- "loess"
  params[["globalAlignmentFdr"]] <- 0.05
  params[["context"]] <- "experiment-wide"
  params[["chromFile"]] <- "mzML"
  fileInfo <- getRunNames(dataPath = dataPath, params = params)
  mzPntrs <- list2env(getMZMLpointers(fileInfo))
  features <- list2env(getFeatures(fileInfo, maxFdrQuery = 0.05, runType = "DIA_Proteomics"))

  precursors <- getPrecursors(fileInfo, oswMerged = TRUE, params[["runType"]], params[["context"]], params[["maxPeptideFdr"]])
  precursors <- precursors[precursors$peptide_id %in% c("7040", "9861", "14383"),]
  peptideIDs <-  c(7040L, 14383L, 9861L)
  peptideScores <- getPeptideScores(fileInfo, peptides = peptideIDs, TRUE, "DIA_Proteomics", "experiment-wide")
  peptideScores <- lapply(peptideIDs, function(pep) dplyr::filter(peptideScores, .data$peptide_id == pep))
  names(peptideScores) <- as.character(peptideIDs)

  peptideScores <- list2env(peptideScores)
  multipeptide <- getMultipeptide(precursors, features)
  prec2chromIndex <- list2env(getChromatogramIndices(fileInfo, precursors, mzPntrs))
  mergeName <- "temp"
  adaptiveRTs <- new.env()
  refRuns <- new.env()

  expect_warning(multipeptide <- getNodeRun(runA="run1", runB="run2", mergeName, dataPath = ".",
                            fileInfo, features, mzPntrs, prec2chromIndex, precursors, params,
                            adaptiveRTs, refRuns, multipeptide, peptideScores, ropenms),
                 "Chromatogram indices for 7040 are missing.")

  expect_identical(ls(mzPntrs), c("run0", "run1", "run2", "temp"))
  expect_is(mzPntrs[["temp"]], "mzRpwiz")
  expect_equal(features$temp[1,], childFeatures()[1,], tolerance = 1e-04)
  expect_equal(features$temp[c(2,3), c(1,3,4,8)],
               data.frame(transition_group_id = c(9719L, 9720L), RT = 2594.85, intensity = c(14.62899, 20.94305),
                          m_score = c(1.041916e-03, 5.692077e-05), row.names = c(2L, 3L)), tolerance = 1e-04)
  expect_identical(fileInfo["temp", "chromatogramFile"], file.path(".", "mzml", "temp.chrom.mzML"))
  expect_identical(fileInfo["temp", "runName"], "temp")
  expect_identical(prec2chromIndex$temp[,"transition_group_id"], c(32L, 4618L, 9719L, 9720L))
  expect_identical(prec2chromIndex$temp[,"chromatogramIndex"], list(rep(NA_integer_, 6), 1:6, 7:12, 13:18))
  expect_equal(adaptiveRTs[["run1_run2"]], 77.0036, tolerance = 1e-04)
  expect_equal(adaptiveRTs[["run2_run1"]], 76.25354, tolerance = 1e-04)
  expect_identical(refRuns[["temp"]], data.frame("var1" = c(2L,1L,1L), "var2" = c("32","4618","9720")))
  expect_equal(multipeptide[["9861"]][5:6, "m_score"], c(1.041916e-03, 5.692077e-05), tolerance = 1e-04)
  expect_equal(peptideScores[["9861"]][4,"pvalue"], 5.603183e-05, tolerance = 1e-04)
  expect_equal(peptideScores[["14383"]][4,"pvalue"], 5.603183e-05, tolerance = 1e-04)


  data(masterXICs_DIAlignR, package="DIAlignR")
  outData <- mzR::chromatograms(mzR::openMSfile(file.path(".", "mzml", "temp.chrom.mzML"), backend = "pwiz"))
  for(i in 1:6){
    expect_equal(outData[[i]][[1]], masterXICs_DIAlignR[[1]][[i]][[1]], tolerance = 1e-04)
    expect_equal(outData[[i]][[2]], masterXICs_DIAlignR[[1]][[i]][[2]], tolerance = 1e-04)
  }
  outData <- readRDS("temp_av.rds", refhook = NULL)
  for(i in 1:3) expect_equal(outData[[2]][,i], masterXICs_DIAlignR[[2]][,i+2], tolerance = 1e-04)
  file.remove("temp_av.rds")
  file.remove(file.path("mzml", "temp.chrom.mzML"))
  unlink("mzml", recursive = TRUE)
})

test_that("test_getChildFeature",{
  data(masterXICs_DIAlignR, package="DIAlignR")
  newXICs <- masterXICs_DIAlignR
  params <- paramsDIAlignR()
  dataPath <- system.file("extdata", package = "DIAlignR")
  fileInfo <- getRunNames(dataPath = dataPath)
  features <- getFeatures(fileInfo, maxFdrQuery = 1.00, runType = "DIA_Proteomics")
  df.ref <- features$run1[features$run1$transition_group_id == 4618L, ]
  df.eXp <- features$run2[features$run2$transition_group_id == 4618L, ]
  outData <- getChildFeature(newXICs[[1]], newXICs[[2]][,3:5], df.ref, df.eXp, params)

  df <- childFeatures()
  expect_equal(outData, df, tolerance = 1e-04)
})

test_that("test_trfrParentFeature",{
  data(masterXICs_DIAlignR, package="DIAlignR")
  newXICs <- masterXICs_DIAlignR
  timeParent <- newXICs[[2]][, c("tAligned.ref", "alignedChildTime")]
  colnames(timeParent) <- c("tAligned", "alignedChildTime")
  params <- paramsDIAlignR()
  df <- data.frame(transition_group_id = 4618L,
                   feature_id = bit64::as.integer64(1:5),
                   RT = c(5238.35, 5395.79, 5519.05, 5022.30, 5054.13),
                   intensity = c(310.01, 94.8768, 15.8803, 74.6494, 115.591),
                   leftWidth = c(5220.758, 5367.551, 5500.687, 5005.685, 5036.41),
                   rightWidth = c(5261.723, 5422.168, 5524.584, 5046.651, 5077.375),
                   peak_group_rank = c(1L, 2L, 3L, 4L, 5L),
                   m_score = c(5.692e-05, 0.201, 0.33, 0.367, 0.368))
  outData <- trfrParentFeature(newXICs[[1]], timeParent, df, params)
  df[,3:6] <- list(c(5237.8, 5393.2, 5514.4, 4997.2, 5026.2),
                   c(229.707813, 29.670797, 2.622945, 24.580870, 15.672424),
                   c(5217.35, 5365.9, 5497.3, 4980.1, 5010.8),
                   c(5261.7, 5420.5, 5521.2, 5019.4, 5048.4))
  expect_equal(outData, df, tolerance = 1e-04)

  #TODO: Think about a test where missing value may occur
})

test_that("test_getChildXICs",{
  dataPath <- system.file("extdata", package = "DIAlignR")
  params <- paramsDIAlignR()
  params[["chromFile"]] <- "mzML"
  fileInfo <- getRunNames(dataPath = dataPath, params = params)
  features <- getFeatures(fileInfo, maxFdrQuery = 1.00, runType = "DIA_Proteomics")
  mzPntrs <- getMZMLpointers(fileInfo)
  precursors <- data.frame(transition_group_id = 4618L, peptide_id = 14383L,
                           sequence = "QFNNTDIVLLEDFQK", charge = 3L,
                           group_label = "14299_QFNNTDIVLLEDFQK/3",
                           transition_ids	= I(list(27706:27711)))
  peptideScores <- getPeptideScores(fileInfo, peptides = 14383L, TRUE, "DIA_Proteomics", "experiment-wide")

  prec2chromIndex <- getChromatogramIndices(fileInfo, precursors, mzPntrs)
  refRun <- data.frame(2, "4618")
  params <- paramsDIAlignR()
  params[["kernelLen"]] <- 0L
  params[["polyOrd"]] <- 4L
  params[["chromFile"]] <- "mzML"
  params[["keepFlanks"]] <- TRUE
  params[["globalAlignment"]] <- "linear"
  outData <- getChildXICs(runA = "run2", runB="run1", fileInfo, features, mzPntrs, precursors,
               prec2chromIndex, refRun, peptideScores, params)
  rm(mzPntrs)
  expData <- masterXICs_DIAlignR
  expect_identical(names(outData[[1]][[1]]), "4618")
  expect_identical(dim(outData[[2]][[1]]), c(204L, 3L))
  expect_equal(outData[[3]], list(ab = 80.58908, ba = 79.22248), tolerance = 1e-04)

  # Fetch 1st peptide's aligned time vectors
  expect_equal(outData[[2]][[1]][1:6,3], c(4963.05, 4966.45, 4969.85, 4973.25,
                                           4976.65, 4980.05), tolerance = 1e-04)
  expect_equal(outData[[2]][[1]][15:20,3], c(5010.85, 5014.25, 5016.8, 5019.35, 5022.75, 5025.3), tolerance = 1e-02)
  expect_equal(outData[[2]][[1]][50:55,3], c(5106.4, 5108.975, 5111.55, 5114.1, 5116.65, 5120.05), tolerance = 1e-02)
  expect_equal(outData[[2]][[1]][199:204,3], c(5567.25, 5570.7, 5574.1, 5576.65, 5579.2, 5582.7), tolerance = 1e-02)

  # Fetch 1st peptide's precursorID 4618 's 2nd fragment-ion
  expect_equal(outData[[1]][[1]][["4618"]][[2]][1:5,2], c(0.1968677, 1.4567888, 1.1024035, 1.0237436, 0.3937056), tolerance = 1e-02)
  expect_equal(outData[[1]][[1]][["4618"]][[2]][100:105,2], c(0.7480956, 0.7185537, 0.5118452, 0.3445192, 0.5610699, 0.50399), tolerance = 1e-02)
  expect_equal(outData[[1]][[1]][["4618"]][[2]][174:177,2], c(0.3149763, 0.4921656, 0.5807825, 0.3543722), tolerance = 1e-02)
})
