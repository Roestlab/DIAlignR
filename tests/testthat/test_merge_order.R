context("Hierarchical clustering.")

test_that("test_nrDesc", {
  tree <- ape::read.tree(text = "((run1:0.5,run2:0.5)master2:0.5,(run3:0.5,run4:0.5)master3:0.5)master1;")
  outData <- nrDesc(tree)
  expect_identical(outData, c(1, 1, 1, 1, 4, 2, 2))
})

test_that("test_getTree", {
  m <- matrix(c(0,1,2,3, 1,0,1.5,1.5, 2,1.5,0,1, 3,1.5,1,0), byrow = TRUE,
              ncol = 4, dimnames = list(c("run1", "run2", "run3", "run4"),
                                        c("run1", "run2", "run3", "run4")))
  distMat <- as.dist(m, diag = FALSE, upper = FALSE)
  expect_message(outData <- getTree(distMat))
  expect_equal(outData,
          ape::read.tree(text = "((run1:0.5,run2:0.5)master2:0.5,(run3:0.5,run4:0.5)master3:0.5)master1;")
               )
})

test_that("test_getNodeIDs", {
  tree <- ape::read.tree(text = "((run1:0.5,run2:0.5)master2:0.5,(run3:0.5,run4:0.5)master3:0.5)master1;")
  outData <- getNodeIDs(tree)
  expData <- c("run1" = 1L, "run2" = 2L, "run3" = 3L, "run4" = 4L,
               "master1" = 5L, "master2" = 6L, "master3" = 7L)
  expect_identical(outData, expData)
})

test_that("test_traverseUp", {
  skip_if_no_pyopenms()
  dataPath <- system.file("extdata", package = "DIAlignR")
  params <- paramsDIAlignR()
  params[["keepFlanks"]] <- TRUE
  params[["XICfilter"]] <- "none"
  params[["globalAlignmentFdr"]] <- 0.05
  fileInfo <- getRunNames(dataPath = dataPath)
  mzPntrs <- list2env(getMZMLpointers(fileInfo))
  features <- list2env(getFeatures(fileInfo, maxFdrQuery = 0.05, runType = "DIA_proteomics"))
  precursors <- data.frame(transition_group_id = 4618L, peptide_id = 14383L,
                           sequence = "QFNNTDIVLLEDFQK", charge = 3L,
                           group_label = "14299_QFNNTDIVLLEDFQK/3",
                           transition_ids	= I(list(27706:27711)))
  prec2chromIndex <- list2env(getChromatogramIndices(fileInfo, precursors, mzPntrs))
  adaptiveRTs <- new.env()
  refRuns <- new.env()

  tree <- ape::read.tree(text = "(run1:7,run2:2)master1;")
  tree <- ape::reorder.phylo(tree, "postorder")
  ropenms <- get_ropenms(condaEnv = envName, useConda=TRUE)
  m <- capture_messages(traverseUp(tree, dataPath, fileInfo, features, mzPntrs, prec2chromIndex, precursors,
                                        params, adaptiveRTs, refRuns, ropenms))
  expect_equal(m, c("run1 + run2 = master1\n",
                    "Created a child run: master1\n",
                    "Created all master runs.\n"))

  expect_setequal(ls(mzPntrs), c("run0", "run1", "run2", "master1"))
  expect_is(mzPntrs[["master1"]], "mzRpwiz")
  expect_equal(features$master1,  data.frame(transition_group_id = 4618L,
                    feature_id = bit64::as.integer64(7675762503084486466),
                    RT = 5237.8, intensity = 229.707813, leftWidth = 5217.35, rightWidth = 5261.7,
                    peak_group_rank = 1L, m_score = 5.692e-05), tolerance = 1e-04)
  expect_identical(fileInfo["master1", "chromatogramFile"], file.path(dataPath, "mzml", "master1.chrom.mzML"))
  expect_identical(fileInfo["master1", "runName"], "master1")
  expect_identical(prec2chromIndex$master1[,"transition_group_id"], 4618L)
  expect_identical(prec2chromIndex$master1[,"chromatogramIndex"][[1]], 1:6)
  expect_equal(adaptiveRTs[["run1_run2"]], 77.82315, tolerance = 1e-04)
  expect_equal(adaptiveRTs[["run2_run1"]], 70.4146, tolerance = 1e-04)
  expect_identical(refRuns[["master1"]], 1L)

  data(masterXICs_DIAlignR, package="DIAlignR")
  outData <- mzR::chromatograms(mzR::openMSfile(file.path(dataPath, "mzml", "master1.chrom.mzML"), backend = "pwiz"))
  for(i in seq_along(outData)){
    expect_equal(outData[[i]][[1]], masterXICs_DIAlignR[[1]][[i]][[1]], tolerance = 1e-04)
    expect_equal(outData[[i]][[2]], masterXICs_DIAlignR[[1]][[i]][[2]], tolerance = 1e-04)
  }
  outData <- readRDS(file.path(dataPath, "master1_av.rds"), refhook = NULL)
  expect_equal(outData[[1]], masterXICs_DIAlignR[[2]], tolerance = 1e-04)
  file.remove(file.path(dataPath, "master1_av.rds"))
  file.remove(file.path(dataPath, "mzml", "master1.chrom.mzML"))
})

test_that("test_traverseDown", {
  skip_if_no_pyopenms()
  dataPath <- system.file("extdata", package = "DIAlignR")
  params <- paramsDIAlignR()
  params[["keepFlanks"]] <- TRUE
  params[["XICfilter"]] <- "none"
  params[["globalAlignmentFdr"]] <- 0.05
  fileInfo <- getRunNames(dataPath = dataPath)
  mzPntrs <- list2env(getMZMLpointers(fileInfo))
  features <- list2env(getFeatures(fileInfo, maxFdrQuery = 0.05, runType = "DIA_proteomics"))
  precursors <- data.frame(transition_group_id = 4618L, peptide_id = 14383L,
                           sequence = "QFNNTDIVLLEDFQK", charge = 3L,
                           group_label = "14299_QFNNTDIVLLEDFQK/3",
                           transition_ids	= I(list(27706:27711)))
  prec2chromIndex <- list2env(getChromatogramIndices(fileInfo, precursors, mzPntrs))
  adaptiveRTs <- new.env()
  refRuns <- new.env()

  tree <- ape::read.tree(text = "(run1:7,run2:2)master1;")
  tree <- ape::reorder.phylo(tree, "postorder")
  ropenms <- get_ropenms(condaEnv = envName, useConda=TRUE)
  traverseUp(tree, dataPath, fileInfo, features, mzPntrs, prec2chromIndex, precursors,
             params, adaptiveRTs, refRuns, ropenms)
  multipeptide <- list2env(getMultipeptide(precursors, features), hash = TRUE)
  df <- multipeptide[["4618"]]
  analytes <- precursors$transition_group_id
  expect_message(traverseDown(tree, dataPath, fileInfo, multipeptide, prec2chromIndex, mzPntrs, analytes,
               adaptiveRTs, refRuns, params),
               ("Mapping peaks from master1 to run1 and run2.\n|run1 has been aligned to master1.\n|run2 has been aligned to master1.\n|master1 run has been propagated to all parents."),
               all = TRUE)
  # Mapping peaks from master1 to run1 and run2.
  # run1 has been aligned to master1.
  # run2 has been aligned to master1.
  # master1 run has been propagated to all parents.
  df$alignment_rank <- 1L; df$alignment_rank[df$run == "run0"] <- NA_integer_
  expect_equal(multipeptide[["4618"]], df)
  file.remove(file.path(dataPath, "master1_av.rds"))
  file.remove(file.path(dataPath, "mzml", "master1.chrom.mzML"))
})

test_that("test_alignToMaster", {
  skip_if_no_pyopenms()
  dataPath <- system.file("extdata", package = "DIAlignR")
  params <- paramsDIAlignR()
  params[["keepFlanks"]] <- TRUE
  params[["XICfilter"]] <- "none"
  params[["globalAlignmentFdr"]] <- 0.05
  fileInfo <- getRunNames(dataPath = dataPath)
  mzPntrs <- list2env(getMZMLpointers(fileInfo))
  features <- list2env(getFeatures(fileInfo, maxFdrQuery = 0.05, runType = "DIA_proteomics"))
  precursors <- data.frame(transition_group_id = 4618L, peptide_id = 14383L,
                           sequence = "QFNNTDIVLLEDFQK", charge = 3L,
                           group_label = "14299_QFNNTDIVLLEDFQK/3",
                           transition_ids	= I(list(27706:27711)))
  prec2chromIndex <- list2env(getChromatogramIndices(fileInfo, precursors, mzPntrs))
  adaptiveRTs <- new.env()
  refRuns <- new.env()
  tree <- ape::reorder.phylo(ape::read.tree(text = "(run1:7,run2:2)master1;"), "postorder")

  ropenms <- get_ropenms(condaEnv = envName, useConda=TRUE)
  traverseUp(tree, dataPath, fileInfo, features, mzPntrs, prec2chromIndex, precursors, params,
    adaptiveRTs, refRuns, ropenms)
  multipeptide <- list2env(getMultipeptide(precursors, features), hash = TRUE)
  analytes <- precursors$transition_group_id
  alignedVecs <- readRDS(file = file.path(dataPath, "master1_av.rds"))
  adaptiveRT <- max(adaptiveRTs[["run1_run2"]], adaptiveRTs[["run2_run1"]])
  multipeptide[["4618"]]$alignment_rank[multipeptide[["4618"]]$run == "master1"] <- 1L
  df <- multipeptide[["4618"]]

  alignToMaster(ref = "master1", eXp = "run1", alignedVecs, 1L, adaptiveRT,
    multipeptide, prec2chromIndex, mzPntrs, fileInfo, analytes, params)
  df$alignment_rank[df$run == "run1"] <- 1L
  expect_equal(multipeptide[["4618"]], df)

  alignToMaster(ref = "master1", eXp = "run2", alignedVecs, 2L, adaptiveRT,
                multipeptide, prec2chromIndex, mzPntrs, fileInfo, analytes, params)
  df$alignment_rank[df$run == "run2"] <- 1L
  expect_equal(multipeptide[["4618"]], df)

  file.remove(file.path(dataPath, "master1_av.rds"))
  file.remove(file.path(dataPath, "mzml", "master1.chrom.mzML"))
})
