context("create sqMass file")

test_that("test_createSqMass",{
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR)
  XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt"]]
  transitionIDs <- list(c(35L, 36L, 37L, 38L, 39L, 410L))
  filename <- "temp.chrom.sqMass"
  createSqMass(filename, XICs, transitionIDs, lossy = TRUE)
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = filename)
  outData <- extractXIC_group2(con, 0:5)
  DBI::dbDisconnect(con)
  file.remove(filename)
  for (i in seq_along(outData)){
    expect_equal(outData[[i]][,1], XICs[[1]][[i]][,1], tolerance = 1e-03) # Time
    expect_equal(outData[[i]][,2], XICs[[1]][[i]][,2], tolerance = 1e-03) # Intensity
  }
})

test_that("test_blobXICs",{
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR)
  XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt"]][["4618"]]
  nativeIds <- 27706:27711
  outData <- blobXICs(XICs, nativeIds, lossy = TRUE)

  for(i in seq(2,12,by = 2)){
    expect_equal(uncompressVec(outData[i-1,1][[1]], 5L), XICs[[i/2]][,1], tolerance = 1e-03)
    expect_equal(uncompressVec(outData[i,1][[1]], 6L), XICs[[i/2]][,2], tolerance = 1e-03)
  }
  expect_identical(colnames(outData), c("DATA", "NATIVE_ID", "DATA_TYPE", "COMPRESSION"))
  expect_identical(outData[,2], rep(as.character(nativeIds), each = 2))
  expect_identical(outData[,3], rep.int(c(2L, 1L), times = 6))
  expect_identical(outData[,4], rep.int(c(5L, 6L), times = 6))

  outData <- blobXICs(XICs, nativeIds, lossy = FALSE)
  for(i in seq(2,12,by = 2)){
    expect_equal(uncompressVec(outData[i-1,1][[1]], 1L), XICs[[i/2]][,1], tolerance = 1e-03)
    expect_equal(uncompressVec(outData[i,1][[1]], 1L), XICs[[i/2]][,2], tolerance = 1e-03)
  }
  expect_identical(colnames(outData), c("DATA", "NATIVE_ID", "DATA_TYPE", "COMPRESSION"))
  expect_identical(outData[,2], rep(as.character(nativeIds), each = 2))
  expect_identical(outData[,3], rep.int(c(2L, 1L), times = 6))
  expect_identical(outData[,4], rep.int(c(1L, 1L), times = 6))
})
