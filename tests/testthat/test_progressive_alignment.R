context("Progressive alignment")

test_that("test_progAlignRuns", {
  dataPath <- system.file("extdata", package = "DIAlignR")
  params <- paramsDIAlignR()
  params[["kernelLen"]] <- 9L
  ropenms <- get_ropenms(condaEnv =  envName)
  params[["context"]] <- "experiment-wide"
  params[["globalAlignment"]] <- "linear"
  BiocParallel::register(BiocParallel::MulticoreParam())
  for(kk in c("sqMass", "mzML")){
    params[["chromFile"]] <- kk
    if(params[["chromFile"]] == "mzML") skip_if_no_pyopenms()
    for(fun in c(lapply, BiocParallel::bplapply)){
      expect_warning(progAlignRuns(dataPath, params = params, outFile = "temp.tsv", ropenms = ropenms, applyFun = fun))
      outData <- read.table("temp.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE)
      expData <- read.table("test3.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE)
      expect_identical(dim(outData), dim(expData))
      expect_identical(colnames(outData), colnames(expData))
      expect_identical(outData[["peptide_id"]], expData[["peptide_id"]])
      expect_identical(outData[["precursor"]], expData[["precursor"]])
      expect_identical(outData[["run"]], expData[["run"]])
      for(i in 4:14){
        print(i)
        expect_equal(outData[[i]], expData[[i]], tolerance = 1e-04)
      }
      file.remove("temp.tsv")
      file.remove(file.path(dataPath, "master.merged.osw"))
      file.remove(file.path(dataPath, "multipeptide.rds"))
      file.remove(list.files(dataPath, pattern = "*_av.rds", full.names = TRUE))
      if(params[["chromFile"]] == "mzML"){
        file.remove(list.files(file.path(dataPath, "mzml"), pattern = "^master[0-9]+\\.chrom\\.mzML$", full.names = TRUE))
      } else if(params[["chromFile"]] == "sqMass"){
        file.remove(list.files(file.path(dataPath, "mzml"), pattern = "^master[0-9]+\\.chrom\\.sqMass$", full.names = TRUE))
      }
    }
  }
})
