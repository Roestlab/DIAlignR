#' Get mzML filenames from osw RUN table.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-14
#' @param dataPath (char) path to mzml and osw directory.
#' @param pattern (char) must be either *.osw or *merged.osw .
#' @return A dataframe with three columns:
#' \item{spectraFile}{(string) as mentioned in RUN table of osw files.}
#' \item{spectraFileID}{(string) ID in RUN table of osw files.}
#' \item{featureFile}{(string) Path to the feature file.}
#' @keywords internal
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' \dontrun{
#' filenamesFromOSW(dataPath, "*.osw")
#' filenamesFromOSW(dataPath, "*merged.osw")
#' }
filenamesFromOSW <- function(dataPath, pattern){
  # Fetch mzML filenames from RUN table.
  query <- "SELECT DISTINCT RUN.FILENAME AS filename, RUN.ID AS ID FROM RUN"
  if(pattern == "*.osw$"){
    message("Looking for .osw files.")
    # Look for .osw files in osw/ directory.
    temp <- list.files(path = file.path(dataPath, "osw"), pattern="*.osw$")
    # Throw an error if no .osw files are found.
    if(length(temp) == 0){return(stop("No .osw files are found."))}
    filenames <- lapply(seq_along(temp), function(i){
      oswName <- file.path(dataPath, "osw", temp[i])
      con <- DBI::dbConnect(RSQLite::SQLite(), dbname = oswName)
      # Fetch mzML filenames from RUN table.
      cbind(tryCatch(expr = DBI::dbGetQuery(con, statement = query),
                     finally = DBI::dbDisconnect(con)), oswName)
    })
    filenames <- do.call(rbind, filenames)
    message(nrow(filenames), " .osw files are found.")
  } else if (pattern == "*merged.osw$") {
    # Look for merged.osw files in osw/ directory.
    message("Looking for merged.osw file.")
    temp <- list.files(path = file.path(dataPath, "osw"), pattern="*merged.osw$")
    # Throw an error if no merged.osw files are found.
    if(length(temp) == 0){return(stop("No merged.osw file is found."))}
    oswName <- file.path(dataPath, "osw", temp[1])
    con <- DBI::dbConnect(RSQLite::SQLite(), dbname = oswName)
    # Fetch mzML filenames from RUN table.
    filenames <- cbind(tryCatch(expr = DBI::dbGetQuery(con, statement = query),
                          finally = DBI::dbDisconnect(con)), oswName)
    message(nrow(filenames), " runs are in ", temp[1], " file")
  } else {
    message("Only .osw and merged.osw files can be read.")
    return(NULL)
  }
  colnames(filenames) <- c("spectraFile", "spectraFileID", "featureFile")
  filenames[["spectraFileID"]] <- as.character(filenames[["spectraFileID"]]) # Convert from integer64 to character.
  filenames[["featureFile"]] <- as.character(filenames[["featureFile"]]) # Convert from factor to character.
  filenames
}

#' Get mzML filenames from the directory.
#'
#' Reads all mzML names avaialble in the directory.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-14
#' @param dataPath (char) Path to mzml and osw directory.
#' @return  A dataframe with two columns:
#' \item{runName}{(string) contain respective mzML names without extension.}
#' \item{chromatogramFile}{(string) Path to the chromatogram file.}
#' @keywords internal
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' \dontrun{
#' filenamesFromMZML(dataPath)
#' }
filenamesFromMZML <- function(dataPath, chromFile){
  if(chromFile == "mzML") p <- ".chrom.mzML$"
  if(chromFile == "sqMass") p <- ".chrom.sqMass$"
  temp <- list.files(path = file.path(dataPath, "mzml"), pattern=p)
  message(length(temp), " ", sub("\\$","",p), " files are found.")
  mzMLfiles <- vapply(temp, function(x) sub(p,"", x), "", USE.NAMES = FALSE)
  output <- data.frame("runName" = mzMLfiles, "chromatogramFile" = file.path(dataPath, "mzml", temp))
  output[["chromatogramFile"]] <- as.character(output[["chromatogramFile"]]) # Convert from factor to character.
  output[["runName"]] <- as.character(output[["runName"]]) # Convert from factor to character.
  output
}

#' Get names of all runs
#'
#' Fetches all osw files, then, keeps only those runs which has corresponding mzML files.
#' mzML file names must match with RUN.FILENAME columns of osw files.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-14
#' @param dataPath (char) Path to mzml and osw directory.
#' @param oswMerged (logical) TRUE for experiment-wide FDR and FALSE for run-specific FDR by pyprophet.
#' @return (dataframe) it has five columns:
#' \item{spectraFile}{(string) as mentioned in RUN table of osw files.}
#' \item{runName}{(string) contain respective mzML names without extension.}
#' \item{spectraFileID}{(string) ID in RUN table of osw files.}
#' \item{featureFile}{(string) Path to the feature file.}
#' \item{chromatogramFile}{(string) Path to the chromatogram file.}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' getRunNames(dataPath = dataPath, oswMerged = TRUE)
#' @export
getRunNames <- function(dataPath, oswMerged = TRUE, params = paramsDIAlignR()){
  # Get filenames from RUN table of osw files.
  if(oswMerged == FALSE){
    filenames <- filenamesFromOSW(dataPath, pattern = "*.osw$")
  } else{
    filenames <- filenamesFromOSW(dataPath, pattern = "*merged.osw$")
  }
  # Get names of mzml files.
  nameCutPattern = "(.*)(/)(.*)" # regex expression to fetch mzML file name from RUN.FILENAME columns of osw files.
  runs <- vapply(filenames[["spectraFile"]], function(x) gsub(nameCutPattern, replacement = "\\3", x), "")
  fileExtn <- strsplit(runs[[1]], "\\.")[[1]][2]
  fileExtn <- paste0(".", fileExtn)
  filenames[["runName"]] <- vapply(runs, function(x) strsplit(x, split = fileExtn)[[1]][1], "")

  mzMLfiles <- filenamesFromMZML(dataPath, params[["chromFile"]])
  # Check if osw files have corresponding mzML file.
  runs <- intersect(filenames[["runName"]], mzMLfiles[["runName"]])
  if(length(runs) != length(filenames[["runName"]])){
    cat("Following files did not have their counterpart in mzml directory\n")
    print(setdiff(filenames[["runName"]], mzMLfiles[["runName"]]))
  }
  if(length(runs) == 0){
    message("Names in RUN table of osw files aren't matching to mzML filenames.")
    message("Check if you have correct file names.")
    return(stop("Name mismatch between osw and mzML files."))
  }
  filenames <- filenames[filenames[["runName"]] %in% runs,]
  filenames <- merge(filenames, mzMLfiles, by = "runName")
  filenames[["chromatogramFile"]] <- as.character(filenames[["chromatogramFile"]])
  filenames[["featureFile"]] <- as.character(filenames[["featureFile"]])
  rownames(filenames) <- paste0("run", 0:(nrow(filenames)-1), "")
  filenames
}


#' Get intersection of runs and fileInfo
#'
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-04-15
#' @param fileInfo (data-frame) output of getRunNames function.
#' @param runs (vector of string) names of mzML files without extension.
#' @return (dataframe) it has five columns:
#' \item{spectraFile}{(string) as mentioned in RUN table of osw files.}
#' \item{runName}{(string) contain respective mzML names without extension.}
#' \item{spectraFileID}{(string) ID in RUN table of osw files.}
#' \item{featureFile}{(string) path to the feature file.}
#' \item{chromatogramFile}{(string) path to the chromatogram file.}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' fileInfo <- getRunNames(dataPath = dataPath, oswMerged = TRUE)
#' runs <- c("hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
#'           "hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt")
#' updateFileInfo(fileInfo, runs)
#' @export
updateFileInfo <- function(fileInfo, runs = NULL){
  if(!is.null(runs)){
    fileInfo <- fileInfo[fileInfo$runName %in% runs,]
    missingRun <- setdiff(runs, fileInfo$runName)
    if(length(missingRun) != 0){
      return(stop(missingRun, " runs are not found."))
    }
  }
  fileInfo
}


addMasterToOSW <- function(dataPath, runs, oswMerged = TRUE){
  df <- data.frame(ID = 1:length(runs), FILENAME = paste(runs, "mzML.gz", sep="."))
  temp <- list.files(path = file.path(dataPath, "osw"), pattern="*merged.osw$", full.names = TRUE)
  newFile <- file.path(dataPath, "master.merged.osw")
  if(file.copy(from = temp, to = newFile)){
    conn <- DBI::dbConnect(RSQLite::SQLite(), newFile)
    DBI::dbExecute(conn,"drop table if exists myTempTable")
    DBI::dbWriteTable(conn,"myTempTable",df)
    DBI::dbExecute(conn,"INSERT INTO RUN (ID, FILENAME) select ID,FILENAME from myTempTable")
    DBI::dbExecute(conn,"drop table if exists myTempTable")
    DBI::dbDisconnect(conn)
  }
}
