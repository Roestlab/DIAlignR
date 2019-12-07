#' Get mzML filenames from osw RUN table.
#'
#' @return A dataframe with single column.
filenamesFromOSW <- function(dataPath, pattern){
  # Fetch mzML filenames from RUN table.
  query <- "SELECT DISTINCT RUN.FILENAME AS filename FROM RUN"
  if(pattern == "*.osw"){
    message("Looking for .osw files.")
    # Look for .osw files in osw/ directory.
    temp <- list.files(path = file.path(dataPath, "osw"), pattern="*.osw")
    # Throw an error if no .osw files are found.
    if(length(temp) == 0){return(stop("No .osw files are found."))}
    filenames <- sapply(seq_along(temp), function(i){
      con <- DBI::dbConnect(RSQLite::SQLite(), dbname = file.path(dataPath, "osw", temp[i]))
      # Fetch mzML filenames from RUN table.
      tryCatch(expr = DBI::dbGetQuery(con, statement = query), finally = DBI::dbDisconnect(con))
    })
    filenames <- as.data.frame(unique(unlist(filenames)))
    colnames(filenames) <-  c("filename")
    # Convert filename column from factor to character
    filenames$filename <- as.character(filenames$filename)
    message(nrow(filenames), " .osw files are found.")
  } else if (pattern == "*merged.osw") {
    message("Looking for merged.osw file.")
    # Look for merged.osw files in osw/ directory.
    temp <- list.files(path = file.path(dataPath, "osw"), pattern="*merged.osw")
    # Throw an error if no merged.osw files are found.
    if(length(temp) == 0){return(stop("No merged.osw file is found."))}
    con <- DBI::dbConnect(RSQLite::SQLite(), dbname = file.path(dataPath, "osw", temp[1]))
    # Fetch mzML filenames from RUN table.
    filenames <- tryCatch(expr = DBI::dbGetQuery(con, statement = query), finally = DBI::dbDisconnect(con))
    message(nrow(filenames), " are in ", temp[1], " file")
  } else {
    message("Only .osw and merged.osw files can be read.")
    filenames <- NULL
  }
  filenames
}

#' Get mzML filenames from the directory.
#'
#' @return A dataframe with single column.
filenamesFromMZML <- function(dataPath){
  temp <- list.files(path = file.path(dataPath, "mzml"), pattern="*.chrom.mzML")
  message(length(temp), " .chrom.mzML files are found.")
  mzMLfiles <- sapply(temp, function(x) strsplit(x, split = ".chrom.mzML")[[1]][1])
  mzMLfiles
}

#' Get names of all runs
#'
#' @return A vector with names of all runs.
#' @export
getRunNames <- function(dataPath, oswMerged = TRUE, nameCutPattern = "(.*)(/)(.*)"){
  # Get filenames from RUN table of osw files.
  if(oswMerged == FALSE){
    filenames <- filenamesFromOSW(dataPath, pattern = "*.osw")
  } else{
    filenames <- filenamesFromOSW(dataPath, pattern = "*merged.osw")
  }
  # Get names of mzml files.
  runs <- sapply(filenames[,"filename"], function(x) gsub(nameCutPattern, replacement = "\\3", x))
  fileExtn <- strsplit(runs[[1]], "\\.")[[1]][2]
  fileExtn <- paste0(".", fileExtn)
  filenames$runs <- sapply(runs, function(x) strsplit(x, split = fileExtn)[[1]][1])

  mzMLfiles <- filenamesFromMZML(dataPath)
  # Check if osw files have corresponding mzML file.
  runs <- intersect(filenames$runs, mzMLfiles)
  if(length(runs) != length(filenames$runs)){
    cat("Following files did not have their counterpart in mzml directory\n")
    print(setdiff(filenames$runs, mzMLfiles))
  }
  if(length(runs) == 0){
    message("Names in RUN table of osw files aren't matching to mzML filenames.")
    message("Check if you have correct file names.")
    return(stop("Name mismatch between osw and mzML files."))
  }
  filenames <- filenames[filenames$runs %in% runs,]
  rownames(filenames) <- paste0("run", 0:(nrow(filenames)-1), "")
  filenames
}
