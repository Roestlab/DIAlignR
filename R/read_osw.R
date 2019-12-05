#' Get a table of analytes, their chromatogram indices and scores.
#'
#' @return A data-frame.
fetchAnalytesInfo <- function(oswName, maxFdrQuery, oswMerged,
                              analytes, filename, runType, analyteInGroupLabel = FALSE){
  # Establish a connection of SQLite file.
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = oswName)
  # Generate a query.
  query <- getQuery(maxFdrQuery, oswMerged, analytes = analytes,
                    filename = filename, runType = runType,
                    analyteInGroupLabel = analyteInGroupLabel)
  # Run query to get peptides, their coordinates and scores.
  analytesInfo <- tryCatch(expr = DBI::dbGetQuery(con, statement = query),
                           finally = DBI::dbDisconnect(con))
  analytesInfo
}

#' Get a table of analytes, their chromatogram indices and scores.
#'
#' @return A data-frame.
#' @export
getOswAnalytes <- function(dataPath, filenames, oswMerged = TRUE,
                           maxFdrQuery = 0.05, runType  = "DIA_proteomics"){
  oswFiles <- list()
  for(i in 1:nrow(filenames)){
    # Get a query to search against the osw files.
    if(oswMerged == TRUE){
      oswName <- list.files(path = file.path(dataPath, "osw"), pattern="*merged.osw")
      oswName <- file.path(dataPath, "osw", oswName[1])
    } else{
      oswName <- paste0(file.path(dataPath, "osw", filenames$runs[i]), ".osw")
    }
    # Establish a connection of SQLite file.
    con <- DBI::dbConnect(RSQLite::SQLite(), dbname = oswName)
    # Generate a query.
    query <- getAnalytesQuery(maxFdrQuery, oswMerged = oswMerged,
                      filename = filenames$filename[i], runType = runType)
    # Run query to get peptides, their coordinates and scores.
    oswAnalytes <- tryCatch(expr = DBI::dbGetQuery(con, statement = query),
                             finally = DBI::dbDisconnect(con))
    oswFiles[[i]] <- oswAnalytes
  }
  names(oswFiles) <- rownames(filenames)
  oswFiles
}

