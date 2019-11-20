#' Get a table of analytes, their chromatogram indices and scores.
#'
#' @return A data-frame.
fetchAnalytesInfo <- function(oswName, maxFdrQuery, oswMerged,
                              peptides, filename, runType){
  # Establish a connection of SQLite file.
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = oswName)
  # Generate a query.
  query <- getQuery(maxFdrQuery, oswMerged, peptides = peptides,
                    filename = filename, runType = runType)
  # Run query to get peptides, their coordinates and scores.
  analytesInfo <- tryCatch(expr = DBI::dbGetQuery(con, statement = query),
                           finally = DBI::dbDisconnect(con))
  analytesInfo
}
