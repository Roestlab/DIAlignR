#' Create an sqMass file
#'
#' Writes a sqMass file having chromatograms and their native IDs.
#'
#' @details
#' - compression is one of 0 = no, 1 = zlib, 2 = np-linear, 3 = np-slof,
#'  4 = np-pic, 5 = np-linear + zlib, 6 = np-slof + zlib, 7 = np-pic + zlib \cr
#' - data_type is one of 0 = mz, 1 = int, 2 = rt \cr
#' - data contains the raw (blob) data for a single data array
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2021) + GPL-3
#' Date: 2021-01-16
#' @import DBI RSQLite
#' @param filename (string) name of the mzML file to be written. Extension should be .chrom.sqMass.
#' @param XICs (list of list of data-frames) list of extracted ion chromatograms of all precursors.
#' @param transitionIDs (list of integer) length must be the same as of XICs.
#' @return (None)
#' @seealso \code{\link{createMZML}, \link{blobXICs}}
#' @examples
#' data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR)
#' XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt"]]
#' XICs <- list(XICs[[1]], XICs[[1]])
#' nativeIds <- list(27706:27711, 1:6)
#' sqName <- "testfile.chrom.sqMass"
#' \dontrun{
#' createSqMass(sqName, XICs, nativeIds)
#' con <- DBI::dbConnect(RSQLite::SQLite(), dbname = sqName)
#' XIC_group <- extractXIC_group2(con, 0:5)
#' DBI::dbDisconnect(con)
#' file.remove(sqName)
#' }
#' @export
createSqMass <- function(filename, XICs, transitionIDs){
  db = DBI::dbConnect(RSQLite::SQLite(), dbname=filename)
  con = DBI::dbConnect(RSQLite::SQLite(), dbname=":memory:")

  query = "CREATE TABLE DATA(SPECTRUM_ID INT, CHROMATOGRAM_ID INT, COMPRESSION INT, DATA_TYPE INT, DATA BLOB NOT NULL);"
  DBI::dbExecute(con, query)
  # chromatogram table
  query = "CREATE TABLE CHROMATOGRAM(ID INT PRIMARY KEY NOT NULL, RUN_ID INT, NATIVE_ID TEXT NOT NULL);"
  DBI::dbExecute(con, query)

  # Convert XICs to compatible format for SQLite
  dfs <- lapply(seq_along(XICs), function(prec){
    if(is.null(XICs[[prec]])) return(NULL) # Skip empty XICs
    blobXICs(XICs[[prec]], transitionIDs[[prec]])
  })
  dfs <- dplyr::bind_rows(dfs)
  n1 <- (nrow(dfs)/2)
  dfs$SPECTRUM_ID <- NA_integer_
  dfs$CHROMATOGRAM_ID <- rep(seq(0, n1-1), each = 2L)

  df <- dfs[, c("SPECTRUM_ID", "CHROMATOGRAM_ID", "COMPRESSION", "DATA_TYPE", "DATA")]
  DBI::dbWriteTable(conn=con, name="DATA", df, append=T, row.names = FALSE)
  idx <- seq.int(from = 1, to = 2*n1, by= 2)
  df <- dfs[idx,c("CHROMATOGRAM_ID", "NATIVE_ID")]
  df$RUN_ID <- 0L
  colnames(df)[1] <- "ID"
  DBI::dbWriteTable(conn=con, name="CHROMATOGRAM", df, append=T, row.names = FALSE)
  # Create indices.
  DBI::dbExecute(con, "CREATE INDEX data_chr_idx ON DATA(CHROMATOGRAM_ID);")

  # Store as sqMass
  RSQLite::sqliteCopyDatabase(con, db)
  DBI::dbDisconnect(db)
  DBI::dbDisconnect(con)
  invisible(NULL)
}


#' Format XICs to blob
#'
#' @details
#' DATA_TYPE is one of 0 = mz, 1 = intensity, 2 = rt \cr
#' COMPRESSION is one of 0 = no, 1 = zlib, 2 = np-linear,
#'  3 = np-slof, 4 = np-pic, 5 = np-linear + zlib,
#'   6 = np-slof + zlib, 7 = np-pic + zlib
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2021) + GPL-3
#' Date: 2021-01-16
#' @import RMSNumpress
#' @param XICs (list) a list of data-frames. Each data frame has elution time and intensity of fragment-ion XIC.
#' @param nativeId (integer) transition ID of the xic.
#' @return (data.frame)
#'
#' @keywords internal
#' @examples
#' data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR)
#' XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt"]][["4618"]]
#' nativeIds <- 27706:27711
#' \dontrun{
#' blobXICs(XICs, nativeIds)
#' }
blobXICs <- function(XICs, nativeIds){
  n1 <- length(XICs)
  # Iterate over each fragment-ion.
  df <- vapply(XICs, function(xic){
    v <- vector(mode = "list", length = 2L)
    v[[1]] <- memCompress(RMSNumpress::encodeLinear(xic[,1], RMSNumpress::optimalLinearFixedPoint(xic[,1])), type = "gzip")
    v[[2]] <- memCompress(RMSNumpress::encodeSlof(xic[,2], RMSNumpress::optimalSlofFixedPoint(xic[,2])), type = "gzip")
    v
  }, vector(mode = "list", length = 2L), USE.NAMES = FALSE)
  df <- as.data.frame(do.call(cbind, list(c(df))))
  colnames(df) <- "DATA"
  df$NATIVE_ID <- rep(as.character(nativeIds), times = 1, each = 2)
  df$DATA_TYPE <- rep(c(2L, 1L), n1) # Time, intensity
  df$COMPRESSION <- rep(c(5L, 6L), n1)
  df
}
