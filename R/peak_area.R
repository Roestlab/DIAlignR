#' Calculates area of a peak in XIC group
#'
#' Retention time from reference run is mapped to experiment run using AlignObj.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-04-13
#'
#' @param XICs (list) list of extracted ion chromatograms of a precursor.
#' @param left (numeric) left boundary of the peak.
#' @param right (numeric) right boundary of the peak.
#' @param integrationType (string) method to ompute the area of a peak contained in XICs. Must be
#'  from "intensity_sum", "trapezoid", "simpson".
#' @param baselineType (string) method to estimate the background of a peak contained in XICs. Must be
#'  from "base_to_base", "vertical_division_min", "vertical_division_max".
#' @param fitEMG (logical) enable/disable exponentially modified gaussian peak model fitting.
#' @param baseSubtraction (logical) TRUE: remove background from peak signal using estimated noise levels.
#' @param transitionIntensity (logical) TRUE: return intensity of each transition, FALSE: return sum of all transitions.
#' @return (numeric)
#' @keywords internal
#' @seealso \code{\link{getMultipeptide}, \link{setAlignmentRank}}
#' @examples
#' data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
#' XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]]
#' \dontrun{
#' calculateIntensity(XICs, 5220, 5261, integrationType = "intensity_sum",
#'  baselineType = "base_to_base", fitEMG = FALSE)
#' }
calculateIntensity <- function(XICs, left, right, integrationType, baselineType,
                               fitEMG = FALSE, baseSubtraction = TRUE, transitionIntensity = FALSE){
  time <- lapply(XICs, `[`, i =, j =1)
  intensityList <- lapply(XICs, `[`, i =, j= 2)
  intensity <- areaIntegrator(time, intensityList, left, right, integrationType, baselineType,
                              fitEMG, baseSubtraction)
  if(transitionIntensity) return (intensity)
  sum(intensity, na.rm = FALSE)
}


newRow <- function(xics, left, right, RT, analyte, run, params){
  intensity <- calculateIntensity(xics, left, right, params[["integrationType"]], params[["baselineType"]],
                                  params[["fitEMG"]], params[["baseSubtraction"]], params[["transitionIntensity"]])
  row <- data.frame("transition_group_id" = analyte, "feature_id" = bit64::NA_integer64_,
                    "RT" = RT, "intensity"= NA_real_, "leftWidth" = left, "rightWidth" = right,
                    "m_score" = NA_real_, "peak_group_rank" = NA_integer_, "run" = run,
                    "alignment_rank" = 1L)
  if(params[["transitionIntensity"]]){
    row[1, "intensity"][[1]] <- list(intensity)
  } else {
    row[1, "intensity"] <- intensity
  }
  row
}

#' Calculates area of peaks in peakTable
#'
#' For the give peak boundary in peakTable, the function extracts raw chromatograms and recalculate intensities.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-05-28
#'
#' @importFrom magrittr %>%
#' @inheritParams alignTargetedRuns
#' @param peakTable (data-frame) usually an output of alignTargetedRuns. Must have these columns: run, precursor, leftWidth, rightWidth.
#' @param dataPath (string) path to xics and osw directory.
#' @param oswMerged (logical) TRUE for experiment-wide FDR and FALSE for run-specific FDR by pyprophet.
#' @return (data-frame)
#' @seealso \code{\link{alignTargetedRuns}, \link{calculateIntensity}}
#' @examples
#' peakTable <- data.frame(precursor = c(1967L, 1967L, 2474L, 2474L),
#'                    run = rep(c("hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt",
#'                    "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"), 2),
#'                    intensity = c(186.166, 579.832, 47.9525, 3.7413),
#'                    leftWidth = c(5001.76, 5025.66, 6441.51, 6516.6),
#'                    rightWidth = c(5076.86, 5121.25, 6475.65, 6554.2), stringsAsFactors = FALSE)
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' newTable <- recalculateIntensity(peakTable, dataPath)
#' @export
recalculateIntensity <- function(peakTable, dataPath = ".", oswMerged = TRUE, params = paramsDIAlignR()){
  runs <- unique(peakTable$run)
  analytes <- unique(peakTable$precursor)
  fileInfo <- getRunNames(dataPath, oswMerged, params)
  fileInfo <- updateFileInfo(fileInfo, runs)

  ######### Get Precursors from the query and respectve chromatogram indices. ######
  precursors <- getPrecursorByID(analytes, fileInfo)

  ######### Collect pointers for each mzML file. #######
  message("Collecting metadata from mzML files.")
  mzPntrs <- getMZMLpointers(fileInfo)
  message("Metadata is collected from mzML files.")

  ############# Get chromatogram Indices of precursors across all runs. ############
  prec2chromIndex <- getChromatogramIndices(fileInfo, precursors, mzPntrs)

  newArea <- list()
  for (run in rownames(fileInfo)){
    newArea[[run]] <- rep(NA_real_, length(analytes))
    runname <- fileInfo[run, "runName"]
    for (i in seq_along(analytes)){
      analyte <- analytes[i]
      df <- dplyr::filter(peakTable, .data$precursor == analyte, .data$run == runname) %>%
        dplyr::select(.data$leftWidth, .data$rightWidth)
      chromIndices <- prec2chromIndex[[run]][["chromatogramIndex"]][[i]]

      # Get XIC_group from reference run. if missing, go to next analyte.
      if(any(is.na(chromIndices))){
        warning("Chromatogram indices for ", analyte, " are missing in ", fileInfo[run, "runName"])
        message("Skipping ", analyte, " in ", fileInfo[run, "runName"], ".")
        next
      } else {
        if(params[["chromFile"]] =="mzML") fetchXIC = extractXIC_group
        if(params[["chromFile"]] =="sqMass") fetchXIC = extractXIC_group2
        XICs <- fetchXIC(mzPntrs[[run]], chromIndices = chromIndices)
      }
      if(params[["smoothPeakArea"]]){
        XICs <- smoothXICs(XICs, type = params[["XICfilter"]], kernelLen = params[["kernelLen"]], polyOrd = params[["polyOrd"]])
      }
      area <- calculateIntensity(XICs, df[1, "leftWidth"], df[1, "rightWidth"],  integrationType = params[["integrationType"]],
                                 baselineType = params[["baselineType"]], fitEMG = params[["fitEMG"]], baseSubtraction = params[["baseSubtraction"]])
      newArea[[run]][i] <- area
    }
  }

  for(mz in mzPntrs){
    if(is(mz)[1] == "SQLiteConnection") DBI::dbDisconnect(mz)
    if(is(mz)[1] == "mzRpwiz") rm(mz)
  }

  newArea <- as.data.frame(do.call(cbind, newArea))
  newArea$precursor <- analytes
  newArea <- tidyr::pivot_longer(newArea, -.data$precursor, names_to = "run",
                                 values_to = "intensity") %>% as.data.frame()
  newArea$run <- fileInfo[newArea$run, "runName"]
  newArea
}



reIntensity <- function(df, XICs, params){
  idx <- which(df[["alignment_rank"]] == 1)
  for(i in idx){
    analyte_chr <- as.character(df$transition_group_id[i])
    area <- calculateIntensity(XICs[[analyte_chr]], df$leftWidth[i], df$rightWidth[i],
                               params[["integrationType"]], params[["baselineType"]], params[["fitEMG"]])
    df$intensity[i] <- area
  }
  df
}
