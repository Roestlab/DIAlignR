#' Converts to OpenMS format
#'
#' MSPLIT returns an OpenSWATH assay library which is not compatible for OpenMS version 3.0
#' This function reformat the library for TargetedFileConverter (OpenMS).
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-06-14
#' @import dplyr
#' @importFrom magrittr %>%
#' @param filename (string) Path to the MSPLIT output library. The file extension must be .txt.
#' @return (None) writes a .tsv file.
#' @keywords internal
msplit_to_openms <- function(filename){
  lib <- read.table(filename, header = TRUE, sep = "\t", numerals = "no.loss",
                    comment.char = "", stringsAsFactors = FALSE)
  # Remove three columns
  lib2 <- dplyr::select(lib, -c(MissedCleavages, Replicates, NrModifications))
  # Rename columns
  lib2 <- dplyr::rename(lib2, NormalizedRetentionTime = Tr_recalibrated,
                        CollisionEnergy = CE, Decoy = decoy, ProteinId = ProteinName,
                        ModifiedPeptideSequence = FullUniModPeptideName,
                        LabelType = GroupLabel)
  lib2$TransitionGroupId <- substr(x = lib2$transition_group_id, start = 1,
                                   stop = regexpr(pattern = "_", lib2$transition_group_id)-1)
  lib2$TransitionGroupId <- as.integer(lib2$TransitionGroupId)
  lib2$TransitionGroupId <- lib2$TransitionGroupId - min(lib2$TransitionGroupId)

  start <- sapply(gregexpr(pattern = "/", lib2$transition_name), `[`, 2) + 1
  stp <- sapply(lib2$transition_name, nchar)
  lib2$ProductCharge <- substr(x = lib2$transition_name, start = start, stop = stp)

  lib2$PeptideGroupLabel <- NA_character_
  lib2$CompoundName <- NA_character_
  lib2$SumFormula <- NA_character_
  lib2$SMILES <- NA_character_
  lib2$Adducts <- NA_character_
  lib2$UniprotId <- NA_character_
  lib2$GeneName <- "NA"
  lib2$FragmentType <- substr(lib2$Annotation, 1, 1)
  lib2$FragmentSeriesNumber <- substr(lib2$Annotation, 2, 3)
  lib2$PrecursorIonMobility <- -1L
  lib2$TransitionId <- seq(0, nrow(lib2)-1)
  lib2$DetectingTransition <- 1L
  lib2$IdentifyingTransition <- 0L
  lib2$QuantifyingTransition <- 1L
  lib2$Peptidoforms <- NA_character_

  lib2 <- lib2 %>% select(PrecursorMz, ProductMz, PrecursorCharge, ProductCharge,	LibraryIntensity,
                          NormalizedRetentionTime, PeptideSequence, ModifiedPeptideSequence, PeptideGroupLabel,
                          LabelType, CompoundName, SumFormula, SMILES, Adducts,	ProteinId, UniprotId,
                          GeneName,	FragmentType,	FragmentSeriesNumber,	Annotation,	CollisionEnergy,
                          PrecursorIonMobility, TransitionGroupId, TransitionId, Decoy, DetectingTransition,
                          IdentifyingTransition, QuantifyingTransition,	Peptidoforms)

  filename <- paste0(substr(filename, 1, nchar(filename)-4), ".tsv")
  write.table(lib2, file = filename, sep = "\t", row.names = FALSE,
              quote = FALSE, na = "")
}


# Library alignment function
libraryAlignment <- function(assayLibs, x, y, alignType = "loess", spanvalue = 0.1){
  df.x <-  assayLibs[[x]] %>%  dplyr::select(.data$tempID, .data$NormalizedRetentionTime)
  df.y <-  assayLibs[[y]] %>% dplyr::select(.data$tempID, .data$NormalizedRetentionTime)
  RUNS_RT <- dplyr::inner_join(df.x, df.y, by = "tempID", suffix = c(".x", ".y"))
  # TODO: Update spanvalue or alignment method if number of points are less.
  if(alignType == "linear"){
    fit <- lm(NormalizedRetentionTime.y ~ NormalizedRetentionTime.x, data = RUNS_RT)
  } else {
    fit <- loess(NormalizedRetentionTime.y ~ NormalizedRetentionTime.x,
                 data = RUNS_RT, span = spanvalue, control=loess.control(surface="direct"))
  }
  # direct surface allows to extrapolate outside of training data boundary while using predict.
  fit
}


#' Get master library from MSPLIT
#'
#' MSPLIT returns an OpenSWATH assay libraries which can be combined to obtain a comprehensive library
#' This library allows to perform false-discovery rate extimation on all samples together and further use
#' alignment algorithm of DIAlignR.
#'
#' In star alignment, a base library is picked among n runs. Other n-1 runs are aligned to it either by linear
#' or loess alignment. Using qvalue IDs are picked from a run. These IDs include mz, library intensities and
#' aligned retention time. IDs are added to the base library that is , finally, written as masterLib.tsv in the working
#' directory.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-06-14
#' @import dplyr
#' @importFrom magrittr %>%
#' @param msplitFilteredPattern (string) pattern of MSPLITfiltered.txt files. Pattern must start with '*'.
#' @param msplitAssayPattern (string) pattern of MSPLITassay_library.tsv files. Pattern must start with '*'.
#' @param spanvalue (numeric) spanvalue for LOESS fit. For targeted proteomics 0.1 could be used.
#' @return (None) writes a masterLib.tsv file.
#' @examples
#' msplitFilteredPattern <- "*_MSPLITfiltered.txt"
#' msplitAssayPattern <- "*_Openswath_assaylib.tsv"
#' \dontrun{
#' msplit_masterLib(msplitFilteredPattern, msplitAssayPattern)
#' }
#' @keywords internal
msplit_masterLib <- function(msplitFilteredPattern, msplitAssayPattern, alignType = "loess", spanvalue = 0.1){
  # Read MSPLITfiltered.txt file. This contains the spectra that were used to build Openswath_assaylib.
  temp <- list.files(pattern = msplitFilteredPattern)
  message("Found ", length(temp), " files with ", msplitFilteredPattern, " pattern.")
  MSPLITfiltered <- list()
  for(i in seq_along(temp)){
    MSPLITfiltered[[i]] <- read.table(temp[i], header = TRUE, sep = "\t", comment.char = "", stringsAsFactors = FALSE)
    # Get Peptide, Mz.1, z.1, experiment RT, library RT, p-value, q-value
    # TODO: Add Unimod in Peptide
    MSPLITfiltered[[i]] <- MSPLITfiltered[[i]][, c(5,6,7,17,18, 34,35)]
    MSPLITfiltered[[i]][, "stat18"] <- as.numeric(sapply(MSPLITfiltered[[i]][, "stat18"], function(s) substr(s, 1, nchar(s)-3)))
    #Loss of some precision because there are 19 digits after decimal in character format.
  }
  temp <-substr(temp, 1, nchar(temp)-nchar(msplitFilteredPattern)+1)
  names(MSPLITfiltered) <- temp

  ###### Read Openswath_assaylib.tsv  #########
  temp <- list.files(pattern = msplitAssayPattern)
  message("Found ", length(temp), " files with ", msplitAssayPattern, " pattern.")
  OpenswathAssaylib <- list()
  for(i in seq_along(temp)){
    assayLib <- read.table(temp[i], header = TRUE, sep = "\t", comment.char = "", stringsAsFactors = FALSE)
    assayLib <- assayLib[assayLib$Decoy == 0, ]
    assayLib <- assayLib[, c("PrecursorMz", "NormalizedRetentionTime", "PeptideSequence", "ModifiedPeptideSequence", "PrecursorCharge")]
    assayLib <- dplyr::distinct(assayLib)
    OpenswathAssaylib[[i]] <- assayLib
  }
  temp <-substr(temp, 1, nchar(temp)-nchar(msplitAssayPattern)+1)
  names(OpenswathAssaylib) <- temp

  ###### Merge MSPLIT outputs. We get q-value for each entry in Openswath_assaylib.tsv  ######
  runs <- intersect(names(OpenswathAssaylib), names(MSPLITfiltered))
  names(runs) <- paste("run", 1:length(OpenswathAssaylib), sep = "")
  # TODO: right_join on ModifiedPeptideSequence
  for (run in runs){
    OpenswathAssaylib[[run]] <- MSPLITfiltered[[run]] %>%
      dplyr::select(Peptide, stat1, z.1, stat19) %>%
      dplyr::right_join(OpenswathAssaylib[[run]], by = c("stat1" = "NormalizedRetentionTime", "Peptide" =  "PeptideSequence", "z.1" = "PrecursorCharge")) %>%
      dplyr::rename(NormalizedRetentionTime = stat1, PeptideSequence = Peptide, PrecursorCharge = z.1, qvalue = stat19) %>%
      dplyr::distinct()
  }

  ###### Get peptide sequence and charge combination that acts as temporary ID for each precursor  ######
  AllPrecursors <- c()
  for(run in runs){
    newIDs <-OpenswathAssaylib[[run]] %>% select(ModifiedPeptideSequence, PrecursorCharge) %>%
      tidyr::unite("tempID", sep="/") %>% .$tempID
    if(any(duplicated(newIDs))){
      warning("Duplication of transition_group_id in run ", run)
    }
    OpenswathAssaylib[[run]]$tempID <- newIDs
    AllPrecursors <- union(AllPrecursors, OpenswathAssaylib[[run]]$tempID)
  }
  message("There are total ", length(AllPrecursors), " precursors.")


  ###### Find precursors that are missing in one or more libraries.  ######
  commonPeps <- AllPrecursors
  for(run in runs){
    commonPeps <- intersect(OpenswathAssaylib[[run]]$tempID, commonPeps)
  }
  message(length(commonPeps), " are in all libraries.")
  rm(MSPLITfiltered)


  # Find the base library based on number of precursors it has
  num_of_pep <- sapply(OpenswathAssaylib, function(x) nrow(x))
  refRunIdx <- unname(which.max(num_of_pep))

  # Create base library from the library with most entries.
  # Append empty rows in the base library which it is missing currently.
  starLib <- OpenswathAssaylib[[refRunIdx]]
  starLib$refRunIdx <- NA_integer_
  starLib <- starLib %>% select(PrecursorMz, NormalizedRetentionTime, PeptideSequence,
                                ModifiedPeptideSequence, PrecursorCharge, tempID,
                                qvalue, refRunIdx)
  peps <- setdiff(AllPrecursors, starLib$tempID)
  peps <- data.frame("PrecursorMz" = rep(NA_real_, length(peps)),
                     "NormalizedRetentionTime" = rep(NA_real_, length(peps)),
                     "tempID" = peps,
                     "PeptideSequence" = rep(NA_character_, length(peps)),
                     "ModifiedPeptideSequence" = rep(NA_character_, length(peps)),
                     "PrecursorCharge" = rep(NA_integer_, length(peps)),
                     "qvalue" = rep(NA_real_, length(peps)),
                     "refRunIdx" = rep(NA_integer_, length(peps)),
                     stringsAsFactors = FALSE)
  starLib <- rbind(starLib, peps)

  # Get a table of qvalues from all runs
  qvalTbl <- lapply(OpenswathAssaylib , function(df) df[,c("tempID", "qvalue")])
  qvalTbl <- Reduce(function(d1, d2) merge(d1, d2, by = "tempID", all = TRUE, no.dups = TRUE), qvalTbl)

  # Get a table of retention times from all runs
  rtTbl <- lapply(OpenswathAssaylib , function(df) df[,c("tempID", "NormalizedRetentionTime")])
  rtTbl <- Reduce(function(d1, d2) merge(d1, d2, by = "tempID", all = TRUE, no.dups = TRUE), rtTbl)

  if(!all(names(OpenswathAssaylib) == runs)){
    stop()
  }
  names(qvalTbl) <- c("tempID", names(runs))
  names(rtTbl) <- c("tempID", names(runs))

  # Set order of rows in qvalue table and RT table
  qvalTbl <- qvalTbl[match(starLib$tempID,qvalTbl$tempID),]
  rownames(qvalTbl) <- qvalTbl[,"tempID"]
  qvalTbl <- qvalTbl[,-c(1)]

  rtTbl <- rtTbl[match(starLib$tempID,rtTbl$tempID),]
  rownames(rtTbl) <- rtTbl[,"tempID"]
  rtTbl <- rtTbl[,-c(1)]

  # Star alignment
  ######## Get pairwise alignments of all libraries to the base library. #######
  fits <- list()
  for(i in 1:length(runs)){
    run <- names(runs)[i]
    fits[[paste0("ref_", run)]] <- libraryAlignment(OpenswathAssaylib, x = runs[[i]],
                                                         y = runs[[refRunIdx]], alignType, spanvalue)
  }

  ###### Get base library retention time.  ######
  # Get the reference library for each analyte based on qvalue
  # Map its retention time to the base library.
  for (i in 1:nrow(starLib)){
    run <- names(which.max(qvalTbl[i,]))
    starLib$refRunIdx[i] <- unname(which.max(qvalTbl[i,]))
    # map retention time to the base library.
    starLib$NormalizedRetentionTime[i] <- predict(fits[[paste0("ref_", run)]], newdata = data.frame(NormalizedRetentionTime.x = rtTbl[i, run]))
    starLib$qvalue[i] <- qvalTbl[i,run]
  }
  if(any(is.na(starLib$NormalizedRetentionTime))){
    warning(which(is.na(starLib$NormalizedRetentionTime)) ,"have missing RT caliberation.")
  }

  # Create new transition_group_id
  starLib$TransitionGroupId <- seq(0, nrow(starLib)-1)

  ###### Add dignostic plots and table ######
  masterAlignSummary(fits, runs)

  ####### Add tempID in all Openswath assay libraries ########
  # Read assay library and update transition_group_id and transition_name.
  temp <- list.files(pattern = msplitAssayPattern)
  OpenswathAssaylib <- list()
  for(i in 1:length(temp)){
    OpenswathAssaylib[[i]] <- read.table(temp[i], header = TRUE, sep = "\t", comment.char = "", stringsAsFactors = FALSE, numerals = "no.loss")
    OpenswathAssaylib[[i]] <- OpenswathAssaylib[[i]][OpenswathAssaylib[[i]]$Decoy == 0, ]
  }
  temp <-substr(temp, 1, nchar(temp)-nchar(msplitAssayPattern)+1)
  names(OpenswathAssaylib) <- temp

  if(!all(names(OpenswathAssaylib) == runs)){
    stop("Names of OpenswathAssaylib are not matching with runs.")
  }

  # Add tempIDs
  for(run in runs){
    # Check if the transition_group_ids are duplicated in the run.
    newIDs <- OpenswathAssaylib[[run]] %>%
      select(PrecursorMz, NormalizedRetentionTime, PeptideSequence,
             ModifiedPeptideSequence, PrecursorCharge) %>% distinct() %>%
      select(ModifiedPeptideSequence, PrecursorCharge) %>%
      tidyr::unite("tempID", sep="/") %>% .$tempID
    if(any(duplicated(newIDs))){
      stop("Duplication of transition_group_id in run ", run)
    }
    newIDs <- OpenswathAssaylib[[run]] %>%
      select(ModifiedPeptideSequence, PrecursorCharge) %>%
      tidyr::unite("tempID", sep="/") %>% .$tempID
    # Add tempID in the assay library
    OpenswathAssaylib[[run]]$tempID <- newIDs
  }

  ######## Start filling the base library ###############
  # Get the master library template
  masterLib <- data.frame("PrecursorMz"=double(), "ProductMz"=double(), "PrecursorCharge"=integer(), "ProductCharge"=integer(),
                          "LibraryIntensity"=double(), "NormalizedRetentionTime" = double(), "PeptideSequence"=character(),  "ModifiedPeptideSequence"=character(),
                          "PeptideGroupLabel"=character(), "LabelType"=character(), "CompoundName"=character(), "SumFormula"=character(),
                          "SMILES"=character(), "Adducts"=character(), "ProteinId"=character(), "UniprotId"=character(), "GeneName"=character(), "FragmentType"=character(),
                          "FragmentSeriesNumber"=integer(), "Annotation"=character(), "CollisionEnergy"=double(), "PrecursorIonMobility" = integer(),
                          "TransitionGroupId"=integer(), TransitionId =integer(), "Decoy"=integer(), "DetectingTransition"=integer(),
                          "IdentifyingTransition"=integer(), "QuantifyingTransition"=integer(), "Peptidoforms"=character(),
                          "tempID"=character(), stringsAsFactors=FALSE)

  # For all analytes get the row from assay library and add to master library.
  # Thus keep all transition mz and intensities from best qvalue.
  # Add the retention time and trantion group ID from the base library.
  for(i in 1:nrow(starLib)){
    run <- runs[[starLib$refRunIdx[i]]]
    df <- OpenswathAssaylib[[run]][OpenswathAssaylib[[run]]$tempID == starLib$tempID[i], ]
    df$TransitionGroupId <- starLib$TransitionGroupId[i]
    df$NormalizedRetentionTime <- starLib$NormalizedRetentionTime[i]
    masterLib <- rbind(masterLib, df)
  }
  masterLib$TransitionId <- seq(0, nrow(masterLib) -1)
  masterLib <- dplyr::select(masterLib, -.data$tempID)
  # Now master library has same format as assay library. It has RT and Transition Group ID from the base library.
  write.table(masterLib, file = "masterLib.tsv", sep = "\t", row.names = FALSE, quote = FALSE, na = "")

  message("Number of precursors from each run")
  k <- table(starLib$refRunIdx)
  names(k) <- runs
  print(k)
}


masterAlignSummary <- function(fits, runs){
  df <- data.frame(runs)
  df$rse <- round(sapply(fits, getRSE), 2)
  df$n <- sapply(fits, function(fit) ifelse(is(fit, "loess"), fit$n, nobs(fit)))
  print(df)
  pdf("run_aligned_master_library.pdf")
  for(i in seq_along(fits)){
    fit <- fits[[i]]
    # Plot the fit. Separate code for loess and linear fit.
    if(is(fit, "loess")){
      plot(fit$x, fit$y, pch=19, cex=0.1, xlab = runs[i], ylab = "master library",
           main = paste0("RSE = ", df$rse[i], ", N = ", df$n[i]))
      j <- order(fit$x)
      lines(fit$x[j], fit$fitted[j],col="green",lwd=2)
    } else{
      plot(fit$model$x, fit$model$y, pch=19, cex=0.1, xlab = runs[i], ylab = "master library",
           main = paste0("RSE = ", df$rse[i], ", N = ", df$n[i]))
      j <- order(fit$model$x)
      lines(fit$model$x[j], fit$fitted.values[j],col="green",lwd=2)
    }
  }
  dev.off()
}


#' Get psuedo iRT spectra
#'
#' This function cuts the retention time axis in multiple zones from which intense and confident spectra
#' are selected for the purspose of buillding psuedo iRT library.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-06-14
#' @import dplyr
#' @importFrom magrittr %>%
#' @param msplitFilteredPattern (string) pattern of MSPLITfiltered.txt files. Pattern must start with '*'.
#' @param maxQval (numeric) upper limit of qvalue for spectra to be considered for psuedo iRT.
#' @param zoneWidth (numeric)
#' @param minPep (integer) minimum number of precursors requires in each retention time zone.
#' @param maxPep (integer) maximum number of precursors requires in each retention time zone.
#' @param minIntensity (numeric) lower limit of ion count for spectra to be considered for psuedo iRT.
#' @return (None) writes *_iRT.txt files.
#' @examples
#' msplitFilteredPattern <- "*_MSPLITfiltered.txt"
#' \dontrun{
#' psuedo_iRTspectra(msplitFilteredPattern)
#' }
#' @keywords internal
psuedo_iRTspectra <- function(msplitFilteredPattern, maxQval = 1e-04, zoneWidth =300,
                     minPep = 5L, maxPep = 10L, minIntensity = 1e06){
  temp <- list.files(pattern = "*MSPLITfiltered.txt")
  message("Found ", length(temp), " files with ", msplitFilteredPattern, " pattern.")

  for (filename in temp){
    message("Reading ", filename)
    MSPLITfiltered <- read.table(filename, header = TRUE, sep = "\t", numerals = "no.loss",
                                 comment.char = "", stringsAsFactors = FALSE)
    # Get Peptide, Mz.1, z.1, ion count, experiment RT, library RT, q-value
    MSPLITfiltered <- MSPLITfiltered[, c(5,6,7,15, 17,18, 35)]
    class(MSPLITfiltered$stat19) <- "numeric"
    if(class(MSPLITfiltered$IonCount) != "numeric"){
      class(MSPLITfiltered$IonCount)
    }
    # Select spectra that have qvalue less than maxQval. Avoid Methionine oxidation spectra.
    # For each precusor get the spectra with maximum ion count.
    # Split retention time zones in 24 zones of 300 seconds length.
    # In each zone, select top 10 intense precursors.
    # Select (minimum 5) precursors that have ion count greater than minIntensity,
    # if not found, lower intense precursors are included.
    # Output precursor's sequence, m/z, charge and retention time.
    MSPLITfiltered <- MSPLITfiltered %>% filter(stat19 < maxQval & !grepl("M+15.99", Peptide, fixed = TRUE)) %>%
      group_by(Peptide, z.1) %>% filter(IonCount == max(IonCount))  %>% ungroup() %>%
      arrange(stat1) %>% mutate(RTzone =  stat1 %/% zoneWidth) %>%
      group_by(RTzone) %>%
      arrange(RTzone, desc(IonCount), .by_group = TRUE) %>% top_n(maxPep, IonCount) %>%
      top_n(max(maxPep - sum(IonCount <= minIntensity), minPep), IonCount) %>% ungroup() %>%
      select(Peptide, Mz = Mz.1, z = z.1, stat1) %>% as.data.frame()
    filename <- sub('\\.txt$', '', filename)
    filename <- paste0(filename, "_iRT.txt")
    write.table(MSPLITfiltered, file = filename, sep = "\t", row.names = FALSE, quote = FALSE)
    message("Filetered ", nrow(MSPLITfiltered), " psuedo-iRT spectra in ", filename)
  }
}


#' Builds assay library and pseudo iRT library
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-06-14
#' @import dplyr
#' @importFrom magrittr %>%
#' @param msplitFilteredPattern (string) pattern of pseudo iRT spectra files. Pattern must start with '*'.
#' @param msplitAssaySuffix (numeric) suffix of MSPLIT assay library.
#' @param masterLibPath (numeric) path of the master library to be converted to sample-specific assay library.
#' @param startTime (integer) start time of the run in second. Usually it is 0.
#' @param stopTime (integer) stop time of the run in seconds. Usually it is 7200 sec (2 hour run).
#' @return (None) writes _norm_iRT_assaylib.tsv and _norm_masterLib.tsv files.
#' @examples
#' psuedoSpectraPattern <- "*_MSPLITfiltered_iRT.txt"
#' msplitAssaySuffix <- "_Openswath_assaylib.tsv"
#' masterLibPath <- "masterLib.tsv"
#' \dontrun{
#' psuedo_iRTlibrary(psuedoSpectraPattern, msplitAssaySuffix, masterLibPath)
#' }
#' @keywords internal
psuedo_iRTlibrary <- function(psuedoSpectraPattern, msplitAssaySuffix, masterLibPath,
                         startTime = 0, stopTime = 7200){
  # Read master library
  masterLib <- read.table(masterLibPath, header = TRUE, sep = "\t", numerals = "no.loss",
                          comment.char = "", stringsAsFactors = FALSE)
  masterDF <- masterLib %>% dplyr::select(PeptideSequence, ModifiedPeptideSequence,
                                          PrecursorCharge, NormalizedRetentionTime)

  temp <- list.files(pattern = psuedoSpectraPattern)
  message("Found ", length(temp), " files containing pseudo iRT spectra.")
  for (filename in temp){
    # Read pseudo iRT spectra files
    message("Reading ", filename)
    iRT <- read.table(filename, header = TRUE, sep = "\t", numerals = "no.loss",
                      comment.char = "", stringsAsFactors = FALSE)
    # TODO: Add Unimod in Peptide

    # Get corresponding MSPLIT OpenSWATH assay library
    run <- substr(filename, 1, nchar(filename)-nchar(psuedoSpectraPattern)+1)
    assayFile <- paste0(run, msplitAssaySuffix)
    message("Reading ", assayFile)
    assayLib <- read.table(assayFile, header = TRUE, sep = "\t", numerals = "no.loss",
                           comment.char = "", stringsAsFactors = FALSE)

    # Filter assay library with pseudo iRT spectra file.
    iRT_assay <- iRT %>% dplyr::select(Peptide, z) %>%
      dplyr::inner_join(assayLib, by = c("Peptide"="ModifiedPeptideSequence", "z"="PrecursorCharge")) %>%
      dplyr::rename(ModifiedPeptideSequence = Peptide, PrecursorCharge = z) %>%
      dplyr::select(PrecursorMz, ProductMz, PrecursorCharge, ProductCharge, LibraryIntensity,
             NormalizedRetentionTime, PeptideSequence, ModifiedPeptideSequence, PeptideGroupLabel,
             LabelType,	CompoundName,	SumFormula,	SMILES,	Adducts, ProteinId,
             UniprotId, GeneName, FragmentType, FragmentSeriesNumber, Annotation,
             CollisionEnergy, PrecursorIonMobility, TransitionGroupId, TransitionId,
             Decoy, DetectingTransition, IdentifyingTransition, QuantifyingTransition,
             Peptidoforms)

    # Normalize each iRT. Generate respective master library

    # Normalize assay library iRT time to 0-100 scale.
    iRT_assay[["NormalizedRetentionTime"]] = round(100*(iRT_assay[["NormalizedRetentionTime"]] - startTime)/stopTime, 3)

    # Select common rows between master library and iRT assay library
    iRT_assayDF <- iRT_assay %>% dplyr::select(PeptideSequence, ModifiedPeptideSequence, PrecursorCharge, NormalizedRetentionTime)
    DF <- dplyr::inner_join(masterDF, iRT_assayDF,
                            by = c("PeptideSequence", "ModifiedPeptideSequence", "PrecursorCharge")) %>% dplyr::distinct()

    # Get a linear transformation (assayRT ~ masterRT)
    fit <- lm(NormalizedRetentionTime.y ~ NormalizedRetentionTime.x, data = DF)
    message(run)
    message("rsq between experimental retention time and indexed retention time is ", summary(fit)$adj.r.squared)

    # Map master library to assay library linearly. Create a copy of master library.
    newVal <- predict(fit, newdata = data.frame(NormalizedRetentionTime.x = masterDF[["NormalizedRetentionTime"]]))
    masterLib2 <- masterLib
    masterLib2[["NormalizedRetentionTime"]] <- round(newVal, 3)

    # Write normalized iRT assay library.
    iRT_assayFile <- paste0(run, "_norm_iRT_assaylib.tsv")
    write.table(iRT_assay, file = iRT_assayFile, sep = "\t", row.names = FALSE,
                quote = FALSE, na = "")
    message(iRT_assayFile, " pseudo iRT library is written.")

    # Write normalized master library.
    masterFile <- paste0(run, "_norm_masterLib.tsv")
    write.table(masterLib2, file = masterFile, sep = "\t", row.names = FALSE, quote = FALSE)
    message(masterFile, " assay library is written.")
  }
}


#' Builds assay library and pseudo iRT library
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-06-25
#' @import dplyr
#' @importFrom magrittr %>%
#' @param masterLibPath (numeric) path of the master library to be converted to sample-specific assay library.
#' @param iRTLibPath (string) pattern of pseudo iRT spectra files. Pattern must start with '*'.
#' @param msplitAssaySuffix (numeric) suffix of MSPLIT assay library.
#' @param lib_lib_align (integer) start time of the run in second. Usually it is 0.
#' @param spanvalue (integer) stop time of the run in seconds. Usually it is 7200 sec (2 hour run).
#' @return (None) writes _norm_iRT_assaylib.tsv and _norm_masterLib.tsv files.
#' @examples
#' masterLibPath <- "masterLib.tsv"
#' iRTLibPath <- "iRTassays.tsv"
#' msplitAssaySuffix <- "_Openswath_assaylib.tsv"
#' \dontrun{
#' iRTAlignment(masterLibPath, iRTLibPath, msplitAssaySuffix, lib_lib_align = "loess")
#' }
#' @keywords internal
iRTAlignment <- function(masterLibPath, iRTLibPath, msplitAssaySuffix, lib_lib_align = "linear", spanvalue = 0.1){
  # Read master library
  masterLib <- read.table(masterLibPath, header = TRUE, sep = "\t", numerals = "no.loss",
                          comment.char = "", stringsAsFactors = FALSE)
  masterDF <- masterLib %>% dplyr::select(ModifiedPeptideSequence, PrecursorCharge,
                                          NormalizedRetentionTime) %>% dplyr::distinct()

  # Read iRT library
  irtLib <- read.table(iRTLibPath, header = TRUE, sep = "\t", numerals = "no.loss",
                          comment.char = "", stringsAsFactors = FALSE)
  irtDF <- irtLib %>% dplyr::select(ModifiedPeptideSequence, PrecursorCharge,
                                       NormalizedRetentionTime)  %>% dplyr::distinct()

  # Align assay library to iRT (linearly) to normalize retention time
  masterNormDF <- masterDF
  temp <- dplyr::inner_join(masterNormDF, irtDF, by = c("ModifiedPeptideSequence", "PrecursorCharge"))
  fit <- lm(NormalizedRetentionTime.y ~ NormalizedRetentionTime.x, data = temp)
  message("RT of the master library is normalized with ", nobs(fit), " iRT peptides.")
  message("Linear fit has adjusted R-squared of ", summary(fit)$adj.r.squared)
  newVal <- predict(fit, newdata = data.frame(NormalizedRetentionTime.x = masterNormDF[["NormalizedRetentionTime"]]))
  masterNormDF[["NormalizedRetentionTime"]] <- round(newVal, 3)
  masterNormDF <-  dplyr::left_join(masterLib, masterNormDF,
                               by = c("ModifiedPeptideSequence", "PrecursorCharge"))
  masterNormDF <- masterNormDF %>% dplyr::select(-NormalizedRetentionTime.x) %>%
    dplyr::rename(NormalizedRetentionTime = NormalizedRetentionTime.y) %>%
    dplyr::relocate(NormalizedRetentionTime, .before = PeptideSequence)

  write.table(masterNormDF, file = "masterLib_norm.tsv", sep = "\t", row.names = FALSE, quote = FALSE, na = "")

  temp <- list.files(pattern = msplitAssaySuffix)
  message("Found ", length(temp), " assay libraries.")
  for (filename in temp){
    # Get MSPLIT OpenSWATH assay library
    run <- substr(filename, 1, nchar(filename)-nchar(msplitAssaySuffix))
    message("Reading ", filename)
    assayLib <- read.table(filename, header = TRUE, sep = "\t", numerals = "no.loss",
                           comment.char = "", stringsAsFactors = FALSE)

    # Join assay library and master library
    assayDF <- assayLib %>% dplyr::select(ModifiedPeptideSequence, PrecursorCharge,
                                          NormalizedRetentionTime) %>% dplyr::distinct()
    assayDF <- dplyr::left_join(masterDF, assayDF,
                            by = c("ModifiedPeptideSequence", "PrecursorCharge"))

    # We need to convert master libray to assay libray. Everything remain same but the retention time.
    # Precursors already present in assay library, have retention time copied over.
    # Other precursors will have (linear/non-lineasrly) modified retention time.
    if(lib_lib_align == "linear"){
      fit <- lm(NormalizedRetentionTime.y ~ NormalizedRetentionTime.x, data = assayDF, na.action = na.omit)
    } else {
      fit <- loess(NormalizedRetentionTime.y ~ NormalizedRetentionTime.x, data = assayDF,
                         span = spanvalue, na.action = na.omit, control=loess.control(surface="direct"))
    }
    idx <- which(is.na(assayDF[["NormalizedRetentionTime.y"]]))
    newVal <- predict(fit, newdata = data.frame(NormalizedRetentionTime.x = assayDF[idx,"NormalizedRetentionTime.x"]))
    assayDF[["NormalizedRetentionTime.y"]][idx] <- round(newVal, 3)
    assayDF <- assayDF %>% select(ModifiedPeptideSequence, PrecursorCharge, NormalizedRetentionTime.y) %>%
      rename(NormalizedRetentionTime = NormalizedRetentionTime.y)
    message("Added ", length(idx), " precursors in ", run, " library.")
    message("Retention times are mapped with a ", lib_lib_align, " method using ", nrow(assayDF) - length(idx), " points.")
    message("Residual Standard error = ", round(getRSE(fit), 3))

    # Align assay library to iRT to normalize retention time
    temp <- dplyr::inner_join(assayDF, irtDF, by = c("ModifiedPeptideSequence", "PrecursorCharge"))
    # Get a linear transformation of assay library to iRT
    fit <- lm(NormalizedRetentionTime.y ~ NormalizedRetentionTime.x, data = temp)
    message("RT of the library is normalized with ", nobs(fit), " iRT peptides.")
    message("Linear fit has adjusted R-squared of ", summary(fit)$adj.r.squared)
    newVal <- predict(fit, newdata = data.frame(NormalizedRetentionTime.x = assayDF[["NormalizedRetentionTime"]]))
    assayDF[["NormalizedRetentionTime"]] <- round(newVal, 3)

    # Copy other columns from master to assay library
    assayDF <-  dplyr::left_join(masterLib, assayDF,
                                    by = c("ModifiedPeptideSequence", "PrecursorCharge"))
    assayDF <- assayDF %>% dplyr::select(-NormalizedRetentionTime.x) %>%
      dplyr::rename(NormalizedRetentionTime = NormalizedRetentionTime.y) %>%
      dplyr::relocate(NormalizedRetentionTime, .before = PeptideSequence)

    filename <- paste0(run,  "_master_norm.tsv")
    write.table(assayDF, file = filename, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  }
}
