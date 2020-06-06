# MCP paper 2019
# Gupta et al., 2019, Molecular & Cellular Proteomics 18, 806 –817
# https://doi.org/10.1074/mcp.TIR118.001132

# Error curve
plotErrorCurve <- function(x, clr, SameGraph, xmax = 120, ...){
  x <- x[!is.na(x)]
  breaks = seq(0, xmax, by=0.5)
  duration.cut = cut(x, breaks, right = FALSE)
  duration.freq = table(duration.cut)
  cumfreq0 = c(0, cumsum(duration.freq))
  if(SameGraph == TRUE){lines(breaks, cumfreq0/length(x), col = clr, ...)}
  else{plot(breaks, cumfreq0/length(x), col = clr, type = "l", ...)}
}

# Error statistics
statPerPair <- function(AlignError, peptides, pw, deltaTime){
  num_of_pep <- sum(!is.na(AlignError))
  MatchedWithinHalfPeak <- sum(abs(AlignError) <= (0.5*pw*deltaTime), na.rm = TRUE)
  MatchedWithin1PW <- sum(abs(AlignError) <= (1*pw*deltaTime), na.rm = TRUE)
  MatchedWithin2PW <- sum(abs(AlignError) <= (2*pw*deltaTime), na.rm = TRUE)
  NApep <- sum(is.na(AlignError))
  Peptides1PeakDistAbv <- paste(peptides[which(abs(AlignError) > (pw*deltaTime))], collapse = " ")
  Peptides2PeakDistAbv <- paste(peptides[which(abs(AlignError) > (2*pw*deltaTime))], collapse = " ")
  ErrorStatPair <- c(MatchedWithinHalfPeak, MatchedWithin1PW, MatchedWithin2PW,
                     MatchedWithinHalfPeak*100/num_of_pep, MatchedWithin1PW*100/num_of_pep,
                     MatchedWithin2PW*100/num_of_pep, NApep)
  return(ErrorStatPair)
}

# Get precursors
mcp_precursors <- function(fileInfo){
  oswName <- unique(fileInfo[["featureFile"]])
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = as.character(oswName))
  query <- "SELECT DISTINCT PRECURSOR.ID AS transition_group_id,
      TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID AS transition_id,
      PEPTIDE.ID AS peptide_id,
      PEPTIDE.MODIFIED_SEQUENCE AS sequence,
      PRECURSOR.CHARGE AS charge,
      PRECURSOR.GROUP_LABEL AS group_label
      FROM PRECURSOR
      INNER JOIN TRANSITION_PRECURSOR_MAPPING ON TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID = PRECURSOR.ID
      INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID = PRECURSOR.ID
      INNER JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID
      ORDER BY transition_group_id, transition_id;"
  # Run query to get peptides, their coordinates and scores.
  precursorsInfo <- tryCatch(expr = { output <- DBI::dbSendQuery(con, statement = query)
    DBI::dbFetch(output)},
    finally = {DBI::dbClearResult(output)
    DBI::dbDisconnect(con)})
  # Each precursor has only one row. tidyr::nest creates a tibble object that is twice as heavy to regular list.
  library(magrittr)
  precursorsInfo <- dplyr::group_by(precursorsInfo, .data$transition_group_id, .data$peptide_id, .data$sequence, .data$charge, .data$group_label) %>%
    dplyr::summarise(transition_ids = base::list(.data$transition_id)) %>% dplyr::ungroup() %>% as.data.frame()
  precursorsInfo
}

# Figure 3
figure3 <- function(linear_err, loess_err_075, loess_err_003, Hybrid_err,
                    analytes, strep0Pairs, strep10Pairs, strep0_10Pairs){
  pdf("Figure3.pdf")
  xmax <- 80; pw <- 9; deltaTime<-3.4
  par(mgp=c(2.5,1,0))

  # Figure 3A
  plotErrorCurve(abs(c(as.matrix(linear_err))), "green", FALSE, xmax, main = "Between all 120 pairs", ylim = c(0, 1.02), xlab = "Retention time difference (in sec)", ylab = "Cumulative fraction of peptides", lwd = 2, cex.lab = 1.3, cex.axis = 1.2)
  plotErrorCurve(abs(c(as.matrix(loess_err_075))), "orange", TRUE, xmax, lwd = 2)
  plotErrorCurve(abs(c(as.matrix(loess_err_003))), "blue", TRUE, xmax, lwd = 2)
  plotErrorCurve(abs(c(as.matrix(Hybrid_err))), "red", TRUE, xmax, lwd = 2)
  grid(col = "darkgray"); abline(v = pw*deltaTime/2, lty = "dashed"); legend("bottomright",  c("Linear alignment", "LOESS, default span = 0.75", "LOESS, optimized span = 0.03", "Chromatogram alignment"), col = c("green", "orange", "blue", "red"), lwd = 2, box.lwd = 1)

  # Figure 3B
  xmax <- 40
  plotErrorCurve(abs(c(as.matrix(loess_err_003[, strep0_10Pairs]))), "blue", FALSE, xmax, lty = 1, lwd = 1.5, main = "Change in fit with pair-run type", xlab = "Retention time difference (in sec)", ylab = "Cumulative fraction of peptides", cex.lab = 1.3, cex.axis = 1.2)
  plotErrorCurve(abs(c(as.matrix(loess_err_003[, strep10Pairs]))), "blue", TRUE, xmax, lty = 2, lwd = 1.5)
  plotErrorCurve(abs(c(as.matrix(loess_err_003[, strep0Pairs]))), "blue", TRUE, xmax, lty = 3, lwd = 1.5)
  plotErrorCurve(abs(c(as.matrix(Hybrid_err[, strep0_10Pairs]))), "red", TRUE, xmax, lty = 1, lwd = 1.5)
  plotErrorCurve(abs(c(as.matrix(Hybrid_err[, strep10Pairs]))), "red", TRUE, xmax, lty = 2, lwd = 1.5)
  plotErrorCurve(abs(c(as.matrix(Hybrid_err[, strep0Pairs]))), "red", TRUE, xmax, lty = 3, lwd = 1.5)
  grid(col = "darkgray"); legend("bottomright",  c("Chromatogram alignment Strep0", "Chromatogram alignment Strep10", "Chromatogram alignment Strep0_Strep10", "LOESS Strep0", "LOESS Strep10", "LOESS Strep0_Strep10"), col = c("red", "red", "red", "blue", "blue", "blue"), lty= c(3:1, 3:1), lwd = 1.5, box.lwd = 1)

  # Figure 3C
  names4ErrorStat <- c("MatchedWithinHalfPeak", "MatchedWithin1PW", "MatchedWithin2PW",
                       "MatchedWithinHalfPeakPct", "MatchedWithin1PeakPct", "MatchedWithin2PeakPct",
                       "MissingPeps")
  ErrorStatHybrid <- apply(Hybrid_err, MARGIN = 2,  statPerPair, analytes, pw, deltaTime)
  rownames(ErrorStatHybrid) <- names4ErrorStat
  ErrorStatHybrid <- as.data.frame(t(ErrorStatHybrid), stringsAsFactors = FALSE)
  ErrorStatHybrid[,1:7] <- lapply(ErrorStatHybrid[,1:7], as.numeric)

  ErrorStatLoess_003 <- apply(loess_err_003, MARGIN = 2,  statPerPair, analytes, pw, deltaTime)
  rownames(ErrorStatLoess_003) <- names4ErrorStat
  ErrorStatLoess_003 <- as.data.frame(t(ErrorStatLoess_003), stringsAsFactors = FALSE)
  ErrorStatLoess_003[,1:7] <- lapply(ErrorStatLoess_003[,1:7], as.numeric)

  hist(ErrorStatHybrid[,"MatchedWithinHalfPeak"], col = "red", ylim = c(0,50), xlim = c(320,430), breaks = seq(320, 430, 10), xlab = "Number of peaks aligned", main = "Histogram of number of aligned peaks in 120 pairs", cex.lab = 1.3, cex.axis = 1.2)
  hist(ErrorStatLoess_003[,"MatchedWithinHalfPeak"], col=rgb(0,0,1,0.5), breaks = seq(320, 430, 10), add = T)
  legend("topright",  c("Chromatrogram\nAlignment", "LOESS"), col = c("red", rgb(0,0,1,0.5)), lwd = 4, box.lwd = 1)

  # Figure 3D
  hist(c(as.matrix(Hybrid_err)), breaks = seq(-290,290,5), xlim = c(-75, 75), col = "red", xlab = "RT difference(sec)",  main = "Histogram of RT difference", cex.lab = 1.3, cex.axis = 1.2, axes = FALSE)
  axis(side = 1, at = seq(-75,75,25))
  axis(side = 2, at = seq(0,25000,5000))
  hist(c(as.matrix(loess_err_003)), col=rgb(0,0,1,0.5), breaks = seq(-290,290,5), add = T)
  grid(col = "darkgray")
  legend("topright",  c("Chromatogram\nAlignment", "LOESS"), col = c("red", rgb(0,0,1,0.5)), lwd = 4, box.lwd = 1)

  dev.off()
}


# Hybrid alignment
hybridAlignment <- function(StrepChromatogram, fileInfo, ObservedRT, precursors, pair_names,
                            globalFits, RSE, analytes, XICfilter = "sgolay", kernelLen = 9, polyOrd = 3,
                            RSEdistFactor = 3.5, alignType = "hybrid", normalization = "mean",
                            simMeasure = "dotProductMasked", goFactor = 0.125, geFactor = 40, cosAngleThresh = 0.3,
                            OverlapAlignment = TRUE, dotProdThresh = 0.96, gapQuantile = 0.5, kerLen = 9L, hardConstrain = FALSE,
                            samples4gradient = 100){
  analytes_chr <- as.character(analytes)
  hybrAlign <- list()
  message("Performing hybrid alignment.")
  for(pair in pair_names){
    ref <- strsplit(pair, split = "_")[[1]][1]
    eXp <- strsplit(pair, split = "_")[[1]][2]
    globalFit <- globalFits[[pair]]
    adaptiveRT <- RSEdistFactor*RSE[[pair]]
    MappedTime <- data.frame(analytes, stringsAsFactors=FALSE)
    colnames(MappedTime) <- c("transition_group_id")
    rownames(MappedTime) <- analytes_chr
    MappedTime$refRT <- ObservedRT[analytes_chr, ref]
    MappedTime$expRT <- ObservedRT[analytes_chr, eXp]
    XICs.A <- StrepChromatogram[[fileInfo[ref,"runName"]]]
    XICs.B <- StrepChromatogram[[fileInfo[eXp,"runName"]]]
    for(analyte in analytes_chr){
      if(!is.na(MappedTime[analyte, "refRT"])){
        XICs.ref.s <- smoothXICs(XICs.A[[analyte]], type = XICfilter, kernelLen = kernelLen, polyOrd = polyOrd)
        XICs.eXp.s <- smoothXICs(XICs.B[[analyte]], type = XICfilter, kernelLen = kernelLen, polyOrd = polyOrd)
        alignedTime  <- DIAlignR:::getMappedRT( ObservedRT[analyte, ref], XICs.ref.s, XICs.eXp.s, globalFit, alignType, adaptiveRT,
                                                normalization, simMeasure, goFactor, geFactor, cosAngleThresh,
                                                OverlapAlignment, dotProdThresh, gapQuantile, kerLen, hardConstrain,
                                                samples4gradient, objType = "light")
        MappedTime[analyte, "expRT_aligned"] <- alignedTime
      }
    }
    MappedTime$hybrid_err <- abs(MappedTime[, "expRT_aligned"] - ObservedRT[, eXp])
    hybrAlign[[pair]] <- MappedTime
    message(paste0("Aligned pair ", pair))
  }
  Hybrid_err  <- sapply(hybrAlign, `[[`, "hybrid_err")
  df <- as.data.frame(round(Hybrid_err, 2))
  group_label <- precursors$group_label[match(analytes, precursors$transition_group_id)]
  df <- cbind(analytes, group_label, df)
  write.table(format(df, scientific = FALSE), file="hybridAlignment.csv",
              sep = "\t", row.names = FALSE, quote = FALSE)
  Hybrid_err
}

#' DIAlignR paper results
#' @keywords internal
#' @examples
#' \dontrun{
#' peptidesFile <- "data/peptides.txt"
#' annotationFile <- "data/SkylineResult500Peptides.csv"
#' DIAlignR:::mcp_2019(annotationFile, peptidesFile)
#' }
mcp_2019 <- function(annotationFile, peptidesFile, dataPath = ".", oswMerged = TRUE,
                     globalAlignment = "loess", globalAlignmentFdr = 0.01, globalAlignmentSpan = 0.1,
                     XICfilter = "sgolay", kernelLen = 9, polyOrd = 3,
                     RSEdistFactor = 3.5, alignType = "hybrid", normalization = "mean",
                     simMeasure = "dotProductMasked", goFactor = 0.125, geFactor = 40, cosAngleThresh = 0.3,
                     OverlapAlignment = TRUE, dotProdThresh = 0.96, gapQuantile = 0.65, kerLen = 9L, hardConstrain = FALSE,
                     samples4gradient = 100){
  peptides <- readLines(peptidesFile)
  fileInfo <- getRunNames(dataPath, oswMerged)
  ManualAnnot500Peptides <- read.table(annotationFile, sep = ",", header = TRUE, na.strings = "#N/A", comment.char = "")

  ############# All run pairs ##############
  pair_names <- vector()
  for (i in 1:(nrow(fileInfo)-1)){
    for (j in (i+1): nrow(fileInfo)){
      pair_names <- c(paste(rownames(fileInfo)[i], rownames(fileInfo)[j], sep = "_"), pair_names)
    }}
  num_of_pair <- length(pair_names)
  strep0Runs <- rownames(fileInfo)[grep("Strep0", fileInfo$runName)]
  strep10Runs <- rownames(fileInfo)[grep("Strep10", fileInfo$runName)]
  strep0Pairs <- c(); strep10Pairs <- c(); strep0_10Pairs <- c()
  for(pair in pair_names){
    k <- sum(strsplit(pair, split = "_")[[1]] %in% strep10Runs)
    if(k == 2){strep10Pairs <- c(pair, strep10Pairs)}
    else if(k == 0){strep0Pairs <- c(pair, strep0Pairs)}
    else{strep0_10Pairs <- c(pair, strep0_10Pairs)}
  }

  ############# Run pairs colors ##############
  pair_color <- c()
  for(pair in pair_names){
    if(pair %in% strep0Pairs)pair_color <- c(pair_color, "strep0Pairs")
    if(pair %in% strep10Pairs)pair_color <- c(pair_color, "strep10Pairs")
    if(pair %in% strep0_10Pairs)pair_color <- c(pair_color, "strep0_10Pairs")
  }


  ######## Modify PTMs in Annotation file ##########
  AminoAcid <- c("A", "R", "D", "N", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  MassDa <- c(71, 156, 114, 115, 103, 129, 128, 57, 137, 113, 113, 128, 131, 147, 97, 87, 101, 186, 163, 99)
  names(MassDa) <- AminoAcid

  Peptide.Modified.Sequence.newlevels <- sapply(levels(ManualAnnot500Peptides$Peptide.Modified.Sequence), function(seq) {
    posStart <- stringi::stri_locate_all(pattern = '[', seq, fixed = TRUE)[[1]][,1]
    posEnd <- stringi::stri_locate_all(pattern = ']', seq, fixed = TRUE)[[1]][,1]
    outseq <- seq
    if(length(posStart) > 0){
      pattern <- "(.*)(\\[)(\\+|\\-)(.*?\\])(.*)"
      for(i in 1:length(posStart)) {
        newVal <- as.numeric(MassDa[substr(seq, posStart[i]-1, posStart[i]-1)]) + as.numeric(substr(seq, posStart[i]+1, posEnd[i]-1))
        outseq <- gsub(pattern, replacement = paste0("\\1[", newVal, "]\\5"), outseq)
      }
    }
    return(outseq)
  })
  Peptide.Modified.Sequence.new <- Peptide.Modified.Sequence.newlevels[ManualAnnot500Peptides$Peptide.Modified.Sequence]
  ManualAnnot500Peptides<- cbind(ManualAnnot500Peptides, Peptide.Modified.Sequence.new)
  transitionID <- paste(ManualAnnot500Peptides$Peptide.Modified.Sequence.new, ManualAnnot500Peptides$Precursor.Charge, sep = "/")
  ManualAnnot500Peptides<- cbind(ManualAnnot500Peptides, transitionID)

  ########## Get manually ObservedRT ##########
  num_of_runs <- nrow(fileInfo)
  # pattern <- "(.*)(_)([[:alpha:]]*)(\\[.*?\\])?(.*)(/)(.*)"
  pattern <- "(.*)(_)(.*)(/)(.*)"
  Peptide.Sequence <- sapply(peptides, function(x) as.character(gsub(pattern, replacement = "\\3", x)))
  # Peptide.Sequence <- sapply(Peptide.Sequence, function(x) gsub("[^[:alpha:] ]", replacement="", x))
  charge <- sapply(peptides, function(x) as.integer(gsub(pattern, replacement = "\\5", x)))

  pattern <- "(hroest_K12080[8,9]_)(.*?)(%PlasmaBiol)(.*?)"
  filenames <- sapply(fileInfo$runName, function(x) gsub(pattern, replacement = "\\2\\4", x))
  names(filenames) <- NULL
  ObservedRT <- cbind.data.frame(Peptide.Sequence, charge, stringsAsFactors= FALSE)
  for (run in filenames){
    RT <- 60*sapply(1:dim(ObservedRT)[1], function(i){
      ManualAnnot500Peptides$Peptide.Retention.Time[
        (ManualAnnot500Peptides$Peptide.Modified.Sequence.new == ObservedRT[i,1]) &
          (ManualAnnot500Peptides$Precursor.Charge == ObservedRT[i,2]) &
          (ManualAnnot500Peptides$Replicate.Name == run)][1]
    })
    ObservedRT <- cbind(ObservedRT, RT)}
  colnames(ObservedRT) <- c("Peptide.Sequence", "charge", rownames(fileInfo))
  ObservedRT["4205_SEASVGFGAAR/2", "run14"] <- 1980.7

  sum(rowSums(is.na(ObservedRT[,3:(num_of_runs+2)])) == num_of_runs)
  naPeptides <- ObservedRT[rowSums(is.na(ObservedRT[,3:(num_of_runs+2)])) == num_of_runs,c(1,2)]

  # 368 rows(transition IDs) don't have retention time value in Skyline output file. Because Hannes randomly selected the 500 peptides for manual validation.
  # We have total 816 transition IDs. Transition ID consists of peptide sequence(w/ modification) and charge state. Total we have 583 unique peptide sequences(w/ modification). Hence, 233 peptide sequences more than one charge state. There are 447 peptide sequences in Skyline result file. With available charge state, the number of transition IDs in Skyline result file is 452. Hence, we expect 816 - 452 = 364 transition IDs to have NA retention time in ObservedRT.
  # We see 368 - 364 = 4 peptides which have NA values for all runs. This means that in the Skyline file these 4 peptides have NA value associated with Peptide Retention time.

  # 4410_QEVNELSR/2
  # 4571_AIEEDGSIEIVTPDHK/3
  # 4815_TINNIIK/2
  # 469_YADSVPLVYEIASIPEK/2
  # These peptides have no retention time associated with them in skyline result file.

  # Add wrong annotated peptides in naPeptides
  incorrectPeps <- c("7118_MDELLDSLNANLDEM[147]SER/3", "6857_SEFPENELWDLTALYK/2",
                     "1288_FILTSDELFELK/3", "20672_TGLLAETFGGQVMETVGIENMIGTLYTEGPK/3",
                     "18091_Q[111]LDLLAHGER/2")
  ObservedRT <- ObservedRT[!(rownames(ObservedRT) %in% rownames(naPeptides)), ]

  ###### Remove peptides that are not in SCORE_PEPTIDE table. ##############
  precursors <- mcp_precursors(fileInfo)
  precursors <- precursors[precursors$group_label %in% rownames(ObservedRT),]
  ObservedRT <- ObservedRT[rownames(ObservedRT)%in%precursors$group_label,]
  rownames(ObservedRT) <- precursors$transition_group_id[match(rownames(ObservedRT), precursors$group_label)]
  rownames(ObservedRT) <- as.character(rownames(ObservedRT))
  analytes <- precursors$transition_group_id
  analytes_chr <- as.character(analytes)

  ########### Collect pointers for each mzML file. #######
  message("Collecting metadata from mzML files.")
  mzPntrs <- getMZMLpointers(fileInfo)
  message("Metadata is collected from mzML files.")

  ############# Get chromatogram Indices of precursors across all runs. ############
  message("Collecting chromatogram indices for all precursors.")
  prec2chromIndex <- getChromatogramIndices(fileInfo, precursors, mzPntrs)

  ########### Fetch all chromatograms. #################
  message("Fetching Extracted-ion chromatograms from runs")
  StrepChromatogram <- getXICs4AlignObj(mzPntrs, fileInfo, fileInfo[, "runName"],
                                        prec2chromIndex, analytes)
  rm(mzPntrs)


  ########## Remove analytes that have annotation outside of chromatograms. #########
  time <- StrepChromatogram[[fileInfo["run0", "runName"]]][[analytes_chr[1]]][[1]][["time"]]
  ChromLength <- time[length(time)] - time[1] # 594

  Chrom_Observed_TimeDiff <- sapply(colnames(ObservedRT)[3: ncol(ObservedRT)],
                                    function(run) {lapply(analytes_chr, function(pep)
                                      return(ObservedRT[pep, run] - StrepChromatogram[[fileInfo[run, "runName"]]][[pep]][[1]][[1,"time"]]))
                                    })
  rownames(Chrom_Observed_TimeDiff) <- analytes_chr

  outofChrom <- (Chrom_Observed_TimeDiff < 0 | Chrom_Observed_TimeDiff > ChromLength)
  outofChromCount <- apply(outofChrom, MARGIN =1, table, useNA = "always")
  outofChromSummary <- matrix(NA, nrow = length(analytes), ncol = 3)
  rownames(outofChromSummary) <- analytes_chr; colnames(outofChromSummary) <- c("TRUE", "FALSE", "NA")
  for (pep in analytes_chr) {
    if("TRUE" %in% names(outofChromCount[[pep]]))
      outofChromSummary[pep, "TRUE"] <- outofChromCount[[pep]][["TRUE"]]
    else
      outofChromSummary[pep, "TRUE"] <- 0
    if("FALSE" %in% names(outofChromCount[[pep]]))
      outofChromSummary[pep, "FALSE"] <- outofChromCount[[pep]][["FALSE"]]
    else
      outofChromSummary[pep, "FALSE"] <- 0
    outofChromSummary[pep, "NA"] <- as.numeric(tail(outofChromCount[[pep]], 1))
  }
  outofChromSummary[outofChromSummary[,"TRUE"] != 0,]
  outChromPeps <- rownames(outofChromSummary[outofChromSummary[,"TRUE"] != 0,])
  Chrom_Observed_TimeDiff[outChromPeps, 1:4]

  for (run in colnames(outofChrom)){
    for(pep in outChromPeps){
      if(!is.na(outofChrom[pep, run]))
        if(outofChrom[pep, run]){
          ObservedRT[pep, run] <- NA
        }
    }
  }
  removePeps <- analytes_chr[rowSums(is.na(ObservedRT[,3:(num_of_runs+2)])) >= (num_of_runs-1)]
  analytes_chr <- analytes_chr[!analytes_chr %in% removePeps]
  analytes <- as.integer(analytes_chr)
  ObservedRT <- ObservedRT[analytes_chr,]
  dim(ObservedRT); length(analytes_chr)
  write.table(format(ObservedRT, scientific = FALSE), file="ObservedRT.csv",
              sep = "\t", row.names = TRUE, quote = FALSE, col.names=NA)

  ######## Setting up peak width ##########
  leftWidth <- matrix(NA, nrow = nrow(ObservedRT), ncol = ncol(ObservedRT)-2)
  rownames(leftWidth) <- rownames(ObservedRT); colnames(leftWidth) <- colnames(ObservedRT[,3:18])
  rightWidth <- matrix(NA, nrow = nrow(ObservedRT), ncol = ncol(ObservedRT)-2)
  rownames(rightWidth) <- rownames(ObservedRT); colnames(rightWidth) <- colnames(ObservedRT[,3:18])
  for(i in 1:nrow(leftWidth)){
    for(j in 1:num_of_runs){
      leftWidth[i,j] <- 60*ManualAnnot500Peptides$Start.Time[
        (ManualAnnot500Peptides$Peptide.Modified.Sequence.new == ObservedRT[i,1]) &
          (ManualAnnot500Peptides$Precursor.Charge == ObservedRT[i,2]) &
          (ManualAnnot500Peptides$Replicate.Name == filenames[j])][1]
      rightWidth[i,j] <- 60*ManualAnnot500Peptides$End.Time[
        (ManualAnnot500Peptides$Peptide.Modified.Sequence.new == ObservedRT[i,1]) &
          (ManualAnnot500Peptides$Precursor.Charge == ObservedRT[i,2]) &
          (ManualAnnot500Peptides$Replicate.Name == filenames[j])][1]
    }
  }

  ######### get all features #############
  features <- getFeatures(fileInfo, maxFdrQuery = 1e-04, runType = "DIA_proteomics")

  ######### get all-pair global alignment ############
  globalFits <- list()
  RSE <- list()
  for(pair in pair_names){
    ref <- strsplit(pair, split = "_")[[1]][1]
    eXp <- strsplit(pair, split = "_")[[1]][2]
    globalFit <- getGlobalAlignment(features, ref, eXp, globalAlignment, globalAlignmentFdr, globalAlignmentSpan)
    globalFits[[pair]] <- globalFit
    RSE[[pair]] <- DIAlignR:::getRSE(globalFit)
  }

  ########### Hybrid RT alignment ###############
  Hybrid_err <- hybridAlignment(StrepChromatogram, fileInfo, ObservedRT, precursors, pair_names, globalFits, RSE, analytes,
                                XICfilter, kernelLen, polyOrd,
                                RSEdistFactor, alignType, normalization,
                                simMeasure, goFactor, geFactor, cosAngleThresh,
                                OverlapAlignment, dotProdThresh, gapQuantile,
                                kerLen, hardConstrain, samples4gradient)


  ############## Linear alignment ##############
  linearAlign <- list()
  message("Performing linear alignment")
  for(pair in pair_names){
    ref <- strsplit(pair, split = "_")[[1]][1]
    eXp <- strsplit(pair, split = "_")[[1]][2]
    globalFit <- getGlobalAlignment(features, ref, eXp, fitType = "linear", globalAlignmentFdr, globalAlignmentSpan)
    MappedTime <- data.frame(analytes, stringsAsFactors=FALSE)
    colnames(MappedTime) <- c("transition_group_id")
    rownames(MappedTime) <- analytes_chr
    MappedTime$refRT <- ObservedRT[analytes_chr, ref]
    MappedTime$expRT <- ObservedRT[analytes_chr, eXp]
    MappedTime[, "expRT_aligned"]  <- stats::predict(globalFit, data.frame("RT.ref" = MappedTime$refRT))
    MappedTime$global_rse <- DIAlignR:::getRSE(globalFit) # Residual Standard Error
    MappedTime$global_err <- MappedTime[, "expRT_aligned"]-  ObservedRT[, eXp]
    linearAlign[[pair]] <- MappedTime
  }
  linear_err  <- sapply(linearAlign, `[[`, "global_err")
  message("Linear alignment completed")

  ######### Loess alignment with deault span ###################
  loessDefault <- list()
  message("Performing loess alignment with default span value 0.75")
  for(pair in pair_names){
    ref <- strsplit(pair, split = "_")[[1]][1]
    eXp <- strsplit(pair, split = "_")[[1]][2]
    globalFit <- getGlobalAlignment(features, ref, eXp, fitType = "loess", globalAlignmentFdr, 0.75)
    MappedTime <- data.frame(analytes, stringsAsFactors=FALSE)
    colnames(MappedTime) <- c("transition_group_id")
    rownames(MappedTime) <- analytes_chr
    MappedTime$refRT <- ObservedRT[analytes_chr, ref]
    MappedTime$expRT <- ObservedRT[analytes_chr, eXp]
    MappedTime[, "expRT_aligned"]  <- stats::predict(globalFit, data.frame("RT.ref" = MappedTime$refRT))
    MappedTime$global_rse <- DIAlignR:::getRSE(globalFit) # Residual Standard Error
    MappedTime$global_err <- MappedTime[, "expRT_aligned"]-  ObservedRT[, eXp]
    loessDefault[[pair]] <- MappedTime
  }
  loess_err_075  <- sapply(loessDefault, `[[`, "global_err")
  message("Loess alignment completed")

  ######## Loess alignment with 0.03 spanvalue ###########
  loessOptimum <- list()
  message("Performing loess alignment with default span value 0.03")
  for(pair in pair_names){
    ref <- strsplit(pair, split = "_")[[1]][1]
    eXp <- strsplit(pair, split = "_")[[1]][2]
    globalFit <- getGlobalAlignment(features, ref, eXp, fitType = "loess", globalAlignmentFdr, 0.03)
    MappedTime <- data.frame(analytes, stringsAsFactors=FALSE)
    colnames(MappedTime) <- c("transition_group_id")
    rownames(MappedTime) <- analytes_chr
    MappedTime$refRT <- ObservedRT[analytes_chr, ref]
    MappedTime$expRT <- ObservedRT[analytes_chr, eXp]
    MappedTime[, "expRT_aligned"]  <- stats::predict(globalFit, data.frame("RT.ref" = MappedTime$refRT))
    MappedTime$global_rse <- DIAlignR:::getRSE(globalFit) # Residual Standard Error
    MappedTime$global_err <- MappedTime[, "expRT_aligned"]-  ObservedRT[, eXp]
    loessOptimum[[pair]] <- MappedTime
  }
  loess_err_003  <- sapply(loessOptimum, `[[`, "global_err")
  message("Loess alignment completed")

  ######### Get figure 3 of the paper ################
  figure3(linear_err, loess_err_075, loess_err_003, Hybrid_err,
          analytes, strep0Pairs, strep10Pairs, strep0_10Pairs)

  names4ErrorStat <- c("MatchedWithinHalfPeak", "MatchedWithin1PW", "MatchedWithin2PW",
                       "MatchedWithinHalfPeakPct", "MatchedWithin1PeakPct", "MatchedWithin2PeakPct",
                       "MissingPeps")
  ErrorStatHybrid <- apply(Hybrid_err, MARGIN = 2,  statPerPair, analytes, pw=9, deltaTime= 3.4)
  rownames(ErrorStatHybrid) <- names4ErrorStat
  ErrorStatHybrid <- as.data.frame(t(ErrorStatHybrid), stringsAsFactors = FALSE)
  ErrorStatHybrid[,1:7] <- lapply(ErrorStatHybrid[,1:7], as.numeric)
  print(colMeans(ErrorStatHybrid))
}


