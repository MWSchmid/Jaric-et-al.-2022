#!/usr/bin/env Rscript

## arguments from commandline
myarg <- commandArgs()
argPos <- grep("--args", myarg, fixed = TRUE)

f.print.message <- function(...) { cat("=== ", format(Sys.time(), "%Y %b %d %X"), paste("===", ...,"\n")) }

suppressPackageStartupMessages({
  library("RColorBrewer")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapperAds.R")
}
)

# some directories
flankSize <- as.numeric(myarg[argPos+1])
sampleAndFileTabFile <- myarg[argPos+2]
geneAnnoFile <- myarg[argPos+3]
rDir <- myarg[argPos+4]

useSamples <- FALSE
if ("--useSamples" %in% myarg) {
  useSamples <- TRUE
}

useLog <- FALSE
if ("--useLog" %in% myarg) {
  useLog <- TRUE
}

hasMarkers <- FALSE
if ("--markers" %in% myarg) {
  hasMarkers <- TRUE
  markerFile <- myarg[which("--markers" == myarg) + 1]
}

#flankSize <- 2000
#sampleAndFileTabFile <- "/media/mwschmid/myData/MWSchmid/Ivana_ATAC/GitIgnore_coveragePlots/sampleAndFileTab.csv"
#geneAnnoFile <- "/media/mwschmid/myData/MWSchmid/Ivana_ATAC/GitIgnore_coveragePlots/forCoveragePlot.txt"
#rDir <- "/media/mwschmid/myData/MWSchmid/Ivana_ATAC/GitIgnore_coveragePlots/both"

#flankSize <- 0
#sampleAndFileTabFile <- "/media/mwschmid/myData/MWSchmid/MarkusStoffel_ATAC/GitIgnore_covPlots/sampleAndFileTabRegion.csv"
#geneAnnoFile <- "/media/mwschmid/myData/MWSchmid/MarkusStoffel_ATAC/GitIgnore_covPlots/regionExample.txt"
#rDir <- "/media/mwschmid/myData/MWSchmid/MarkusStoffel_ATAC/GitIgnore_covPlots"

#flankSize <- 1000
#markerFile <- "/media/mwschmid/myData/MWSchmid/MarkusStoffelDiffSplice/report/kostasKO/INTERSECT_markersForCovPlot.csv"
#sampleAndFileTabFile <- "/media/mwschmid/myData/MWSchmid/MarkusStoffelDiffSplice/GitIgnore_covPlots/INTERSECT_sampleAndFileTabEachTranscript.csv"
#geneAnnoFile <- "/media/mwschmid/myData/MWSchmid/MarkusStoffelDiffSplice/GitIgnore_covPlots/INTERSECT_forCoveragePlotEachTranscript.txt"
#rDir <- "/media/mwschmid/myData/MWSchmid/MarkusStoffelDiffSplice/GitIgnore_covPlots"

########################################################################
# load data
if (hasMarkers) {
  f.print.message("using markers from:", markerFile)
  markerCols <- c("ensemblID", "symbol", "markerChrom", "markerStart", "markerEnd", "markerType", "markerDirection")
  f.print.message("expecting the following columns\n", paste0(markerCols, collapse= '\n'))
  markers <- read.csv(markerFile, stringsAsFactors = FALSE)
  f.print.message("found the following columns\n", paste0(colnames(markers), collapse= '\n'))
  colnames(markers)[1:length(markerCols)] <- markerCols
  if ("simpleMarkerType" %in% colnames(markers)) {
    f.print.message("Using simple markers instead.")
    markers$markerType <- markers$simpleMarkerType
  }
  f.print.message("Choosing colors.")
  markerTypes <- unique(markers$markerType)
  markerCols <- brewer.pal(length(markerTypes), "Set3"); names(markerCols) <- markerTypes
  markerPchs <- rep(16, length(markerTypes)); names(markerPchs) <- markerTypes
  f.open.figure(file.path(rDir), "0_colorLegend.svg", TRUE, height = 15, width = 10)
  par(mfrow = c(1, 1))
  f.plot.legend(markerCols, markerPchs)
  f.close.figure()
}

geneAnno <- read.table(geneAnnoFile, header = FALSE, stringsAsFactors = FALSE, col.names = c("gene", "transcript", "chrom", "strand", "start", "end", "exon_starts", "exon_ends"))
rownames(geneAnno) <- paste0(geneAnno$gene, "_", geneAnno$transcript)
geneAnno$exon_starts <- as.character(geneAnno$exon_starts)
geneAnno$exon_ends <- as.character(geneAnno$exon_ends)
sampleAndFileTab <- read.csv(sampleAndFileTabFile, stringsAsFactors = FALSE, header = TRUE, row.names = 1)
geneToMaxCov <- list()
genesAndSampleCoverage <- list()
for (curSample in rownames(sampleAndFileTab)) {
  cat(curSample, '\n')
  curFile <- sampleAndFileTab[curSample,"file"]
  allLines <- readLines(curFile)
  for (curLine in allLines) { # well, there's no header [2:length(allLines)]
    fields <- unlist(strsplit(curLine, "\t"))
    gene <- fields[1]
    cat(gene, '\n')
    if (!(gene %in% names(genesAndSampleCoverage))) {
      genesAndSampleCoverage[[gene]] <- list()
      geneToMaxCov[[gene]] <- 0
    }
    if (useLog) {
      cov <- log2(as.numeric(unlist(strsplit(fields[4], ';')))+1)
    } else {
      cov <- as.numeric(unlist(strsplit(fields[4], ';')))
    }
    maxCov <- max(cov)
    if (maxCov > geneToMaxCov[[gene]]) {
      geneToMaxCov[[gene]] <- maxCov
    }
    genesAndSampleCoverage[[gene]][[curSample]] <- cov
    #if (gene == "Pcsk2") {
    #  break
    #}
  }
}

genesAndGroupCoverage <- list()
temp <- unique(sampleAndFileTab[,c("group", "lineType", "lineColor")])
groupToCol <- temp[,c("lineType", "lineColor")]
rownames(groupToCol) <- temp[, "group"]
for (gene in names(genesAndSampleCoverage)) {
  genesAndGroupCoverage[[gene]] <- list()
  for (curGroup in sampleAndFileTab$group) {
    curSamples <- rownames(sampleAndFileTab)[sampleAndFileTab$group == curGroup]
    if (sum(sampleAndFileTab$group == curGroup) > 1) {
      temp <- do.call("cbind", genesAndSampleCoverage[[gene]][curSamples])
      genesAndGroupCoverage[[gene]][[curGroup]] <- rowMeans(temp)
    } else {
      genesAndGroupCoverage[[gene]][[curGroup]] <- unlist(genesAndSampleCoverage[[gene]][curSamples])
    }
  }
}

########################################################################
# plot each gene separately
for (curGene in names(genesAndSampleCoverage)) {
  isFirst <- TRUE
  #fileName <- paste0(curGene, ".png")
  #png(file.path(rDir, fileName), height = 1000, width = 5000, pointsize = 20)
  curAnno <- subset(geneAnno, gene == curGene)
  numTrans <- nrow(curAnno)
  bottomExpand <- numTrans
  if (hasMarkers) {
    curMarkerDat <- subset(markers, symbol == curGene)
    bottomExpand <- bottomExpand + length(unique(curMarkerDat$markerType))
  }
  maxCov <- geneToMaxCov[[gene]]
  f.print.message("plotting", curGene)
  fileName <- paste0(curGene, ".svg")
  svg(file.path(rDir, fileName), height = 6, width = 20)
  if (useSamples) {
    for (curSample in names(genesAndSampleCoverage[[curGene]])) {
      curCov <- genesAndSampleCoverage[[curGene]][[curSample]]
      if (isFirst) {
        isFirst <- FALSE
        plot(1:length(curCov), curCov, ylim = c((-maxCov/10)*bottomExpand, geneToMaxCov[[curGene]]), main = curGene, ylab = "coverage", xlab = paste0("gene and ", round(flankSize/1000, 1), " kb flanking regions"), type = 'l', lty = sampleAndFileTab[curSample,"lineType"], lwd = 3, col = sampleAndFileTab[curSample,"lineColor"], bty = 'n', las = 1)
        # this was for Plotting all equal with Felipe
        #plot(1:length(curCov), curCov, ylim = c(yAbsMin, 10), main = curGene, ylab = "coverage", xlab = paste0("gene and ", round(flankSize/1000, 1), " kb flanking regions"), type = 'l', lty = sampleAndFileTab[curSample,"lineType"], lwd = 3, col = sampleAndFileTab[curSample,"lineColor"], bty = 'n', las = 1)
      } else {
        lines(1:length(curCov), curCov, lty = sampleAndFileTab[curSample,"lineType"], lwd = 3, col = sampleAndFileTab[curSample,"lineColor"])
      }
    }
  } else {
    for (curGroup in names(genesAndGroupCoverage[[curGene]])) {
      curCov <- genesAndGroupCoverage[[curGene]][[curGroup]]
      if (isFirst) {
        isFirst <- FALSE
        plot(1:length(curCov), curCov, ylim = c((-maxCov/10)*bottomExpand, geneToMaxCov[[curGene]]), main = curGene, ylab = "coverage", xlab = paste0("gene and ", round(flankSize/1000, 1), " kb flanking regions"), type = 'l', lty = groupToCol[curGroup,"lineType"], lwd = 3, col = groupToCol[curGroup,"lineColor"], bty = 'n', las = 1)
        # this was for Plotting all equal with Felipe
        #plot(1:length(curCov), curCov, ylim = c(yAbsMin, 10), main = curGene, ylab = "coverage", xlab = paste0("gene and ", round(flankSize/1000, 1), " kb flanking regions"), type = 'l', lty = groupToCol[curGroup,"lineType"], lwd = 3, col = groupToCol[curGroup,"lineColor"], bty = 'n', las = 1)
      } else {
        lines(1:length(curCov), curCov, lty = groupToCol[curGroup,"lineType"], lwd = 3, col = groupToCol[curGroup, "lineColor"])
      }
    }
    
  }
  # indicate exons with lines
  uniStrand <- unique(curAnno$strand)
  if (length(uniStrand) != 1) {
    f.print.message(curGene, "has transcripts on different strands")
    isPlus <- (curAnno$strand[1] == '+') # it should be like the first in the table
  } else {
    isPlus <- (uniStrand == '+')
  }
  for (i in 1:numTrans) {
    exonStartEnds <- cbind(as.numeric(unlist(strsplit(curAnno[i, "exon_starts"], ','))), as.numeric(unlist(strsplit(curAnno[i, "exon_ends"], ','))))
    if (isPlus) {
      exonStartEnds <- exonStartEnds - curAnno$start[1] + flankSize
    } else {
      exonStartEnds <- abs(exonStartEnds - curAnno$end[1]) + flankSize
      temp <- cbind(exonStartEnds[,2], exonStartEnds[,1])
      exonStartEnds <- temp
    }
    yPos <- (-maxCov/10)*(i-1)-maxCov/20
    apply(exonStartEnds, 1, function(x) lines(x, c(yPos,yPos), lwd = 5, col = "black"))
  }
  # add markers if given
  if (hasMarkers) {
    markerTypes <- unique(curMarkerDat$markerType)
    for (i in 1:length(markerTypes)) {
      subDat <- subset(curMarkerDat, markerType == markerTypes[i])
      yPos <- (-maxCov/10)*(numTrans+i-1)-maxCov/20
      markerStartEnds <- subDat[,c("markerStart", "markerEnd")]
      if (isPlus) {
        markerStartEnds <- markerStartEnds - curAnno$start[1] + flankSize
      } else {
        markerStartEnds <- abs(markerStartEnds - curAnno$end[1]) + flankSize
        temp <- cbind(markerStartEnds[,2], markerStartEnds[,1])
        markerStartEnds <- temp
      }
      apply(markerStartEnds, 1, function(x) lines(x, c(yPos,yPos), lwd = 5, col = markerCols[markerTypes[i]]))
    }
  }
  invisible(dev.off())
}

