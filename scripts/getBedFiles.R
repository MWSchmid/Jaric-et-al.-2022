#!/usr/bin/env Rscript

rm(list=ls())
baseDir <- "/media/mwschmid/myData/MWSchmid/Ivana_ATAC"
infileName <- file.path(baseDir, "GitIgnore_results", "union.5.BP.bothTimePoints.Rdata")
outfilePrefix <- file.path(baseDir, "GitIgnore_results", "bedFiles", "union.5.BP.bothTimePoints")
rm(list=ls())

## arguments from commandline
myarg <- commandArgs()
argPos <- grep("--args", myarg, fixed = TRUE)
infileName <- myarg[argPos+1]
outfilePrefix <- myarg[argPos+2]

# rest
suppressPackageStartupMessages({
  options(java.parameters = "-Xmx12g")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapperAds.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper_overviewBandwidthFix.R")
})

##############################################################################
### load file
load(infileName) # sampleTab, deResults, (deResultsNoBatch), myData, myNormData, (design, myCont)

##############################################################################
### collect significant peaks
contrastList <- list(
  oto = grep("^oto", names(deResults), value = TRUE),
  otm = grep("^otm", names(deResults), value = TRUE),
  time = grep("^time", names(deResults), value = TRUE),
  lrt = grep("^LRT", names(deResults), value = TRUE)
)

for (curSet in names(contrastList)) {
  if (length(contrastList[[curSet]]) > 0) {
    tempPeakList <- list(
      sigPeaks_LFC_1.0 = Reduce(union, lapply(deResults[contrastList[[curSet]]], function(x) x$get_significant_entries(0.01, 0.01, 1.0)$any)),
      sigPeaks_LFC_1.5 = Reduce(union, lapply(deResults[contrastList[[curSet]]], function(x) x$get_significant_entries(0.01, 0.01, 1.5)$any)),
      sigPeaks_LFC_2.0 = Reduce(union, lapply(deResults[contrastList[[curSet]]], function(x) x$get_significant_entries(0.01, 0.01, 2.0)$any))
    )
    tempBedList <- lapply(tempPeakList, function(x) as.data.frame(do.call("rbind", strsplit(x, '_')), stringsAsFactors = FALSE))
    for (curPeaks in names(tempBedList)) {
      if (nrow(tempBedList[[curPeaks]]) == 0) { next }
      colnames(tempBedList[[curPeaks]]) <- c("chrom", "start", "end")
      tempBedList[[curPeaks]]$name <- paste0(curSet, '_', curPeaks, '_', 1:nrow(tempBedList[[curPeaks]]))
      write.table(tempBedList[[curPeaks]], paste0(outfilePrefix, '.', curSet, '.', curPeaks, ".bed"), row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
    }
  }
}

