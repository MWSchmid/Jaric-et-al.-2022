#!/usr/bin/env Rscript

rm(list=ls())
baseDir <- "/media/mwschmid/myData/MWSchmid/Ivana_ATAC"
infileName <- file.path(baseDir, "GitIgnore_results", "union.5.BP.onlyTimePointOne.Rdata")
outfilePrefix <- file.path(baseDir, "GitIgnore_results", "union.5.BP.onlyTimePointOne")
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
load(infileName) # sampleTab, deResults, deResultsNoBatch, myData, myNormData

##############################################################################
### check correlation in lfc and number of sig genes
deResults <- deResults[names(deResultsNoBatch)]

geneNumbers <- matrix(NA, nrow = length(names(deResults)), ncol = 5, dimnames = list(names(deResults), c("union", "percCommon", "percSpecificNoBatch", "percSpecificWithBatch", "lfcCorrelation")))
for (curCont in names(deResults)) {
  sigPeaksNoBatch <- deResultsNoBatch[[curCont]]$get_significant_entries(0.01, 0.01, 0.5)
  sigPeaks <- deResults[[curCont]]$get_significant_entries(0.01, 0.01, 0.5)
  for (curSet in names(sigPeaksNoBatch)) {
    merged <- union(sigPeaksNoBatch[[curSet]], sigPeaks[[curSet]])
    common <- intersect(sigPeaksNoBatch[[curSet]], sigPeaks[[curSet]])
    specificNoBatch <- setdiff(sigPeaksNoBatch[[curSet]], common)
    specificWithBatch <- setdiff(sigPeaks[[curSet]], common)
    total <- length(merged)
    sameRows <- intersect(rownames(deResultsNoBatch[[curCont]]$table), rownames(deResults[[curCont]]$table))
    lfcCorrelation <- cor(cbind(deResultsNoBatch[[curCont]]$table[sameRows, "logFC"], deResults[[curCont]]$table[sameRows, "logFC"]), method = "pearson")
    geneNumbers[curCont,] <- c(total, round(c(length(common), length(specificNoBatch), length(specificWithBatch))/total*100, 2), lfcCorrelation[1,2])
  }
}

write.csv(geneNumbers, paste0(outfilePrefix, ".batchVsNoBatch.csv"), quote = FALSE)

