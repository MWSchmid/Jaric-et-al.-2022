#!/usr/bin/env Rscript

rm(list=ls())
baseDir <- "/media/mwschmid/myData/MWSchmid/Ivana_ATAC"
infileName <- file.path(baseDir, "GitIgnore_results", "union.5.BP.bothTimePoints.Rdata")
outfilePrefix <- file.path(baseDir, "GitIgnore_results", "geneSets", "union.5.BP.bothTimePoints")
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
### sort the contrasts
allContrasts <- c(grep("^oto", names(deResults), value = TRUE),
                  grep("^otm", names(deResults), value = TRUE),
                  grep("^time", names(deResults), value = TRUE),
                  grep("^LRT", names(deResults), value = TRUE))

deResults <- deResults[allContrasts]
if (exists("deResultsNoBatch")) {
  deResultsNoBatch <- deResultsNoBatch[allContrasts]
}

##############################################################################
### write an excel table only with the overview page
f.write.DEGtabs.to.workbook(deResults, dirname(outfilePrefix), basename(outfilePrefix), onlyOverviewSheet = TRUE); gc()
if (exists("deResultsNoBatch")) {
  f.write.DEGtabs.to.workbook(deResultsNoBatch, dirname(outfilePrefix), paste0(basename(outfilePrefix), ".noBatch"), onlyOverviewSheet = TRUE); gc()
}
