#!/usr/bin/env Rscript

rm(list=ls())
baseDir <- "/media/mwschmid/myData/MWSchmid/Ivana_ATAC"
infileName <- file.path(baseDir, "GitIgnore_results", "union.3.BP.onlyTimePointOne.Rdata")
outfilePrefix <- file.path(baseDir, "report", "forFigureEditing", "boxplotCandidateTP1")
peakID <- "chr4_154654980_154655784"
infileName <- file.path(baseDir, "GitIgnore_results", "union.3.BP.onlyTimePointTwo.Rdata")
outfilePrefix <- file.path(baseDir, "report", "forFigureEditing", "boxplotCandidateTP2")
peakID <- "chr4_154654980_154655784"
rm(list=ls())

## arguments from commandline
myarg <- commandArgs()
argPos <- grep("--args", myarg, fixed = TRUE)
infileName <- myarg[argPos+1]
outfilePrefix <- myarg[argPos+2]
peakID <- myarg[argPos+3]

# rest
suppressPackageStartupMessages({
  options(java.parameters = "-Xmx12g")
  library("RColorBrewer")
  library("DESeq2")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapperAds.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper_overviewBandwidthFix.R")
})

f.replace.labels.with.lab <- function(x) {
  if (grepl("Hannover1", x)) { return("LAB1") }
  if (grepl("Hannover2", x)) { return("LAB2") }
  if (grepl("Bern", x)) { return("LAB3") }
  if (grepl("Munster", x)) { return("LAB4") }
  if (grepl("Zurich", x)) { return("LAB5") }
  return("ERROR")
}

##############################################################################
### load file
load(infileName)

forPlot <- data.frame(
  accVal = myNormData[peakID,rownames(sampleTab)],
  time = sampleTab$time,
  lab = factor(sapply(sampleTab$lab, f.replace.labels.with.lab), levels = paste0("LAB", 1:5)),
  color = sampleTab$color,
  stringsAsFactors = FALSE
)
forPlot <- forPlot[order(forPlot$lab),]

##############################################################################
### boxplot
forColor <- unique(forPlot[,c("lab", "color")])$color
pdf(paste0(outfilePrefix, "_boxplot.pdf"), height = 5, width = 3)
boxplot(accVal ~ lab, data = forPlot, las = 1, xlab = "", ylab = "Accessibility (log2(normalized counts+1))", col = forColor)
points(jitter(as.numeric(as.factor(forPlot$lab))), forPlot$accVal, col = "black", pch = 16, cex = 0.5)
invisible(dev.off())


