#!/usr/bin/env Rscript

rm(list=ls())
baseDir <- "/media/mwschmid/myData/MWSchmid/Ivana_ATAC"
infileName <- file.path(baseDir, "GitIgnore_results", "geneSets", "union.3.BP.onlyTimePointOne.minLFC_0.5.numberSigPeaks.csv.gz")
outfilePrefix <- file.path(baseDir, "report", "forFigureEditing", "peakCounts_TP1")
infileName <- file.path(baseDir, "GitIgnore_results", "geneSets", "union.3.BP.onlyTimePointTwo.minLFC_0.5.numberSigPeaks.csv.gz")
outfilePrefix <- file.path(baseDir, "report", "forFigureEditing", "peakCounts_TP2")
rm(list=ls())

## arguments from commandline
myarg <- commandArgs()
argPos <- grep("--args", myarg, fixed = TRUE)
infileName <- myarg[argPos+1]
outfilePrefix <- myarg[argPos+2]

# rest
suppressPackageStartupMessages({
  options(java.parameters = "-Xmx12g")
  library("RColorBrewer")
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
myData <- read.csv(gzfile(infileName), header = TRUE, stringsAsFactors = FALSE, row.names = 1)
otoRows <- grep("^oto", rownames(myData), value = TRUE)
otmRows <- grep("^otm", rownames(myData), value = TRUE)

##############################################################################
### barplot otm
otmData <- t(as.matrix(myData[otmRows, c("up", "down")]))
colnames(otmData) <- sapply(colnames(otmData), f.replace.labels.with.lab)
otmData <- otmData[,sort(colnames(otmData))]
upDownColors <- c("slateblue3", "orangered3")
pdf(paste0(outfilePrefix, "_otm.pdf"), width = 5, height = 6)
barplot(otmData, col = upDownColors, las = 1)
legend("topright", legend = rownames(otmData), cex = 1, bty = 'n', pch = 15, pt.cex = 2, col = upDownColors)
invisible(dev.off())

##############################################################################
### matrix with sig peaks
allLabs <- paste0("LAB", 1:5)
otoData <- matrix(NA, nrow = length(allLabs), ncol = length(allLabs), dimnames = list(allLabs, allLabs))
for (curRow in otoRows) {
  temp <- unlist(strsplit(gsub("oto_T[[:digit:]]_", "", curRow), "_vs_"))
  leftLab <- f.replace.labels.with.lab(temp[1])
  rightLab <- f.replace.labels.with.lab(temp[2])
  upInLeftLab <- myData[curRow, "up"]
  upInRightLab <- myData[curRow, "down"]
  otoData[leftLab, rightLab] <- upInLeftLab
  otoData[rightLab, leftLab] <- upInRightLab
}
otoDataForPlot <- t(otoData[rev(rownames(otoData)),])
pdf(paste0(outfilePrefix, "_oto.pdf"), width = 5, height = 6)
f.image.with.text(rownames(otoDataForPlot), colnames(otoDataForPlot), otoDataForPlot, "downregulated", "upregulated", "", useLog = TRUE, col = brewer.pal(9, "Greens")[1:8])#heat.colors(31))
invisible(dev.off())
# from left to right is upregulated in the labelled lab
# from top to bottom is downregulated