#!/usr/bin/env Rscript

rm(list=ls())
baseDir <- "/media/mwschmid/myData/MWSchmid/Ivana_ATAC"
infileName <- file.path("/home/marc/tempIJA", "union.3.BP.csv")
outfilePrefix <- file.path(baseDir, "GitIgnore_peaks_merged", "union.3.BP.peakChar")
rm(list=ls())

## arguments from commandline
myarg <- commandArgs()
argPos <- grep("--args", myarg, fixed = TRUE)
infileName <- myarg[argPos+1]
outfilePrefix <- myarg[argPos+2]

# rest
suppressPackageStartupMessages({
  options(java.parameters = "-Xmx12g")
  library("ChIPseeker")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapperAds.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper_overviewBandwidthFix.R")
})

##############################################################################
### load file
myData <- read.csv(infileName, header = TRUE, row.names = 1, stringsAsFactors = FALSE)

#########################################################################################
# plot dist to TSS, plut minus 100 kb
totPeaks <- nrow(myData)
forPlot <- subset(myData, abs(distToTSS) < 100001)
percWithin100kb <- round(nrow(forPlot)/totPeaks*100, 2)
distToTSSdensity <- density(forPlot$distToTSS)
upstreamQuantiles <- quantile(forPlot$distToTSS[forPlot$distToTSS<0])
downstreamQuantiles <- quantile(forPlot$distToTSS[forPlot$distToTSS>0])
svg(paste0(outfilePrefix, ".distToTSS.svg"), height = 4, width = 6)
plot(distToTSSdensity, xlab = "Distance to TSS", ylab = "Peak density", main = paste0(percWithin100kb, " % are within 100 kb of the TSS"), lwd = 2)
abline(v = upstreamQuantiles)
abline(v = downstreamQuantiles)
invisible(dev.off())

#########################################################################################
# check if the simpleAnno column is presend
if (!("simpleAnno" %in% colnames(myData))) {
  f.print.message("No simpleAnno available, exiting...")
  quit("no", 0)
}

#########################################################################################
# stacked bar plots with percentages
annotationOrder <- c("USTSS2", "USTSS1", "DSTSS1", "DSTSS2", "5pUTR", "firstExon", "firstIntron", "exon", "intron", "3pUTR", "DS1", "DS2", "distalIntergenic")
forColor <- c("#BCBDDC", "#756BB1", "#3182BD", "#9ECAE1", "#719579", "#F06C45", "#E78E79", "#FE982C", "#EFBEB4", "#F16767", "#F0F0F0", "#BDBDBD", "#636363")
names(forColor) <- annotationOrder
featureCounts <- table(myData$simpleAnno)
featureCounts <- featureCounts[annotationOrder[annotationOrder %in% names(featureCounts)]]
featurePerc <- round(featureCounts/sum(featureCounts)*100, 2)
forPlot <- matrix(featurePerc, ncol = 1, nrow = length(featurePerc), dimnames = list(names(featurePerc), "allPeaks"))
svg(paste0(outfilePrefix, ".genomicFeatures.svg"), height = 6, width = 8)
par(mfrow = c(1,2))
barplot(forPlot, col = forColor[rownames(forPlot)], las = 1, ylab = "% peaks in genomic feature")
plot(c(0,0), type = "n", xaxt = "n", yaxt = "n", main = "", xlab = "", ylab = "", bty = "n")
legend("center", legend = rev(names(forColor)), fill = rev(forColor))
invisible(dev.off())

#summary(myData$distToTSS[myData$simpleAnno == "prom_1"])
#summary(myData$distToTSS[myData$simpleAnno == "prom_2"])
