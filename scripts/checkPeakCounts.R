rm(list=ls())
baseDir <- "/media/mwschmid/myData/MWSchmid/Ivana_ATAC"
scriptDir <- file.path(baseDir, "scripts")
designFile <- file.path(baseDir, "report", "sampleTab.csv")
infileName <- file.path(baseDir, "GitIgnore_results", "peakCounts.csv")
outfilePrefix <- file.path(baseDir, "report", "peakCounts")
rm(list=ls())

## arguments from commandline
myarg <- commandArgs()
argPos <- grep("--args", myarg, fixed = TRUE)
designFile <- myarg[argPos+1]
infileName <- myarg[argPos+2]
outfilePrefix <- myarg[argPos+3]

# rest
suppressPackageStartupMessages({
  options(java.parameters = "-Xmx12g")
  library("RColorBrewer")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapperAds.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper_overviewBandwidthFix.R")
})

## a function for printing
f.print.message <- function(...) { cat("=== ", format(Sys.time(), "%Y %b %d %X"), paste("===", ...,"\n")) }

rDir <- dirname(outfilePrefix)
filePrefixNoPath <- basename(outfilePrefix)

#########################################################################################
# load data
sampleTab <- read.csv(designFile, stringsAsFactors = FALSE, row.names = 1)
myData <- read.csv(infileName, stringsAsFactors = FALSE, header = TRUE)
metrics <- c("numPeaks")
myData$time <- gsub("\\.[[:alnum:]]{1,11}$", "", myData$timePointAndLab)
myData$lab <- gsub("[[:alnum:]]{2}\\.", "", myData$timePointAndLab)

#########################################################################################
# test for differences, use different peak types ar replicate
matrixCols <- paste0(c("percSS", "pValue"), '_', rep(c("time", "lab", "peakInAtLeastXSamples"), each = 2))
allTests <- matrix(NA, ncol = length(matrixCols), nrow = length(metrics), dimnames = list(metrics, matrixCols))
for (curMetric in metrics) {
  res <- anova(lm(terms(as.formula(paste0(curMetric, '~ time + lab + peakInAtLeastXSamples')), keep.order = TRUE), data = myData))
  allTerms <- rownames(res)[1:3]
  allTests[curMetric, paste0("percSS_", allTerms)] <- round(res[allTerms, "Sum Sq"]/sum(res[, "Sum Sq"])*100, 2)
  allTests[curMetric, paste0("pValue_", allTerms)] <- round(res[allTerms, "Pr(>F)"], 5)
}
write.csv(allTests, paste0(outfilePrefix, ".anova.csv"))

#########################################################################################
# plots, use predefined colors and pchs
#allGroups <- unique(myData$timePointAndLab)
#forColor <- brewer.pal(length(allGroups), "Paired"); names(forColor) <- paste0(rep(c("T1.", "T2."), 5), rep(c("Bern", "Hannover1", "Hannover2", "Munster", "Zurich"), each = 2))
#forPch <- c("T1" = 15, "T2" = 16)
#myData$color <- forColor[myData$timePointAndLab]
#myData$pch <- forPch[myData$time]

#########################################################################################
### PCA, doesn't make sense with one dimension

quit("no", 0)

pValCols <- grep("^pValue", colnames(allTests), value = TRUE)
if (sum(allTests[,pValCols] < 0.01) == 0) {
  f.print.message("no differences at all")
  quit("no", 0)
}

# before filter
sigMetricsBefore <- grep("BeforeFilter", sigMetrics, value = TRUE)
pcaRes <- prcomp(myData[,sigMetricsBefore])
eigenvalues <- pcaRes$sd^2 # https://stats.stackexchange.com/questions/9500/why-do-the-r-functions-princomp-and-prcomp-give-different-eigenvalues
percExpl <- round(eigenvalues/sum(eigenvalues)*100, 1)
curColor <- sampleTab[rownames(myData), "color"]
curPch <- sampleTab[rownames(myData), "pch"]
pdf(paste0(outfilePrefix, ".beforeFilter.PCA.pdf"), height = 6, width = 6)
#par(mfrow=c(2,1))
plot(pcaRes$x[,1:2], pch=curPch, col=curColor, bg=curColor, cex = 2, main = paste0("PC1: ", percExpl[1], " PC2:", percExpl[2]))
#barplot(eigenvalues, main="PCA Eigenvalues", col="black", xlab = "", ylab = "")
invisible(dev.off())

# after filter
sigMetricsAfter <- grep("AfterFilter", sigMetrics, value = TRUE)
pcaRes <- prcomp(myData[,sigMetricsAfter])
eigenvalues <- pcaRes$sd^2 # https://stats.stackexchange.com/questions/9500/why-do-the-r-functions-princomp-and-prcomp-give-different-eigenvalues
percExpl <- round(eigenvalues/sum(eigenvalues)*100, 1)
curColor <- sampleTab[rownames(myData), "color"]
curPch <- sampleTab[rownames(myData), "pch"]
pdf(paste0(outfilePrefix, ".afterFilter.PCA.pdf"), height = 6, width = 6)
#par(mfrow=c(2,1))
plot(pcaRes$x[,1:2], pch=curPch, col=curColor, bg=curColor, cex = 2, main = paste0("PC1: ", percExpl[1], " PC2:", percExpl[2]))
#barplot(eigenvalues, main="PCA Eigenvalues", col="black", xlab = "", ylab = "")
invisible(dev.off())




