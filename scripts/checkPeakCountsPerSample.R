rm(list=ls())
baseDir <- "/media/mwschmid/myData/MWSchmid/Ivana_ATAC"
scriptDir <- file.path(baseDir, "scripts")
designFile <- file.path(baseDir, "report", "sampleTab.csv")
infileName <- file.path(baseDir, "GitIgnore_results", "peakCountsPerSample.csv")
outfilePrefix <- file.path(baseDir, "report", "peakCountsPerSample")
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
myData$time <- sampleTab[myData$sample,"time"]
myData$lab <- sampleTab[myData$sample,"lab"]
myData$batchWithin <- sampleTab[myData$sample,"batch"]
metrics <- c("numPeaks")

#########################################################################################
# test for differences
matrixCols <- paste0(c("percSS", "pValue"), '_', rep(c("peakType", "time", "lab", "batchWithin"), each = 2))
allTests <- matrix(NA, ncol = length(matrixCols), nrow = length(metrics), dimnames = list(metrics, matrixCols))
for (curMetric in metrics) {
  res <- anova(lm(terms(as.formula(paste0(curMetric, '~ peakType + time + batchWithin + lab')), keep.order = TRUE), data = myData))
  allTerms <- rownames(res)[1:4]
  allTests[curMetric, paste0("percSS_", allTerms)] <- round(res[allTerms, "Sum Sq"]/sum(res[, "Sum Sq"])*100, 2)
  allTests[curMetric, paste0("pValue_", allTerms)] <- round(res[allTerms, "Pr(>F)"], 5)
}
write.csv(allTests, paste0(outfilePrefix, ".anova.csv"))

#########################################################################################
# plots, use predefined colors and pchs
#sampleTab$G <- with(sampleTab, paste(time, lab, sep = '_'))
#allGroups <- unique(sampleTab$G)
#forColor <- brewer.pal(length(allGroups), "Paired"); names(forColor) <- paste0(rep(c("T1_", "T2_"), 5), rep(c("Bern", "Hannover1", "Hannover2", "Munster", "Zurich"), each = 2))
#forPch <- c("T1" = 15, "T2" = 16)
#sampleTab$color <- forColor[sampleTab$G]
#sampleTab$pch <- forPch[sampleTab$time]

#########################################################################################
### PCA
pValCols <- grep("^pValue", colnames(allTests), value = TRUE)
allTests[rowSums(is.na(allTests[,pValCols])) == 3, pValCols] <- 1
sigMetrics <- rownames(allTests)[rowSums(allTests[,pValCols] < 0.01) > 0]

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




