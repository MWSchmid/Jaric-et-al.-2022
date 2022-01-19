#!/usr/bin/env Rscript

rm(list=ls())
baseDir <- "/media/mwschmid/myData/MWSchmid/Ivana_ATAC"
infileName <- file.path(baseDir, "GitIgnore_results", "union.3.BP.bothTimePoints.Rdata")
outfilePrefix <- file.path(baseDir, "GitIgnore_results", "forCircos")
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

binning <- 1e6
allLabs <- paste0("LAB", 1:5)

##############################################################################
### TP1
infileNameTP1 <- gsub("bothTimePoints", "onlyTimePointOne", infileName)
load(infileNameTP1)
avePerLab <- f.summarize.columns(myNormData, data.frame(sample = rownames(sampleTab), group = sapply(sampleTab$lab, f.replace.labels.with.lab)), mean)[,paste0("LAB", 1:5)]
refColors <- unique(sampleTab[,c("lab", "color")])
labColors <- apply(col2rgb(refColors$color), 2, function(x) paste0(x, collapse = ','))
cat(paste0(paste0("color = ", labColors), collapse = '\n'), '\n')

chromStartEnd <- as.data.frame(do.call("rbind", strsplit(rownames(avePerLab), "_")), stringsAsFactors = FALSE)
colnames(chromStartEnd) <- c("chrom", "start", "end")
chromStartEnd$start <- as.numeric(chromStartEnd$start)
chromStartEnd$end <- as.numeric(chromStartEnd$end)
chromStartEnd$bin <- round(chromStartEnd$end/binning)
binnedData <- aggregate(avePerLab, by = list(chrom = chromStartEnd$chrom, bin = chromStartEnd$bin), mean)
binnedData <- binnedData[order(binnedData$chrom, binnedData$bin), ]
binnedData$chrom <- gsub("chr", "mm", binnedData$chrom)
maskThis <- which.max(binnedData$LAB1)
binnedData[maskThis, allLabs] <- colMeans(binnedData[c(maskThis-1, maskThis+1), allLabs])
for (curLab in allLabs) {
  temp <- cbind(binnedData[,c("chrom")], format(binnedData[,c("bin")]*binning, scientific = FALSE), format(binnedData[,c("bin")]*binning+binning, scientific = FALSE), binnedData[,curLab])
  colnames(temp) <- c("#chr", "start", "end", "value")
  write.table(temp, paste0(outfilePrefix, "_TP1_", curLab, ".hist"), row.names = FALSE, quote = FALSE, sep = ' ')
}

##############################################################################
### TP2
infileNameTP2 <- gsub("bothTimePoints", "onlyTimePointTwo", infileName)
load(infileNameTP2)
avePerLab <- f.summarize.columns(myNormData, data.frame(sample = rownames(sampleTab), group = sapply(sampleTab$lab, f.replace.labels.with.lab)), mean)[,paste0("LAB", 1:5)]
refColors <- unique(sampleTab[,c("lab", "color")])
labColors <- apply(col2rgb(refColors$color), 2, function(x) paste0(x, collapse = ','))
cat(paste0(paste0("color = ", labColors), collapse = '\n'), '\n')

chromStartEnd <- as.data.frame(do.call("rbind", strsplit(rownames(avePerLab), "_")), stringsAsFactors = FALSE)
colnames(chromStartEnd) <- c("chrom", "start", "end")
chromStartEnd$start <- as.numeric(chromStartEnd$start)
chromStartEnd$end <- as.numeric(chromStartEnd$end)
chromStartEnd$bin <- round(chromStartEnd$end/binning)
binnedData <- aggregate(avePerLab, by = list(chrom = chromStartEnd$chrom, bin = chromStartEnd$bin), mean)
binnedData <- binnedData[order(binnedData$chrom, binnedData$bin), ]
binnedData$chrom <- gsub("chr", "mm", binnedData$chrom)
maskThis <- which.max(binnedData$LAB1)
binnedData[maskThis, allLabs] <- colMeans(binnedData[c(maskThis-1, maskThis+1), allLabs])
for (curLab in allLabs) {
  temp <- cbind(binnedData[,c("chrom")], format(binnedData[,c("bin")]*binning, scientific = FALSE), format(binnedData[,c("bin")]*binning+binning, scientific = FALSE), binnedData[,curLab])
  colnames(temp) <- c("#chr", "start", "end", "value")
  write.table(temp, paste0(outfilePrefix, "_TP2_", curLab, ".hist"), row.names = FALSE, quote = FALSE, sep = ' ')
}


