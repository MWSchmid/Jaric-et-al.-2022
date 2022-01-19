#!/usr/bin/env Rscript

rm(list=ls())
baseDir <- "/media/mwschmid/myData/MWSchmid/Ivana_ATAC"
infileName <- file.path(baseDir, "GitIgnore_results", "union.3.BP.bothTimePoints.Rdata")
outfilePrefix <- file.path(baseDir, "report", "forFigureEditing", "featurePlot")
chipAnnoInfile <- "/home/marc/tempIJA/union.3.BP.csAnnot.csv"
rm(list=ls())

## arguments from commandline
myarg <- commandArgs()
argPos <- grep("--args", myarg, fixed = TRUE)
infileName <- myarg[argPos+1]
chipAnnoInfile <- myarg[argPos+2]
outfilePrefix <- myarg[argPos+3]

# rest
suppressPackageStartupMessages({
  options(java.parameters = "-Xmx12g")
  library("RColorBrewer")
  library("DESeq2")
  library("plotrix")
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
### get significant peaks from TP1 and TP2
infileNameTP1 <- gsub("bothTimePoints", "onlyTimePointOne", infileName)
load(infileNameTP1)
anySigT1 <- Reduce(union, lapply(deResults[grep("^ot", names(deResults), value = TRUE)], function(x) x$get_significant_entries(0.01, 0.01, 0.5)$any))
infileNameTP2 <- gsub("bothTimePoints", "onlyTimePointTwo", infileName)
load(infileNameTP2)
anySigT2 <- Reduce(union, lapply(deResults[grep("^ot", names(deResults), value = TRUE)], function(x) x$get_significant_entries(0.01, 0.01, 0.5)$any))

##############################################################################
### load common file
load(infileName)
csAnno <- read.csv(chipAnnoInfile, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
rownames(csAnno) <- paste0("chr", rownames(csAnno))

##############################################################################
### get accessibility and overlaps (copy-paste from openAndClose.R)
f.get.most.accessible <- function(x, percThresh = 0.025) {
  upperBound <- quantile(x, 1-percThresh)
  out <- (x > upperBound)
  return(out)
}

f.get.least.accessible <- function(x, percThresh = 0.025) {
  lowerBound <- quantile(x, percThresh)
  out <- (x < lowerBound)
  return(out)
}

f.get.most.and.least.accessible <- function(x, percThresh = 0.025) {
  lowerBound <- quantile(x, percThresh)
  upperBound <- quantile(x, 1-percThresh)
  out <- (x < lowerBound) | (x > upperBound)
  return(out)
}

extremes <- apply(myNormData, 2, f.get.most.and.least.accessible)
mostAccess <- apply(myNormData, 2, f.get.most.accessible)
leastAccess <- apply(myNormData, 2, f.get.least.accessible)
occurences <- rowSums(extremes)
occMostAccess <- rowSums(mostAccess)
occLeastAccess <- rowSums(leastAccess)

sitesForSpiderGraph <- list(
  all = rownames(myNormData),
  diffT1 = anySigT1,
  diffT2 = anySigT2,
  open = rownames(mostAccess)[occMostAccess > 0],
  close = rownames(leastAccess)[occLeastAccess > 0]
)

##############################################################################
### get the annotation and do the spider plots, requires cs annot
csAnno <- subset(csAnno, !grepl("^prom", simpleAnno))# skip prom 1 and 2, it's almost nothing

labs <- c("Intergenic", "TSS-2kb", "TSS-1kb", "TSS+1kb", "TSS+2kb", "5'UTR", "First exon", "First intron", "Exon", "Intron", "3'UTR", "TES+1kb", "TES+2kb")
names(labs) <- c("distalIntergenic", "USTSS2", "USTSS1", "DSTSS1", "DSTSS2", "5pUTR", "firstExon", "firstIntron", "exon", "intron", "3pUTR", "DS1", "DS2")

##############################################################################
# collect data
numbersForSpiderGraph <- lapply(sitesForSpiderGraph, function(x) table(csAnno[intersect(x, rownames(csAnno)), "simpleAnno"]))
percentForSpiderGraph <- lapply(numbersForSpiderGraph, function(x) round(x/sum(x)*100, 2))
out <- do.call("rbind", percentForSpiderGraph)[,names(labs)]
colnames(out) <- labs
write.csv(t(out), paste0(outfilePrefix, "_percentagesPerFeature.csv"))

#head(csAnno)
#head(sitesForSpiderGraph$all)

##############################################################################
# merge
forPlot <- do.call("rbind", percentForSpiderGraph[c("all", "open", "close")])[,names(labs)]
radialLimit <- pretty(c(0, max(forPlot)))
colForPlot <- c("black", "slateblue2", "firebrick2")
pdf(paste0(outfilePrefix, "_openClose.pdf"), height = 8, width = 8)
radial.plot(forPlot, labels = labs, main = "", radial.lim = radialLimit, 
            start = pi/2, line.col = colForPlot, poly.col = rep(NA, length(colForPlot)),
            rp.type = "p", show.grid.labels = 3, lwd = 3, point.symbols = 16,
            point.col = colForPlot, show.centroid = FALSE)
invisible(dev.off())

forPlot <- do.call("rbind", percentForSpiderGraph[c("all", "diffT1", "diffT2")])[,names(labs)]
radialLimit <- pretty(c(0, max(forPlot)))
colForPlot <- c("black", "darkgoldenrod3", "darkorchid3")
pdf(paste0(outfilePrefix, "_anyDifference.pdf"), height = 8, width = 8)
radial.plot(forPlot, labels = labs, main = "", radial.lim = radialLimit, 
            start = pi/2, line.col = colForPlot, poly.col = rep(NA, length(colForPlot)),
            rp.type = "p", show.grid.labels = 3, lwd = 3, point.symbols = 16,
            point.col = colForPlot, show.centroid = FALSE)
invisible(dev.off())



