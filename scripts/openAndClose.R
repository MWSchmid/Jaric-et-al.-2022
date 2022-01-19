#!/usr/bin/env Rscript

rm(list=ls())
baseDir <- "/media/mwschmid/myData/MWSchmid/Ivana_ATAC"
infileName <- file.path(baseDir, "GitIgnore_results", "union.3.BP.bothTimePoints.Rdata")
chipAnnoInfile <- "/home/marc/tempIJA/union.3.BP.csAnnot.csv"
homeAnnoInfile <- "/home/marc/tempIJA/union.3.BP.homer.csv"
outfilePrefix <- file.path(baseDir, "GitIgnore_results", "openAndClose", "union.3.BP.bothTimePoints")
rm(list=ls())

## arguments from commandline
myarg <- commandArgs()
argPos <- grep("--args", myarg, fixed = TRUE)
infileName <- myarg[argPos+1]
chipAnnoInfile <- myarg[argPos+2]
homeAnnoInfile <- myarg[argPos+3]
outfilePrefix <- myarg[argPos+4]

# rest
suppressPackageStartupMessages({
  options(java.parameters = "-Xmx12g")
  library("gplots")
  library("vegan")
  library("Rtsne")
  #source("/media/mwschmid/myData/MWSchmid/Development/R/lm_wrapper.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapperAds.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper_overviewBandwidthFix.R")
})

#ensembl <- f.connect.to.ensembl(martDataSet = "mmusculus_gene_ensembl", martName = "ENSEMBL_MART_ENSEMBL")

rDir <- dirname(outfilePrefix)
outfilePrefixNoPath <- basename(outfilePrefix)

##############################################################################
### load file
load(infileName) # sampleTab, deResults, (deResultsNoBatch), myData, myNormData, (design, myCont)
csAnno <- read.csv(chipAnnoInfile, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
hoAnno <- read.csv(homeAnnoInfile, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
hoAnno$Chr <- gsub("^chr", "", hoAnno$Chr)
rownames(hoAnno) <- paste0(hoAnno$Chr, '_', hoAnno$Start-1, '_', hoAnno$End)

##############################################################################
### get accessibility and overlaps
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
pdf(paste0(outfilePrefix, "_occurences.pdf"), height = 5, width = 15)
par(mfrow = c(1, 3))
hist(occurences[occurences>0], main = "Both extremes", xlab = "Number of samples", ylab = "Number of sites")
hist(occMostAccess[occMostAccess>0], main = "Most accessible", xlab = "Number of samples", ylab = "Number of sites", col = "slateblue2")
hist(occLeastAccess[occLeastAccess>0], main = "Least accessible", xlab = "Number of samples", ylab = "Number of sites", col = "firebrick2")
invisible(dev.off())

##############################################################################
### get accessibility and overlaps, gene-based
hoAnnoOnlyTSS <- subset(hoAnno, abs(DistanceToTSS) < 2000 & nchar(GeneID) > 3)
rownames(hoAnnoOnlyTSS) <- paste0("chr", rownames(hoAnnoOnlyTSS))
extremeInGenes <- intersect(rownames(extremes), rownames(hoAnnoOnlyTSS))
extremesGenes <- aggregate(extremes[extremeInGenes, ], by = list(gene = hoAnnoOnlyTSS[extremeInGenes, "GeneID"]), function(x) sum(x) > 0)
rownames(extremesGenes) <- extremesGenes$gene
extremesGenes <- extremesGenes[,colnames(extremes)]

mostAccessInGenes <- intersect(rownames(mostAccess), rownames(hoAnnoOnlyTSS))
mostAccessGenes <- aggregate(mostAccess[mostAccessInGenes, ], by = list(gene = hoAnnoOnlyTSS[mostAccessInGenes, "GeneID"]), function(x) sum(x) > 0)
rownames(mostAccessGenes) <- mostAccessGenes$gene
mostAccessGenes <- mostAccessGenes[,colnames(mostAccess)]

leastAccessInGenes <- intersect(rownames(leastAccess), rownames(hoAnnoOnlyTSS))
leastAccessGenes <- aggregate(leastAccess[leastAccessInGenes, ], by = list(gene = hoAnnoOnlyTSS[leastAccessInGenes, "GeneID"]), function(x) sum(x) > 0)
rownames(leastAccessGenes) <- leastAccessGenes$gene
leastAccessGenes <- leastAccessGenes[,colnames(leastAccess)]

occurencesGenes <- rowSums(extremesGenes)
occMostAccessGenes <- rowSums(mostAccessGenes)
occLeastAccessGenes <- rowSums(leastAccessGenes)
pdf(paste0(outfilePrefix, "_occurencesGenes.pdf"), height = 5, width = 15)
par(mfrow = c(1, 3))
hist(occurencesGenes[occurencesGenes>0], main = "Both extremes", xlab = "Number of samples", ylab = "Number of genes")
hist(occMostAccessGenes[occMostAccessGenes>0], main = "Most accessible", xlab = "Number of samples", ylab = "Number of genes", col = "slateblue2")
hist(occLeastAccessGenes[occLeastAccessGenes>0], main = "Least accessible", xlab = "Number of samples", ylab = "Number of genes", col = "firebrick2")
invisible(dev.off())

##############################################################################
### Test distances with PERMANOVA and visualize distances with t-sna
sampleTab_T1 <- subset(sampleTab, time == "T1")
sampleTab_T2 <- subset(sampleTab, time == "T2")
allTests <- list()
allPlots <- list()

sitesForTestAndPlot <- list(
  all = rownames(myNormData),
  extreme = rownames(extremes)[occurences > 0],
  open = rownames(mostAccess)[occMostAccess > 0],
  close = rownames(leastAccess)[occLeastAccess > 0]
)

for (curSet in names(sitesForTestAndPlot)) {
  sitesForPlot <- sitesForTestAndPlot[[curSet]]
  forTest <- t(myNormData[sitesForPlot,rownames(sampleTab)])
  forTest_T1 <- t(myNormData[sitesForPlot,rownames(sampleTab_T1)])
  forTest_T2 <- t(myNormData[sitesForPlot,rownames(sampleTab_T2)])
  allTests[[curSet]] <- adonis(terms(forTest ~ time + lab, keep.order = TRUE), data = sampleTab, method = "manhattan", permutations = 9999)
  allTests[[paste0(curSet, "_T1_with_batch")]] <- adonis(terms(forTest_T1 ~ batch + lab, keep.order = TRUE), data = sampleTab_T1, method = "manhattan", permutations = 9999)
  allTests[[paste0(curSet, "_T2_with_batch")]] <- adonis(terms(forTest_T2 ~ batch + lab, keep.order = TRUE), data = sampleTab_T2, method = "manhattan", permutations = 9999)
  allPlots[[curSet]] <- dist(forTest, method = "manhattan")
  allPlots[[paste0(curSet, "_T1")]] <- dist(forTest_T1, method = "manhattan")
  allPlots[[paste0(curSet, "_T2")]] <- dist(forTest_T2, method = "manhattan")
}

##############################################################################
### Write test results
allSets <- names(allTests)
#allEntries <- paste0(rep(c("time", "lab", "batch"), each = 3), '_', rep(c("R2", "percSS", "P")))
allEntries <- paste0(rep(c("time", "lab", "batch"), each = 2), '_', rep(c("percSS", "P")))
testResultsCombined <- matrix(NA, nrow = length(allSets), ncol = length(allEntries), dimnames = list(allSets, allEntries))
for (curSet in names(allTests)) {
  #curSet <- names(allTests)[1]
  mod <- as.matrix(allTests[[curSet]]$aov.tab)
  percSS <- round(mod[1:(nrow(mod)-2),"SumsOfSqs"]/mod["Total","SumsOfSqs"]*100, 2)
  temp <- mod[1:(nrow(mod)-2),]
  testResultsCombined[curSet, paste0(names(percSS), "_percSS")] <- percSS
  #testResultsCombined[curSet, paste0(rownames(temp), "_R2")] <- temp[,"R2"]
  testResultsCombined[curSet, paste0(rownames(temp), "_P")] <- temp[,"Pr(>F)"]
}
write.csv(testResultsCombined, paste0(outfilePrefix, ".permanovas.csv"))

##############################################################################
### Visualize distances with t-sne
for (curSet in names(sitesForTestAndPlot)) {
  pdf(paste0(outfilePrefix, ".dist.plot.tsne.", curSet, ".pdf"), height = 5, width = 15)
  par(mfrow = c(1,3))
  for (curPlot in grep(curSet, names(allPlots), value = TRUE)) {
    distVis <- Rtsne(allPlots[[curPlot]], is_distance = TRUE, perplexity = 7)
    plot(distVis$Y, asp = 1, pch = 16, col = sampleTab[colnames(as.matrix(allPlots[[curPlot]])), "color"], cex = 2, cex.axis = 1.25, cex.lab = 1.25, cex.main = 1.5, xlab = "t-SNE dimension 1", ylab = "t-SNE dimension 2", main = curPlot, las = 1)
  }
  invisible(dev.off())
}

############################################################################################################################################################
### the same as above but without the samples that are far off
toRemove <- c("T14", "T15", "S22")
sampleTabNoOutlier <- sampleTab[setdiff(rownames(sampleTab), toRemove), ]

##############################################################################
### Test distances with PERMANOVA and visualize distances with t-sna
sampleTabNoOutlier_T1 <- subset(sampleTabNoOutlier, time == "T1")
sampleTabNoOutlier_T2 <- subset(sampleTabNoOutlier, time == "T2")
allTests <- list()
allPlots <- list()

sitesForTestAndPlot <- list(
  all = rownames(myNormData),
  extreme = rownames(extremes)[occurences > 0],
  open = rownames(mostAccess)[occMostAccess > 0],
  close = rownames(leastAccess)[occLeastAccess > 0]
)

for (curSet in names(sitesForTestAndPlot)) {
  sitesForPlot <- sitesForTestAndPlot[[curSet]]
  forTest <- t(myNormData[sitesForPlot,rownames(sampleTabNoOutlier)])
  forTest_T1 <- t(myNormData[sitesForPlot,rownames(sampleTabNoOutlier_T1)])
  forTest_T2 <- t(myNormData[sitesForPlot,rownames(sampleTabNoOutlier_T2)])
  allTests[[curSet]] <- adonis(terms(forTest ~ time + lab, keep.order = TRUE), data = sampleTabNoOutlier, method = "manhattan", permutations = 9999)
  allTests[[paste0(curSet, "_T1_with_batch")]] <- adonis(terms(forTest_T1 ~ batch + lab, keep.order = TRUE), data = sampleTabNoOutlier_T1, method = "manhattan", permutations = 9999)
  allTests[[paste0(curSet, "_T2_with_batch")]] <- adonis(terms(forTest_T2 ~ batch + lab, keep.order = TRUE), data = sampleTabNoOutlier_T2, method = "manhattan", permutations = 9999)
  allPlots[[curSet]] <- dist(forTest, method = "manhattan")
  allPlots[[paste0(curSet, "_T1")]] <- dist(forTest_T1, method = "manhattan")
  allPlots[[paste0(curSet, "_T2")]] <- dist(forTest_T2, method = "manhattan")
}

##############################################################################
### Write test results
allSets <- names(allTests)
#allEntries <- paste0(rep(c("time", "lab", "batch"), each = 3), '_', rep(c("R2", "percSS", "P")))
allEntries <- paste0(rep(c("time", "lab", "batch"), each = 2), '_', rep(c("percSS", "P")))
testResultsCombined <- matrix(NA, nrow = length(allSets), ncol = length(allEntries), dimnames = list(allSets, allEntries))
for (curSet in names(allTests)) {
  #curSet <- names(allTests)[1]
  mod <- as.matrix(allTests[[curSet]]$aov.tab)
  percSS <- round(mod[1:(nrow(mod)-2),"SumsOfSqs"]/mod["Total","SumsOfSqs"]*100, 2)
  temp <- mod[1:(nrow(mod)-2),]
  testResultsCombined[curSet, paste0(names(percSS), "_percSS")] <- percSS
  #testResultsCombined[curSet, paste0(rownames(temp), "_R2")] <- temp[,"R2"]
  testResultsCombined[curSet, paste0(rownames(temp), "_P")] <- temp[,"Pr(>F)"]
}
write.csv(testResultsCombined, paste0(outfilePrefix, ".permanovas.extremeSamplesExcluded.csv"))

##############################################################################
### Visualize distances with t-sne
for (curSet in names(sitesForTestAndPlot)) {
  pdf(paste0(outfilePrefix, ".dist.plot.tsne.extremeSamplesExcluded.", curSet, ".pdf"), height = 5, width = 15)
  par(mfrow = c(1,3))
  for (curPlot in grep(curSet, names(allPlots), value = TRUE)) {
    distVis <- Rtsne(allPlots[[curPlot]], is_distance = TRUE, perplexity = 7)
    plot(distVis$Y, asp = 1, pch = 16, col = sampleTabNoOutlier[colnames(as.matrix(allPlots[[curPlot]])), "color"], cex = 2, cex.axis = 1.25, cex.lab = 1.25, cex.main = 1.5, xlab = "t-SNE dimension 1", ylab = "t-SNE dimension 2", main = curPlot, las = 1)
  }
  invisible(dev.off())
}

##############################################################################
### Plot a heatmap, order samples according to time point and lab
labToID <- paste0("lab", 1:5); names(labToID) <- c("Hannover1", "Hannover2", "Bern", "Munster", "Zurich")
sampleTab$labCrypt <- labToID[sampleTab$lab]
sampleTab <- sampleTab[with(sampleTab, order(time, labCrypt)),]
sitesForPlot <- rownames(extremes)[occurences > 0]
toPlot <- myNormData[sitesForPlot, rownames(sampleTab)]

png(paste0(outfilePrefix, ".absVal.png"), height = 2400, width = 1800)
heatmap.2(toPlot, col = f.blackblueyellowredpinkNICE(51), trace="none", scale = "none", margins = c(5,5), dendrogram = "row", Colv = FALSE)
invisible(dev.off())
png(paste0(outfilePrefix, ".scaledVal.png"), height = 2400, width = 1800)
heatmap.2(toPlot, col = f.blueblackyellow(51), trace="none", scale = "row", margins = c(5,5), dendrogram = "row", Colv = FALSE)
invisible(dev.off())
toPlotLimited <- toPlot
toPlotLimited <- t(scale(t(toPlotLimited)))
toPlotLimited[toPlotLimited < (-3)] <- (-3)
toPlotLimited[toPlotLimited > 3] <- 3
png(paste0(outfilePrefix, ".scaledVal.limited.png"), height = 2400, width = 1800)
heatmap.2(toPlotLimited, col = f.blueblackyellow(51), trace="none", scale = "none", margins = c(5,5), dendrogram = "row", Colv = FALSE)
invisible(dev.off())

# the size was a bit trial and error to get rid of the white lines
png(paste0(outfilePrefix, ".scaledVal.limited.open.png"), height = 200+2400*length(sitesForTestAndPlot$open)/length(sitesForTestAndPlot$extreme), width = 1800)
heatmap.2(toPlotLimited[sitesForTestAndPlot$open,], col = f.blueblackyellow(51), trace="none", scale = "none", margins = c(5,5), dendrogram = "row", Colv = FALSE)
invisible(dev.off())

png(paste0(outfilePrefix, ".scaledVal.limited.close.png"), height = 200+2200*length(sitesForTestAndPlot$close)/length(sitesForTestAndPlot$extreme), width = 1800)
heatmap.2(toPlotLimited[sitesForTestAndPlot$close,], col = f.blueblackyellow(51), trace="none", scale = "none", margins = c(5,5), dendrogram = "row", Colv = FALSE)
invisible(dev.off())

# figure sizes
#5979/(5979+18460)*110
#18460/(5979+18460)*110

byTab <- data.frame(sample = rownames(sampleTab), group = sampleTab$G__, stringsAsFactors = FALSE)
toPlot <- f.summarize.columns(toPlot, byTab, mean)
png(paste0(outfilePrefix, ".averaged.absVal.png"), height = 2400, width = 1800)
heatmap.2(toPlot, col = f.blackblueyellowredpinkNICE(51), trace="none", scale = "none", margins = c(5,5), dendrogram = "row", Colv = FALSE)
invisible(dev.off())
png(paste0(outfilePrefix, ".averaged.scaledVal.png"), height = 2400, width = 1800)
heatmap.2(toPlot, col = f.blueblackyellow(51), trace="none", scale = "row", margins = c(5,5), dendrogram = "row", Colv = FALSE)
invisible(dev.off())

######################################################################
### annotation summary
featureOrder <- c("USTSS2", "USTSS1", "DSTSS1", "DSTSS2", "5pUTR", "firstExon", "firstIntron", "exon", "intron", "3pUTR", "DS1", "DS2", "distalIntergenic")
annotSummary <- list(
  mostOrLeast = table(csAnno[gsub("chr", "", rownames(extremes)[occurences > 0]), "simpleAnno"])[featureOrder],
  most = table(csAnno[gsub("chr", "", rownames(mostAccess)[occMostAccess > 0]), "simpleAnno"])[featureOrder],
  least = table(csAnno[gsub("chr", "", rownames(leastAccess)[occLeastAccess > 0]), "simpleAnno"])[featureOrder],
  baseline = table(csAnno$simpleAnno)[featureOrder]
)
annotSummary <- do.call("cbind", annotSummary)
annotSummaryPerc <- t(round(t(annotSummary)/colSums(annotSummary)*100, 2))
write.csv(annotSummary, paste0(outfilePrefix, ".annotSummaryNumber.csv"))
write.csv(annotSummaryPerc, paste0(outfilePrefix, ".annotSummaryPercentage.csv"))

##############################################################################
### Do the same, but averaged
myNormDataAveraged <- f.summarize.columns(myNormData[,byTab$sample], byTab, mean)
extremes <- apply(myNormDataAveraged, 2, f.get.most.and.least.accessible)
mostAccess <- apply(myNormDataAveraged, 2, f.get.most.accessible)
leastAccess <- apply(myNormDataAveraged, 2, f.get.least.accessible)
occurences <- rowSums(extremes)
occMostAccess <- rowSums(mostAccess)
occLeastAccess <- rowSums(leastAccess)
pdf(paste0(outfilePrefix, ".zzz.perGroup_occurences.pdf"), height = 5, width = 15)
par(mfrow = c(1, 3))
hist(occurences[occurences>0], main = "Both extremes", xlab = "Number of samples", ylab = "Number of sites")
hist(occMostAccess[occMostAccess>0], main = "Most accessible", xlab = "Number of samples", ylab = "Number of sites", col = "slateblue2")
hist(occLeastAccess[occLeastAccess>0], main = "Least accessible", xlab = "Number of samples", ylab = "Number of sites", col = "firebrick2")
invisible(dev.off())

##############################################################################
### Plot a heatmap, order samples according to time point and lab
sitesForPlot <- rownames(extremes)[occurences > 0]
toPlot <- myNormData[sitesForPlot, rownames(sampleTab)]

png(paste0(outfilePrefix, ".zzz.perGroup.absVal.png"), height = 2400, width = 1800)
heatmap.2(toPlot, col = f.blackblueyellowredpinkNICE(51), trace="none", scale = "none", margins = c(5,5), dendrogram = "row", Colv = FALSE)
invisible(dev.off())
png(paste0(outfilePrefix, ".zzz.perGroup.scaledVal.png"), height = 2400, width = 1800)
heatmap.2(toPlot, col = f.blueblackyellow(51), trace="none", scale = "row", margins = c(5,5), dendrogram = "row", Colv = FALSE)
invisible(dev.off())

byTab <- data.frame(sample = rownames(sampleTab), group = sampleTab$G__, stringsAsFactors = FALSE)
toPlot <- f.summarize.columns(toPlot, byTab, mean)
png(paste0(outfilePrefix, ".zzz.perGroup.averaged.absVal.png"), height = 2400, width = 1800)
heatmap.2(toPlot, col = f.blackblueyellowredpinkNICE(51), trace="none", scale = "none", margins = c(5,5), dendrogram = "row", Colv = FALSE)
invisible(dev.off())
png(paste0(outfilePrefix, ".zzz.perGroup.averaged.scaledVal.png"), height = 2400, width = 1800)
heatmap.2(toPlot, col = f.blueblackyellow(51), trace="none", scale = "row", margins = c(5,5), dendrogram = "row", Colv = FALSE)
invisible(dev.off())

######################################################################
### annotation summary
featureOrder <- c("USTSS2", "USTSS1", "DSTSS1", "DSTSS2", "5pUTR", "firstExon", "firstIntron", "exon", "intron", "3pUTR", "DS1", "DS2", "distalIntergenic")
annotSummary <- list(
  mostOrLeast = table(csAnno[gsub("chr", "", rownames(extremes)[occurences > 0]), "simpleAnno"])[featureOrder],
  most = table(csAnno[gsub("chr", "", rownames(mostAccess)[occMostAccess > 0]), "simpleAnno"])[featureOrder],
  least = table(csAnno[gsub("chr", "", rownames(leastAccess)[occLeastAccess > 0]), "simpleAnno"])[featureOrder],
  baseline = table(csAnno$simpleAnno)[featureOrder]
)
annotSummary <- do.call("cbind", annotSummary)
annotSummaryPerc <- t(round(t(annotSummary)/colSums(annotSummary)*100, 2))
write.csv(annotSummary, paste0(outfilePrefix, ".zzz.perGroup.annotSummaryNumber.csv"))
write.csv(annotSummaryPerc, paste0(outfilePrefix, ".zzz.perGroup.annotSummaryPercentage.csv"))








