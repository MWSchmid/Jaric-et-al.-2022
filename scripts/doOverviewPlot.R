rm(list=ls())
baseDir <- "/media/mwschmid/myData/MWSchmid/Ivana_ATAC"
scriptDir <- file.path(baseDir, "scripts")
designFile <- file.path(baseDir, "report", "sampleTab.csv")
infileName <- file.path("/home/marc/tempIJA", "union.3.BP.counts.txt")
outfilePrefix <- file.path(baseDir, "GitIgnore_results", "union.3.BP.plots")
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
  library("DESeq2")
  library("gtools")
  library("vegan")
  library("pheatmap")
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
temp <- read.table(infileName, sep = '\t', stringsAsFactors = FALSE, header = TRUE, row.names=1)
# replace sample colnames
samplePattern <- "S[[:digit:]]{1,2}\\.deDup\\.bam|T[[:digit:]]{1,2}\\.deDup\\.bam"
sampleMatches <- regexpr(samplePattern, colnames(temp)[6:ncol(temp)])
sampleMatchesStr <- regmatches(colnames(temp)[6:ncol(temp)], sampleMatches)
sampleMatchesStr <- gsub("\\.deDup\\.bam", '', sampleMatchesStr, fixed = FALSE)
colnames(temp)[6:ncol(temp)] <- sampleMatchesStr

# add rownames
rownamesToRealNames <- rownames(temp)
rownamesToRealNames <- cbind(with(temp, paste(Chr, Start, End, sep = '_')), rownamesToRealNames)
# remove duplicate entries
temp <- temp[,c("Chr", "Start", "End", rownames(sampleTab))]
before <- nrow(temp)
temp <- unique(temp)
after <- nrow(temp)
if (before != after) { cat("##### removed", before-after, "non-unique rows\n")}
# some changes for the annotation later on
#temp$Chr <- gsub("^chr", "", temp$Chr)
rownames(temp) <- with(temp, paste(Chr, Start, End, sep = '_'))
myData <- temp[,rownames(sampleTab)]; rm(temp)

#########################################################################################
# design and contrasts
sampleTab$G <- with(sampleTab, paste(time, lab, sep = '_'))
sampleTab$G__ <- factor(sampleTab$G) # don't specify levels, we anyway do the contrasts by hand
formulaString <- "~0+G__"
allTimePoints <- unique(sampleTab$time)
allLabs <- unique(sampleTab$lab)
design <- model.matrix(formula(formulaString), data = sampleTab, contrasts.arg = NULL)

#########################################################################################
# remove entries with very low values (<5 in all four samples) - no filter on variance
myData <- f.strip.data(myData, minVal = 5, minTimes = 3, lowerVarQuantileToRemove = 0, colsToStrip = rownames(sampleTab))
write.csv(myData, paste0(outfilePrefix, ".raw.counts.csv"))

#########################################################################################
# normalize the data with DESeq2
myNormData <- f.normalize.counts.DESeq(myData, sampleTab, formulaString)
f.generic.correlation.matrix.alt.col(myNormData, rDir, paste0(filePrefixNoPath, "_norm_pearson"), "pearson", FALSE)
write.csv(myData, paste0(outfilePrefix, ".normalized.counts.csv"))
#quit("no", 0)

#########################################################################################
# do some overview plot
byTab <- data.frame(sample = rownames(sampleTab), group = sampleTab$G__, stringsAsFactors = FALSE)
myAverageNormData <- f.summarize.columns(myNormData, byTab, mean)

#########################################################################################
# correlation matrices with side cols and values
toPlot <- cor(myNormData)
annoColor <- unique(sampleTab[, c("G", "color")])
forColor <- annoColor$color; names (forColor) <- annoColor$G
annoColorList <- list(Group = forColor)
annoCol <- sampleTab[rownames(toPlot), c("G", "G__")]; annoCol$G__ <- NULL
colnames(annoCol) <- c("Group")
pdf(paste0(outfilePrefix, ".sampleCorrelation.pdf"), width=12, height=10)
grid::grid.draw(pheatmap(toPlot, color = colorRampPalette(brewer.pal(n = 7, name ="YlGnBu"))(100), annotation_col = annoCol, annotation_colors = annoColorList)$gtable)
dev.off()
#rendered <- pheatmap(toPlot, color = colorRampPalette(brewer.pal(n = 7, name ="YlGnBu"))(100), annotation_col = annoCol, annotation_colors = annoColorList)
#save_pheatmap_pdf(rendered, paste0(outfilePrefix, ".sampleCorrelation.pdf"), 12, 10)
forBreaks <- range(pretty(range(toPlot)))
breaksForAllPeaks <- seq(forBreaks[1], forBreaks[2], length.out = 21)

varVec <- apply(myNormData, 1, var)
lowerBound <- quantile(varVec, 0.9)
highVarData <- myNormData[varVec>=lowerBound,]
toPlot <- cor(highVarData)
pdf(paste0(outfilePrefix, ".sampleCorrelation.onlyHighVarPeaks.pdf"), width=12, height=10)
grid::grid.draw(pheatmap(toPlot, color = colorRampPalette(brewer.pal(n = 7, name ="YlGnBu"))(100), annotation_col = annoCol, annotation_colors = annoColorList)$gtable)
dev.off()
#rendered <- pheatmap(toPlot, color = colorRampPalette(brewer.pal(n = 7, name ="YlGnBu"))(100), annotation_col = annoCol, annotation_colors = annoColorList)
#save_pheatmap_pdf(rendered, paste0(outfilePrefix, ".sampleCorrelation.onlyHighVarPeaks.pdf"), 12, 10)
forBreaks <- range(pretty(range(toPlot)))
breaksForVarPeaks <- seq(forBreaks[1], forBreaks[2], length.out = 21)

for (curTimepoint in unique(sampleTab$time)) {
  subSampleTab <- subset(sampleTab, time == curTimepoint)
  
  toPlot <- cor(myNormData[,rownames(subSampleTab)])
  annoColor <- unique(subSampleTab[, c("G", "color")])
  forColor <- annoColor$color; names (forColor) <- annoColor$G
  annoColorList <- list(Group = forColor)
  annoCol <- subSampleTab[rownames(toPlot), c("G", "G__")]; annoCol$G__ <- NULL
  colnames(annoCol) <- c("Group")
  pdf(paste0(outfilePrefix, ".sampleCorrelation.", curTimepoint, ".pdf"), width=12, height=10)
  grid::grid.draw(pheatmap(toPlot, color = colorRampPalette(brewer.pal(n = 7, name ="YlGnBu"))(20), annotation_col = annoCol, annotation_colors = annoColorList, breaks = breaksForAllPeaks)$gtable)
  dev.off()
  
  varVec <- apply(myNormData[,rownames(subSampleTab)], 1, var)
  lowerBound <- quantile(varVec, 0.9)
  highVarDataWithin <- myNormData[varVec>=lowerBound, rownames(subSampleTab)]
  toPlot <- cor(highVarDataWithin)
  pdf(paste0(outfilePrefix, ".sampleCorrelation.", curTimepoint, ".onlyHighVarPeaks.pdf"), width=12, height=10)
  grid::grid.draw(pheatmap(toPlot, color = colorRampPalette(brewer.pal(n = 7, name ="YlGnBu"))(20), annotation_col = annoCol, annotation_colors = annoColorList, breaks = breaksForVarPeaks)$gtable)
  dev.off()
}

quit("no", 0)

#########################################################################################
# and a PCA, use pre-defined colors
#allGroups <- unique(sampleTab$G)
#forColor <- brewer.pal(length(allGroups), "Paired"); names(forColor) <- paste0(rep(c("T1_", "T2_"), 5), rep(c("Bern", "Hannover1", "Hannover2", "Munster", "Zurich"), each = 2))
#forPch <- c("T1" = 15, "T2" = 16)
#sampleTab$color <- forColor[sampleTab$G]
#sampleTab$pch <- forPch[sampleTab$time]
pcaRes <- princomp(myNormData)
eigenvalues <- pcaRes$sd^2 # https://stats.stackexchange.com/questions/9500/why-do-the-r-functions-princomp-and-prcomp-give-different-eigenvalues
percExpl <- round(eigenvalues/sum(eigenvalues)*100, 1)
curColor <- sampleTab[colnames(myNormData), "color"]
curPch <- sampleTab[colnames(myNormData), "pch"]
f.open.figure(rDir, paste0(filePrefixNoPath, "_PCA.svg"), TRUE, height = 10, width = 5)
par(mfrow=c(2,1))
plot(pcaRes$loadings[,1:2], pch=curPch, col=curColor, bg=curColor, cex = 2, main = paste0("PC1: ", percExpl[1], " PC2:", percExpl[2]))
barplot(eigenvalues, main="PCA Eigenvalues", col="black", xlab = "", ylab = "")
f.close.figure()

subGenes <- t(myNormData)
RDA2 <- rda(subGenes ~ G__, sampleTab)
axes <- summary(RDA2)$site
temp <- summary(RDA2)$cont
curColor <- sampleTab[colnames(myNormData), "color"]
percExpl <- round(temp$importance[,1:2]["Proportion Explained",]*100, 2)
f.open.figure(rDir, paste0(filePrefixNoPath, "_RDA.svg"), TRUE, height = 10, width = 5)
par(mfrow=c(2,1))
plot(axes[,1:2], pch=rep(16, length(curColor)), col=curColor, bg=curColor, cex = 2, main = paste0("RDA1: ", percExpl[1], " RDA2:", percExpl[2])) # mono hist is a dot, mix hist is a triangle
barplot(temp$importance["Eigenvalue",], main="PCA Eigenvalues", col="black", xlab = "", ylab = "")
f.close.figure()

pcaRes <- princomp(highVarData)
eigenvalues <- pcaRes$sd^2 # https://stats.stackexchange.com/questions/9500/why-do-the-r-functions-princomp-and-prcomp-give-different-eigenvalues
percExpl <- round(eigenvalues/sum(eigenvalues)*100, 1)
curColor <- sampleTab[colnames(highVarData), "color"]
curPch <- sampleTab[colnames(highVarData), "pch"]
f.open.figure(rDir, paste0(filePrefixNoPath, "_PCA_highVar.svg"), TRUE, height = 10, width = 5)
par(mfrow=c(2,1))
plot(pcaRes$loadings[,1:2], pch=curPch, col=curColor, bg=curColor, cex = 2, main = paste0("PC1: ", percExpl[1], " PC2:", percExpl[2]))
barplot(eigenvalues, main="PCA Eigenvalues", col="black", xlab = "", ylab = "")
f.close.figure()

subGenes <- t(highVarData)
RDA2 <- rda(subGenes ~ G__, sampleTab)
axes <- summary(RDA2)$site
temp <- summary(RDA2)$cont
curColor <- sampleTab[colnames(highVarData), "color"]
percExpl <- round(temp$importance[,1:2]["Proportion Explained",]*100, 2)
f.open.figure(rDir, paste0(filePrefixNoPath, "_RDA_highVar.svg"), TRUE, height = 10, width = 5)
par(mfrow=c(2,1))
plot(axes[,1:2], pch=rep(16, length(curColor)), col=curColor, bg=curColor, cex = 2, main = paste0("RDA1: ", percExpl[1], " RDA2:", percExpl[2])) # mono hist is a dot, mix hist is a triangle
barplot(temp$importance["Eigenvalue",], main="PCA Eigenvalues", col="black", xlab = "", ylab = "")
f.close.figure()

#########################################################################################
# PCA/RDA with time point 1
subSampleTab <- subset(sampleTab, time == "T1")
subNormData <- myNormData[,rownames(subSampleTab)]
pcaRes <- princomp(subNormData)
eigenvalues <- pcaRes$sd^2 # https://stats.stackexchange.com/questions/9500/why-do-the-r-functions-princomp-and-prcomp-give-different-eigenvalues
percExpl <- round(eigenvalues/sum(eigenvalues)*100, 1)
curColor <- sampleTab[colnames(subNormData), "color"]
curPch <- sampleTab[colnames(subNormData), "pch"]
f.open.figure(rDir, paste0(filePrefixNoPath, "_onlyTimePointOne_PCA.svg"), TRUE, height = 10, width = 5)
par(mfrow=c(2,1))
plot(pcaRes$loadings[,1:2], pch=curPch, col=curColor, bg=curColor, cex = 2, main = paste0("PC1: ", percExpl[1], " PC2:", percExpl[2]))
barplot(eigenvalues, main="PCA Eigenvalues", col="black", xlab = "", ylab = "")
f.close.figure()

subGenes <- t(subNormData)
RDA2 <- rda(subGenes ~ G__, subSampleTab)
axes <- summary(RDA2)$site
temp <- summary(RDA2)$cont
curColor <- sampleTab[colnames(subNormData), "color"]
percExpl <- round(temp$importance[,1:2]["Proportion Explained",]*100, 2)
f.open.figure(rDir, paste0(filePrefixNoPath, "_onlyTimePointOne_RDA.svg"), TRUE, height = 10, width = 5)
par(mfrow=c(2,1))
plot(axes[,1:2], pch=rep(16, length(curColor)), col=curColor, bg=curColor, cex = 2, main = paste0("RDA1: ", percExpl[1], " RDA2:", percExpl[2])) # mono hist is a dot, mix hist is a triangle
barplot(temp$importance["Eigenvalue",], main="PCA Eigenvalues", col="black", xlab = "", ylab = "")
f.close.figure()

varVec <- apply(subNormData, 1, var)
lowerBound <- quantile(varVec, 0.9)
highVarData <- subNormData[varVec>=lowerBound,]

pcaRes <- princomp(highVarData)
eigenvalues <- pcaRes$sd^2 # https://stats.stackexchange.com/questions/9500/why-do-the-r-functions-princomp-and-prcomp-give-different-eigenvalues
percExpl <- round(eigenvalues/sum(eigenvalues)*100, 1)
curColor <- sampleTab[colnames(highVarData), "color"]
curPch <- sampleTab[colnames(highVarData), "pch"]
f.open.figure(rDir, paste0(filePrefixNoPath, "_onlyTimePointOne_PCA_highVar.svg"), TRUE, height = 10, width = 5)
par(mfrow=c(2,1))
plot(pcaRes$loadings[,1:2], pch=curPch, col=curColor, bg=curColor, cex = 2, main = paste0("PC1: ", percExpl[1], " PC2:", percExpl[2]))
barplot(eigenvalues, main="PCA Eigenvalues", col="black", xlab = "", ylab = "")
f.close.figure()

subGenes <- t(highVarData)
RDA2 <- rda(subGenes ~ G__, subSampleTab)
axes <- summary(RDA2)$site
temp <- summary(RDA2)$cont
curColor <- sampleTab[colnames(highVarData), "color"]
percExpl <- round(temp$importance[,1:2]["Proportion Explained",]*100, 2)
f.open.figure(rDir, paste0(filePrefixNoPath, "_onlyTimePointOne_RDA_highVar.svg"), TRUE, height = 10, width = 5)
par(mfrow=c(2,1))
plot(axes[,1:2], pch=rep(16, length(curColor)), col=curColor, bg=curColor, cex = 2, main = paste0("RDA1: ", percExpl[1], " RDA2:", percExpl[2])) # mono hist is a dot, mix hist is a triangle
barplot(temp$importance["Eigenvalue",], main="PCA Eigenvalues", col="black", xlab = "", ylab = "")
f.close.figure()

#########################################################################################
# PCA/RDA with time point 2
subSampleTab <- subset(sampleTab, time == "T2")
subNormData <- myNormData[,rownames(subSampleTab)]
pcaRes <- princomp(subNormData)
eigenvalues <- pcaRes$sd^2 # https://stats.stackexchange.com/questions/9500/why-do-the-r-functions-princomp-and-prcomp-give-different-eigenvalues
percExpl <- round(eigenvalues/sum(eigenvalues)*100, 1)
curColor <- sampleTab[colnames(subNormData), "color"]
curPch <- sampleTab[colnames(subNormData), "pch"]
f.open.figure(rDir, paste0(filePrefixNoPath, "_onlyTimePointTwo_PCA.svg"), TRUE, height = 10, width = 5)
par(mfrow=c(2,1))
plot(pcaRes$loadings[,1:2], pch=curPch, col=curColor, bg=curColor, cex = 2, main = paste0("PC1: ", percExpl[1], " PC2:", percExpl[2]))
barplot(eigenvalues, main="PCA Eigenvalues", col="black", xlab = "", ylab = "")
f.close.figure()

subGenes <- t(subNormData)
RDA2 <- rda(subGenes ~ G__, subSampleTab)
axes <- summary(RDA2)$site
temp <- summary(RDA2)$cont
curColor <- sampleTab[colnames(subNormData), "color"]
percExpl <- round(temp$importance[,1:2]["Proportion Explained",]*100, 2)
f.open.figure(rDir, paste0(filePrefixNoPath, "_onlyTimePointTwo_RDA.svg"), TRUE, height = 10, width = 5)
par(mfrow=c(2,1))
plot(axes[,1:2], pch=rep(16, length(curColor)), col=curColor, bg=curColor, cex = 2, main = paste0("RDA1: ", percExpl[1], " RDA2:", percExpl[2])) # mono hist is a dot, mix hist is a triangle
barplot(temp$importance["Eigenvalue",], main="PCA Eigenvalues", col="black", xlab = "", ylab = "")
f.close.figure()

varVec <- apply(subNormData, 1, var)
lowerBound <- quantile(varVec, 0.9)
highVarData <- subNormData[varVec>=lowerBound,]

pcaRes <- princomp(highVarData)
eigenvalues <- pcaRes$sd^2 # https://stats.stackexchange.com/questions/9500/why-do-the-r-functions-princomp-and-prcomp-give-different-eigenvalues
percExpl <- round(eigenvalues/sum(eigenvalues)*100, 1)
curColor <- sampleTab[colnames(highVarData), "color"]
curPch <- sampleTab[colnames(highVarData), "pch"]
f.open.figure(rDir, paste0(filePrefixNoPath, "_onlyTimePointTwo_PCA_highVar.svg"), TRUE, height = 10, width = 5)
par(mfrow=c(2,1))
plot(pcaRes$loadings[,1:2], pch=curPch, col=curColor, bg=curColor, cex = 2, main = paste0("PC1: ", percExpl[1], " PC2:", percExpl[2]))
barplot(eigenvalues, main="PCA Eigenvalues", col="black", xlab = "", ylab = "")
f.close.figure()

subGenes <- t(highVarData)
RDA2 <- rda(subGenes ~ G__, subSampleTab)
axes <- summary(RDA2)$site
temp <- summary(RDA2)$cont
curColor <- sampleTab[colnames(highVarData), "color"]
percExpl <- round(temp$importance[,1:2]["Proportion Explained",]*100, 2)
f.open.figure(rDir, paste0(filePrefixNoPath, "_onlyTimePointTwo_RDA_highVar.svg"), TRUE, height = 10, width = 5)
par(mfrow=c(2,1))
plot(axes[,1:2], pch=rep(16, length(curColor)), col=curColor, bg=curColor, cex = 2, main = paste0("RDA1: ", percExpl[1], " RDA2:", percExpl[2])) # mono hist is a dot, mix hist is a triangle
barplot(temp$importance["Eigenvalue",], main="PCA Eigenvalues", col="black", xlab = "", ylab = "")
f.close.figure()

#########################################################################################
# legend colors
forLegend <- unique(sampleTab[,c("color", "pch", "G")])
f.open.figure(rDir, paste0(filePrefixNoPath, "_legendColorAndPch.svg"), TRUE, height = 10, width = 5)
f.plot.legend(forLegend$color, forLegend$pch, forLegend$G)
f.close.figure()

#########################################################################################
# do scatters at the end, may take long
f.do.some.overview.alt.col(myAverageNormData, rDir, paste0(filePrefixNoPath, "_average_norm"), skipScatters = FALSE)

