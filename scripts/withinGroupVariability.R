#!/usr/bin/env Rscript

rm(list=ls())
baseDir <- "/media/mwschmid/myData/MWSchmid/Ivana_ATAC"
infileName <- file.path(baseDir, "GitIgnore_results", "union.5.BP.bothTimePoints.Rdata")
outfilePrefix <- file.path(baseDir, "GitIgnore_results", "union.5.BP.withinGroupVariability")
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
load(infileName) # sampleTab, deResults, myData, myNormData, design, myCont (sampleTab was added later on)
sampleTab <- read.csv("/media/mwschmid/myData/MWSchmid/Ivana_ATAC/report/sampleTab.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
sampleTab <- sampleTab[colnames(myNormData),]
sampleTab$group <- paste0(sampleTab$time, '_', sampleTab$lab)
allGroups <- unique(sampleTab$group)

##############################################################################
### use manhattan distances
distMat <- as.matrix(dist(t(myNormData), "manhattan"))
diag(distMat) <- NA
distMat[lower.tri(distMat)] <- NA
#curGroup <- allGroups[1]
forDistanceTest <- list()
for (curGroup in allGroups) {
  curSamples <- rownames(sampleTab)[sampleTab$group == curGroup]
  temp <- as.vector(distMat[curSamples, curSamples])
  temp <- temp[!is.na(temp)]
  forDistanceTest[[curGroup]] <- cbind(temp, rep(curGroup, length(temp)))
}
forDistanceTest <- as.data.frame(do.call("rbind", forDistanceTest), stringsAsFactors = FALSE)
colnames(forDistanceTest) <- c("dist", "group")
forDistanceTest$dist <- as.numeric(forDistanceTest$dist)
forDistanceTest$time <- gsub("_[[:alnum:]]{1,11}$", "", forDistanceTest$group)
forDistanceTest$lab <- gsub("^T[[:digit:]]{1}_", "", forDistanceTest$group)
write.csv(forDistanceTest, paste0(outfilePrefix, ".distances.csv"))
write.csv(distMat, paste0(outfilePrefix, ".pairwiseDistances.csv"))

##############################################################################
### do the test
res <- as.matrix(anova(lm(dist ~ time + lab, data = forDistanceTest)))
write.csv(res, paste0(outfilePrefix, ".anova.csv"))

##############################################################################
### plot the values by group
pValTime <- res["time", "Pr(>F)"]
pValLab <- res["lab", "Pr(>F)"]
pValTime <- ifelse(pValTime < 0.00001, "< 0.00001", paste0("= ", pValTime))
pValLab <- ifelse(pValLab < 0.00001, "< 0.00001", paste0("= ", pValLab))
forMain <- paste0("Ptime ", pValTime, ", Plab ", pValLab) 
forDistanceTest$group <- factor(forDistanceTest$group, levels = paste0(rep(c("T1_", "T2_"), each = 5), rep(c("Hannover1", "Hannover2", "Bern", "Munster", "Zurich"), 2)))
svg(paste0(outfilePrefix, ".boxplot.svg"), height = 5, width = 6)
par(oma = c(3,3,5,1), mar = c(3,3,5,1))
boxplot(forDistanceTest$dist ~ forDistanceTest$group, ylab = "Manhattan distance", xlab = "", las = 1, horizontal = FALSE, main = forMain)
jitterX <- jitter(rep(0, nrow(forDistanceTest)), 10)
points(as.numeric(forDistanceTest$group)+jitterX, forDistanceTest$dist, col = "black", pch = 16, cex = 0.5)
invisible(dev.off())

