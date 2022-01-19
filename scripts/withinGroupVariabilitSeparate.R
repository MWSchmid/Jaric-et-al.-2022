#!/usr/bin/env Rscript

rm(list=ls())
baseDir <- "/media/mwschmid/myData/MWSchmid/Ivana_ATAC"
infileName <- file.path(baseDir, "GitIgnore_results", "union.5.BP.bothTimePoints.Rdata")
outfilePrefix <- file.path(baseDir, "GitIgnore_results", "union.5.BP.withinGroupVariabilitySeparate")
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
write.csv(forDistanceTest, paste0(outfilePrefix, ".distances.csv"), row.names = FALSE)
distMatT1 <- distMat[grep("^S", rownames(distMat), value = TRUE), grep("^S", colnames(distMat), value = TRUE)]
distMatT2 <- distMat[grep("^T", rownames(distMat), value = TRUE), grep("^T", colnames(distMat), value = TRUE)]
write.csv(distMatT1, paste0(outfilePrefix, ".pairwiseDistances.T1.csv"))
write.csv(distMatT2, paste0(outfilePrefix, ".pairwiseDistances.T2.csv"))

##############################################################################
### do the test
for (curTimePoint in unique(forDistanceTest$time)) {
  subDat <- subset(forDistanceTest, time == curTimePoint)
  res <- as.matrix(anova(lm(dist ~ lab, data = subDat)))
  write.csv(res, paste0(outfilePrefix, ".anova.", curTimePoint, ".csv"))
  pValLab <- res["lab", "Pr(>F)"]
  pValLab <- ifelse(pValLab < 0.00001, "< 0.00001", paste0("= ", pValLab))
  forMain <- paste0("P ", pValLab) 
  subDat$group <- factor(subDat$group, levels = paste0(curTimePoint, "_", c("Hannover1", "Hannover2", "Bern", "Munster", "Zurich")))
  svg(paste0(outfilePrefix, ".boxplot.", curTimePoint, ".svg"), height = 5, width = 3)
  par(oma = c(3,3,5,1), mar = c(3,3,5,1))
  boxplot(subDat$dist ~ subDat$group, ylab = "Manhattan distance", xlab = "", las = 1, horizontal = FALSE, main = forMain)
  jitterX <- jitter(rep(0, nrow(subDat)), 10)
  points(as.numeric(subDat$group)+jitterX, subDat$dist, col = "black", pch = 16, cex = 0.5)
  invisible(dev.off())
}


