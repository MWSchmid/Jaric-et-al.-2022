#!/usr/bin/env Rscript

rm(list=ls())
baseDir <- "/media/mwschmid/myData/MWSchmid/Ivana_ATAC"
infileName <- file.path(baseDir, "GitIgnore_results", "union.5.BP.bothTimePoints.Rdata")
chipAnnoInfile <- "/home/marc/tempIJA/union.5.BP.csAnnot.csv"
homeAnnoInfile <- "/home/marc/tempIJA/union.5.BP.homer.csv"
outfilePrefix <- file.path(baseDir, "GitIgnore_results", "geneSets", "union.5.BP.bothTimePoints")
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
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapperAds.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper_overviewBandwidthFix.R")
})

ensembl <- f.connect.to.ensembl(martDataSet = "mmusculus_gene_ensembl", martName = "ENSEMBL_MART_ENSEMBL")

##############################################################################
### load file
load(infileName) # sampleTab, deResults, (deResultsNoBatch), myData, myNormData, (design, myCont)
csAnno <- read.csv(chipAnnoInfile, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
hoAnno <- read.csv(homeAnnoInfile, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
hoAnno$Chr <- gsub("^chr", "", hoAnno$Chr)
rownames(hoAnno) <- paste0(hoAnno$Chr, '_', hoAnno$Start-1, '_', hoAnno$End)

##############################################################################
### collect significant peaks
allContrasts <- c(grep("^oto", names(deResults), value = TRUE),
                  grep("^otm", names(deResults), value = TRUE),
                  grep("^time", names(deResults), value = TRUE),
                  grep("^LRT", names(deResults), value = TRUE))
numberSigPeaks <- matrix(0, nrow = length(allContrasts), ncol = 3, dimnames = list(allContrasts, c("any", "up", "down")))
out <- list()
for (curCont in allContrasts) {
  #sigPeakList <- deResults[[curCont]]$get_significant_entries(0.01, 0.01, 1)
  sigPeakList <- deResults[[curCont]]$get_significant_entries(0.01, 0.01, 0.5)
  #sigPeakList <- deResults[[curCont]]$get_significant_entries(0.01, 0.01, 0.0)
  for (curSet in names(sigPeakList)) {
    sigPeaks <- sigPeakList[[curSet]]
    numberSigPeaks[curCont, curSet] <- length(sigPeaks)
    if (length(sigPeaks) == 0) {
      f.print.message("no significant peaks in", curCont, curSet)
      next
    }
    out[[paste0(curCont, '_', curSet)]] <- cbind(rep(curCont, length(sigPeaks)), rep(curSet, length(sigPeaks)), gsub("^chr", "", sigPeaks), round(deResults[[curCont]]$table[sigPeaks, "logFC"],3))
  }
}
out <- as.data.frame(do.call("rbind", out), stringsAsFactors = FALSE)
colnames(out) <- c("contrast", "direction", "peak", "logFC")
out$contrast <- as.character(out$contrast)
out$direction <- as.character(out$direction)
out$peak <- as.character(out$peak)
out$logFC <- as.numeric(out$logFC)

##############################################################################
### add annotation info
out$chipSeekerAnno <- csAnno[out$peak, "simpleAnno"]
out$chipSeekerDistToTSS <- csAnno[out$peak, "distToTSS"]
out$chipSeekerGene <- csAnno[out$peak, "geneID"]
out$homerGene <- hoAnno[out$peak, "GeneID"]
out$homerDistToTSS <- hoAnno[out$peak, "DistanceToTSS"]

##############################################################################
### add gene description
uniqueGenes <- unique(c(out$homerGene, out$chipSeekerGene))
geneDescr <- unique(getBM(attributes=c("ensembl_gene_id", "external_gene_name", "strand", "description"), filters="ensembl_gene_id", values=uniqueGenes, mart=ensembl))
rownames(geneDescr) <- geneDescr$ensembl_gene_id
out$chipSeekerName <- geneDescr[out$chipSeekerGene, "external_gene_name"]
out$chipSeekerDescription <- geneDescr[out$chipSeekerGene, "description"]
out$homerName <- geneDescr[out$homerGene, "external_gene_name"]
out$homerDescription <- geneDescr[out$homerGene, "description"]

##############################################################################
### write output and subsets
write.csv(numberSigPeaks, paste0(outfilePrefix, ".numberSigPeaks.csv"))
outOrder <- c("peak", "contrast", "direction", "logFC", "homerDistToTSS", "homerGene", "homerName", "homerDescription", "chipSeekerAnno", "chipSeekerDistToTSS", "chipSeekerGene", "chipSeekerName", "chipSeekerDescription")
write.csv(out[,outOrder], paste0(outfilePrefix, ".allGenes.csv"), row.names = FALSE)
outLess <- subset(out, (direction == "any") & ((abs(chipSeekerDistToTSS) < 2000) | (abs(homerDistToTSS) < 2000)))
write.csv(outLess[,outOrder], paste0(outfilePrefix, ".onlyCloseToTSS.csv"), row.names = FALSE)



