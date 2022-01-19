#!/usr/bin/env Rscript

rm(list=ls())
fileType <- "chipseeker"
maxP <- 0.05
enrichFile <- "/media/mwschmid/myData/MWSchmid/Ivana_ATAC/GitIgnore_results/enrichment/union.3.BP.bothTimePoints.oto_T1.sigPeaks_LFC_1.0.chipSeeker_KEGG.csv"
outfileName <- "/media/mwschmid/myData/MWSchmid/Ivana_ATAC/GitIgnore_results/enrichment/union.3.BP.bothTimePoints.oto_T1.sigPeaks_LFC_1.0.chipSeeker_KEGG.png"
rm(list=ls())

## arguments from commandline
myarg <- commandArgs()
argPos <- grep("--args", myarg, fixed = TRUE)

fileType <- myarg[argPos+1] # interminer, chipseeker, topgo
maxP <- as.numeric(myarg[argPos+2])
enrichFile <- myarg[argPos+3]
outfileName <- myarg[argPos+4]

fileType <- tolower(fileType)

suppressPackageStartupMessages({
  library("gplots")
  #library("RNAseqWrapper")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapperAds.R")
})
f.print.message <- function(...) { cat("=== ", format(Sys.time(), "%Y %b %d %X"), paste("===", ...,"\n")) }

f.print.message(enrichFile)

if (fileType == "interminer") {
  pValEntry <- "pValue"
  symbolEntry <- "symbol"
  termEntry <- "description"
  temp <- read.csv(enrichFile, header = TRUE, stringsAsFactors = FALSE, row.names = 1)
} else if (fileType == "chipseeker") {
  pValEntry <- "pvalue"
  symbolEntry <- "symbol"
  termEntry <- "Description"
  temp <- read.csv(enrichFile, header = TRUE, stringsAsFactors = FALSE, row.names = 1)
} else if (fileType == "topgo") {
  pValEntry <- "weight"
  symbolEntry <- "genesInTermSymbol"
  termEntry <- "Term"
  temp <- read.table(enrichFile, header = TRUE, stringsAsFactors = FALSE, sep = '\t', row.names = 1)
  # make sure the p-values are correct
  temp[[pValEntry]][temp[[pValEntry]] == "< 1e-30"] <- 1e-30
  if (!is.numeric(temp[[pValEntry]])) {
    before <- as.numeric(temp[[pValEntry]][,1])
    temp[[pValEntry]] <- matrix(as.numeric(temp[[pValEntry]]), ncol = ncol(temp$Term), dimnames = dimnames(temp$Term))
    after <- as.numeric(temp[[pValEntry]][,1])
    if (sum(before == after) != length(before)) {
      f.print.message("this should not happen")
    }
  }
} else {
  f.print.message("!!!!! file type not supported, valid are interminer, chipseeker, or topgo!!!!!\n")
  quit("no", 0)
}

# keep only significant entries
temp <- temp[temp[[pValEntry]] < maxP,c(pValEntry, symbolEntry, termEntry)]; rownames(temp) <- NULL
colnames(temp) <- c("pVal", "genes", "term")

# get a matrix with gene vs entry, filled with P-value
allGenes <- unique(unlist(sapply(temp$genes, function(x) unlist(strsplit(x, ";")))))
oneGenePerRow <- lapply(allGenes, function(x) cbind(temp[grep(x, temp$genes),], x))
oneGenePerRow <- do.call("rbind", oneGenePerRow)
colnames(oneGenePerRow)[colnames(oneGenePerRow) == "x"] <- "gene"
pValMat <- f.long.to.wide(oneGenePerRow, "gene", "term", "pVal")
pValMat[is.na(pValMat)] <- 1
pValMat <- (-log10(pValMat))

# sort the columns by sum
sigCounts <- colSums(pValMat>0)
pValMat <- pValMat[,names(sort(sigCounts))]
dim(pValMat)

# do the plot
if ((nrow(pValMat) > 0)& (ncol(pValMat)>0)) {
  f.open.figure(dirname(outfileName), basename(outfileName), TRUE, height = 10+0.05*nrow(pValMat), width = 10+0.05*ncol(pValMat))
  par(mar = c(0,0,0,0), oma = c(10,0,0,0))
  heatmap.2(as.matrix(pValMat), col = colorRampPalette(c("white", "darkorchid4"))(11), trace="none", scale = "none", margins = c(15,15), dendrogram = "row", Colv = FALSE)
  #heatmap.2(as.matrix(pValMat), col = f.blackredyellow(11), trace="none", scale = "none", margins = c(15,15), dendrogram = "row", Colv = FALSE)
  f.close.figure()
} else {
  f.print.message("skipping plots, no entries left.")
}

# vars for testing
#enrichFile <- "/media/mwschmid/myData/MWSchmid/MarkusStoffelDiffSplice/report/enrichment/PATHWAY_manuallyCurated.csv"
#enrichFile <- "/media/mwschmid/myData/MWSchmid/MarkusStoffelDiffSplice/report/enrichment/KEGG_manuallyCurated.csv"
#enrichFile <- "/media/mwschmid/myData/MWSchmid/MarkusStoffelDiffSplice/report/enrichment/REACTOME_manuallyCurated.csv"
#enrichFile <- "/media/mwschmid/myData/MWSchmid/MarkusStoffelDiffSplice/report/enrichment/manuallyCurated_BP_minNodeSize_1.txt"
#rDir <- "/media/mwschmid/myData/MWSchmid/MarkusStoffelDiffSplice/report/figures"
#outfileName <- "tescht.svg"

