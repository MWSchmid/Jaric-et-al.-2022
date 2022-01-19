#!/usr/bin/env Rscript

rm(list=ls())
baseDir <- "/media/mwschmid/myData/MWSchmid/Ivana_ATAC"
scriptDir <- file.path(baseDir, "scripts")
designFile <- file.path(baseDir, "report", "sampleTab.csv")
infileName <- file.path("/home/marc/tempIJA", "union.5.BP.counts.txt")
outfilePrefix <- file.path(baseDir, "GitIgnore_results", "union.5.BP.onlyTimePointOne")
curTime <- "T1"
rm(list=ls())

## arguments from commandline
myarg <- commandArgs()
argPos <- grep("--args", myarg, fixed = TRUE)
curTime <- myarg[argPos+1]
designFile <- myarg[argPos+2]
infileName <- myarg[argPos+3]
outfilePrefix <- myarg[argPos+4]

# rest
suppressPackageStartupMessages({
  options(java.parameters = "-Xmx12g")
  library("XLConnect")
  library("biomaRt")
  library("GenomicRanges")
  library("GenomicFeatures")
  library("limma")
  library("DESeq2")
  library("gtools")
  source("/media/mwschmid/myData/MWSchmid/Development/R/lm_wrapper.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapperAds.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper_overviewBandwidthFix.R")
})

## a function for printing
f.print.message <- function(...) { cat("=== ", format(Sys.time(), "%Y %b %d %X"), paste("===", ...,"\n")) }


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
sampleTab <- subset(sampleTab, time == curTime)
myData <- myData[,rownames(sampleTab)]

#########################################################################################
# remove entries with very low values (<5 in all four samples) - no filter on variance
myData <- f.strip.data(myData, minVal = 5, minTimes = 3, lowerVarQuantileToRemove = 0, colsToStrip = rownames(sampleTab))

#########################################################################################
# normalize and models, use DESeq2 directly to include the batch effect in the group comparisons
formulaString <- "~batchWithinBlock+lab"
dds <- DESeqDataSetFromMatrix(countData = myData, colData = sampleTab, design = formula(formulaString))
dds <- DESeq(dds, minReplicatesForReplace=Inf)

myNormData <- log2(DESeq2::counts(dds, normalized = TRUE)+1)

deResults <- list()

#########################################################################################
# one to one comparisons
allLabs <- unique(sampleTab$lab)
for (i in 1:(length(allLabs)-1)) {
  labA <- allLabs[i]
  for (j in (i+1):length(allLabs)) {
    labB <- allLabs[j]
    curResLabel <- paste0("oto_", curTime, "_", labA, "_vs_", labB)
    mainRes <- results(dds, contrast=c("lab", labA, labB), cooksCutoff=FALSE)
    mainResDf <- data.frame(generic = mainRes$baseMean, logFC = mainRes$log2FoldChange, pVal = mainRes$pvalue, adjP = mainRes$padj, row.names = rownames(mainRes))
    mainResDEGtab <- c.DEGtab(tool = "DESeq", method = "DESeq_default", table = mainResDf, isPairwise = FALSE, pairOrCont = curResLabel)
    deResults[[curResLabel]] <- mainResDEGtab
  }
}

#########################################################################################
# one to many comparisons, requires new fits and variables
for (curLab in allLabs) {
  curResLabel <- paste0("otm_", curTime, "_", curLab, "_vs_others")
  sampleTab[[curResLabel]] <- "other"
  sampleTab[[curResLabel]][sampleTab$lab == curLab] <- curLab
  curFormulaString <- paste0("~0+batchWithinBlock+", curResLabel)
  dds <- DESeqDataSetFromMatrix(countData = myData, colData = sampleTab, design = formula(curFormulaString))
  dds <- DESeq(dds, minReplicatesForReplace=Inf)
  mainRes <- results(dds, contrast=c(curResLabel, curLab, "other"), cooksCutoff=FALSE)
  mainResDf <- data.frame(generic = mainRes$baseMean, logFC = mainRes$log2FoldChange, pVal = mainRes$pvalue, adjP = mainRes$padj, row.names = rownames(mainRes))
  mainResDEGtab <- c.DEGtab(tool = "DESeq", method = "DESeq_default", table = mainResDf, isPairwise = FALSE, pairOrCont = curResLabel)
  deResults[[curResLabel]] <- mainResDEGtab
}

#########################################################################################
# Global F-tests, use batch effect, without outlier refitting
formulaString <- "~ batchWithinBlock + lab"
formulaStringReduced <- "~ batchWithinBlock"
dds <- DESeqDataSetFromMatrix(countData = myData, colData = sampleTab, design = formula(formulaString))
dds <- DESeq(dds, minReplicatesForReplace=Inf, test = "LRT", reduced = as.formula(formulaStringReduced))
mainRes <- results(dds, cooksCutoff=FALSE)
mainResDf <- data.frame(generic = mainRes$baseMean, logFC = mainRes$log2FoldChange, pVal = mainRes$pvalue, adjP = mainRes$padj, row.names = rownames(mainRes))
mainResDEGtab <- c.DEGtab(tool = "DESeq", method = "DESeq_default", table = mainResDf, isPairwise = FALSE, pairOrCont = "LRT: lab")
deResults[["LRT_LAB"]] <- mainResDEGtab

#########################################################################################
# normalize and models, use DESeq2 directly to include the batch effect in the group comparisons
formulaString <- "~lab"
dds <- DESeqDataSetFromMatrix(countData = myData, colData = sampleTab, design = formula(formulaString))
dds <- DESeq(dds, minReplicatesForReplace=Inf)
deResultsNoBatch <- list()

#########################################################################################
# one to one comparisons
allLabs <- unique(sampleTab$lab)
for (i in 1:(length(allLabs)-1)) {
  labA <- allLabs[i]
  for (j in (i+1):length(allLabs)) {
    labB <- allLabs[j]
    curResLabel <- paste0("oto_", curTime, "_", labA, "_vs_", labB)
    mainRes <- results(dds, contrast=c("lab", labA, labB), cooksCutoff=FALSE)
    mainResDf <- data.frame(generic = mainRes$baseMean, logFC = mainRes$log2FoldChange, pVal = mainRes$pvalue, adjP = mainRes$padj, row.names = rownames(mainRes))
    mainResDEGtab <- c.DEGtab(tool = "DESeq", method = "DESeq_default", table = mainResDf, isPairwise = FALSE, pairOrCont = curResLabel)
    deResultsNoBatch[[curResLabel]] <- mainResDEGtab
  }
}

#########################################################################################
# one to many comparisons, requires new fits and variables
for (curLab in allLabs) {
  curResLabel <- paste0("otm_", curTime, "_", curLab, "_vs_others")
  sampleTab[[curResLabel]] <- "other"
  sampleTab[[curResLabel]][sampleTab$lab == curLab] <- curLab
  curFormulaString <- paste0("~", curResLabel)
  dds <- DESeqDataSetFromMatrix(countData = myData, colData = sampleTab, design = formula(curFormulaString))
  dds <- DESeq(dds, minReplicatesForReplace=Inf)
  mainRes <- results(dds, contrast=c(curResLabel, curLab, "other"), cooksCutoff=FALSE)
  mainResDf <- data.frame(generic = mainRes$baseMean, logFC = mainRes$log2FoldChange, pVal = mainRes$pvalue, adjP = mainRes$padj, row.names = rownames(mainRes))
  mainResDEGtab <- c.DEGtab(tool = "DESeq", method = "DESeq_default", table = mainResDf, isPairwise = FALSE, pairOrCont = curResLabel)
  deResultsNoBatch[[curResLabel]] <- mainResDEGtab
}

#########################################################################################
# Global F-tests, use batch effect, without outlier refitting
formulaString <- "~ lab"
formulaStringReduced <- "~ 1"
dds <- DESeqDataSetFromMatrix(countData = myData, colData = sampleTab, design = formula(formulaString))
dds <- DESeq(dds, minReplicatesForReplace=Inf, test = "LRT", reduced = as.formula(formulaStringReduced))
mainRes <- results(dds, cooksCutoff=FALSE)
mainResDf <- data.frame(generic = mainRes$baseMean, logFC = mainRes$log2FoldChange, pVal = mainRes$pvalue, adjP = mainRes$padj, row.names = rownames(mainRes))
mainResDEGtab <- c.DEGtab(tool = "DESeq", method = "DESeq_default", table = mainResDf, isPairwise = FALSE, pairOrCont = "LRT: lab")
deResultsNoBatch[["LRT_LAB"]] <- mainResDEGtab

#########################################################################################
# store results
save(sampleTab, deResults, deResultsNoBatch, myData, myNormData, file = paste0(outfilePrefix, ".Rdata"))
#f.write.DEGtabs.to.workbook(deResults, dirname(outfilePrefix), basename(outfilePrefix), onlySigRawP = TRUE); gc()
#for (curCont in names(deResults)) {
#  sigPeaks <- deResults[[curCont]]$get_significant_entries(0.01, 0.01, 0.5)$any
#  if (length(sigPeaks) < 10) {
#    f.print.message(curCont, "has less than 10 significant peaks, skipping plot.")
#  } else {
#    f.generic.correlation.matrix(myNormData[sigPeaks,], dirname(outfilePrefix), paste0(basename(outfilePrefix), "_sigPeaks_", curCont, "_norm_pearson"), "pearson", FALSE)
#  }
#}
anySigPeak <- unlist(lapply(deResults, function(x) x$get_significant_entries(0.01, 0.01, 0.5)$any))
f.generic.correlation.matrix.alt.col(myNormData[anySigPeak,], dirname(outfilePrefix), paste0(basename(outfilePrefix), "_anySigPeak_norm_pearson"), "pearson", FALSE)

