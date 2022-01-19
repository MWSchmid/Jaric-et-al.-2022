#!/usr/bin/env Rscript

rm(list=ls())
baseDir <- "/media/mwschmid/myData/MWSchmid/Ivana_ATAC"
scriptDir <- file.path(baseDir, "scripts")
designFile <- file.path(baseDir, "report", "sampleTab.csv")
infileName <- file.path("/home/marc/tempIJA", "union.5.BP.counts.txt")
outfilePrefix <- file.path(baseDir, "GitIgnore_results", "union.5.BP.bothTimePoints")
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
  library("XLConnect")
  library("biomaRt")
  library("GenomicRanges")
  library("GenomicFeatures")
  library("limma")
  library("DESeq2")
  library("gtools")
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
# remove entries with very low values (<5 in all four samples) - no filter on variance
myData <- f.strip.data(myData, minVal = 5, minTimes = 3, lowerVarQuantileToRemove = 0, colsToStrip = rownames(sampleTab))

#########################################################################################
# a grouping variable
sampleTab$G <- with(sampleTab, paste(time, lab, sep = '_'))
sampleTab$G__ <- factor(sampleTab$G) # don't specify levels, we anyway do the contrasts by hand
formulaString <- "~ batch + G__"

#########################################################################################
# normalize the data with DESeq2
myNormData <- f.normalize.counts.DESeq(myData, sampleTab, formulaString)

#########################################################################################
# models, use DESeq2 directly to include the batch effect
deResults <- list()

# labs, global F-test, can't use time because not full rank otherwise
formulaStringReduced <- "~ batch"
dds <- DESeqDataSetFromMatrix(countData = myData, colData = sampleTab, design = formula(formulaString))
ddsLRT <- DESeq(dds, minReplicatesForReplace=Inf, test = "LRT", reduced = as.formula(formulaStringReduced))
mainRes <- results(ddsLRT, cooksCutoff=FALSE)
mainResDf <- data.frame(generic = mainRes$baseMean, logFC = mainRes$log2FoldChange, pVal = mainRes$pvalue, adjP = mainRes$padj, row.names = rownames(mainRes))
mainResDEGtab <- c.DEGtab(tool = "DESeq", method = "DESeq_default", table = mainResDf, isPairwise = FALSE, pairOrCont = "LRT: group")
deResults[["LRT_LAB"]] <- mainResDEGtab

dds <- DESeq(dds, minReplicatesForReplace=Inf)

# one to one comparisons
allLabs <- unique(sampleTab$lab)
for (curTime in unique(sampleTab$time)) {
  for (i in 1:(length(allLabs)-1)) {
    labA <- allLabs[i]
    labAwithTime <- paste0(curTime, '_', labA)
    for (j in (i+1):length(allLabs)) {
      labB <- allLabs[j]
      labBwithTime <- paste0(curTime, '_', labB)
      mainRes <- results(dds, contrast=c("G__", labAwithTime, labBwithTime), cooksCutoff=FALSE)
      mainResDf <- data.frame(generic = mainRes$baseMean, logFC = mainRes$log2FoldChange, pVal = mainRes$pvalue, adjP = mainRes$padj, row.names = rownames(mainRes))
      mainResDEGtab <- c.DEGtab(tool = "DESeq", method = "DESeq_default", table = mainResDf, isPairwise = FALSE, pairOrCont = paste0("oto_", curTime, "_", labA, "_vs_", labB))
    }
  }
}

# one to many comparisons
  allLabs <- unique(sampleTab$lab)
  analysisSets <- list()
  analysisSetsWithBatch <- list()
  for (curLab in allLabs) {
    curCont <- paste0("otm_T1_", curLab, "_vs_others")
    analysisSets[[curCont]] <- list(formulaString = paste0("~ 0 + ", curCont))
    analysisSetsWithBatch[[curCont]] <- list(formulaString = paste0("~ 0 + batchWithinBlock + ", curCont))
    sampleTab[[curCont]] <- "low"
    sampleTab[[curCont]][sampleTab$lab == curLab] <- "high"
  }
  
  otmRes <- f.run.DESeq.models.test.only.last.term.in.parallel.highest.vs.lowest(analysisSets, myData, sampleTab, 1)
  otmResWithBatch <- f.run.DESeq.models.test.only.last.term.in.parallel.highest.vs.lowest(analysisSets, myData, sampleTab, 1)








allTimePoints <- unique(sampleTab$time)
allLabs <- unique(sampleTab$lab)
design <- model.matrix(formula(formulaString), data = sampleTab, contrasts.arg = NULL)

# Define contrasts
# one to one within time point
# one to many within time point
# T2 vs T1 given lab
# T2 vs T1 over all lavs
forContrastMaking <- list()
for (curTime in allTimePoints) {
  for (curLab in allLabs) {
    otherLabs <- setdiff(allLabs, curLab)
    curTimeAndLab <- paste0(curTime, '_', curLab)
    curTimeAndOtherLabs <- paste0(paste0(curTime, '_', otherLabs), collapse = '|')
    forContrastMaking[[paste0("otm_", curTime, '_', curLab, "_vs_others")]] <- f.formulate.simple.contrast(curTimeAndLab, curTimeAndOtherLabs, design)
    for (otherLab in otherLabs) {
      otherTimeAndLab <- paste0(curTime, '_', otherLab)
      forContrastMaking[[paste0("oto_", curTimeAndLab, "_vs_", otherLab)]] <- f.formulate.simple.contrast(curTimeAndLab, otherTimeAndLab, design, fixed = TRUE)
    }
  }
}
for (curLab in allLabs) {
  forContrastMaking[[paste0("time_", curLab, "_T2_vs_T1")]] <- f.formulate.simple.contrast(paste0("T2_", curLab), paste0("T1_", curLab), design, fixed = TRUE)
}
forContrastMaking[[paste0("time_all_T2_vs_T1")]] <- f.formulate.simple.contrast("T2_", "T1_", design, fixed = TRUE)

# check for empty fields/cases
forContrastMakingCleaned <- list()
for (curCont in names(forContrastMaking)) {
  searchForVoid <- grep("()", forContrastMaking[[curCont]], fixed = TRUE, value = TRUE)
  if (length(searchForVoid)>0) {
    f.print.message(paste0("Missing fields for contrast: ", curCont)) # none by now
    #f.print.message(paste0("Original formula: ", forContrastMaking[[curCont]]))
  } else {
    forContrastMakingCleaned[[curCont]] <- forContrastMaking[[curCont]]
  }
}
forContrastMakingCleaned <- unlist(forContrastMakingCleaned)
myCont <- makeContrasts(contrasts=forContrastMakingCleaned, levels = design) # the contrast= is important
colnames(myCont) <- names(forContrastMakingCleaned)
write.csv(myCont, paste0(outfilePrefix, "_contrastsForManualCheck.csv"))



#########################################################################################
# store results
save(sampleTab, deResults, myData, myNormData, file = paste0(outfilePrefix, ".Rdata"))
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
f.generic.correlation.matrix(myNormData[anySigPeak,], dirname(outfilePrefix), paste0(basename(outfilePrefix), "_anySigPeak_norm_pearson"), "pearson", FALSE)








