#!/usr/bin/env Rscript

rm(list=ls())
enrichPath <- "/home/marc/tempIJA"
rm(list=ls())

## arguments from commandline
myarg <- commandArgs()
argPos <- grep("--args", myarg, fixed = TRUE)
enrichPath <- myarg[argPos + 1]
rDir <- enrichPath

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
  library("gplots")
  library("ggplot2")
  #library("igraph")
  #library("poweRlaw")
  #library("org.Mm.eg.db")
  #source("/media/mwschmid/myData/MWSchmid/Development/R/OTU_functions.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapperAds.R")
})

## a function for printing
f.print.message <- function(...) { cat("=== ", format(Sys.time(), "%Y %b %d %X"), paste("===", ...,"\n")) }

## some vars
filePattern <- list()
gsubPattern <- list()
fileHeader <- list()
setOrder <- list()
dataDirs <- list()
maxPlist <- list()
minSizeList <- list()

thresholdsToCheck <- c("0.0", "0.5", "1.0", "1.5", "2.0")
minNodeSizeToCheck <- c(1, 5)
for (curCategory in c("BP", "MF", "CC")) {
  for (curThresh in thresholdsToCheck) {
    for (minNodeSize in minNodeSizeToCheck) {
      curKey <- paste0(curCategory, "_minLFC_", curThresh, "_minNodeSize_", minNodeSize)
      #"ALL_union.2.BP.bothTimePoints.time.sigPeaks_LFC_1.0.homer_BP_minNodeSize_1.txt.gz"
      filePattern[[curKey]] <- paste0("ALL_union\\.3\\.BP\\.bothTimePoints\\.ot.*\\.sigPeaks_LFC_", curThresh, "\\.homer_", curCategory, "_minNodeSize_", minNodeSize)
      gsubPattern[[curKey]] <- paste0("ALL_union\\.3\\.BP\\.bothTimePoints\\.|\\.sigPeaks_LFC_", curThresh, "\\.homer_", curCategory, "_minNodeSize_", minNodeSize)
      # this down here is always the same
      fileHeader[[curKey]] <- c("Term", "Annotated", "Significant", "Expected", "weight", "genesInTermSymbol", "genesInTermEnsembl")
      setOrder[[curKey]] <- c("oto_T1", "oto_T2", "otm_T1", "otm_T1_Hannover1", "otm_T1_Hannover2", "otm_T1_Munster", "otm_T1_Zurich", "otm_T1_Bern", "otm_T2", "otm_T2_Hannover1", "otm_T2_Hannover2", "otm_T2_Munster", "otm_T2_Zurich", "otm_T2_Bern")
      dataDirs[[curKey]] <- enrichPath
      maxPlist[[curKey]] <- 0.05
      minSizeList[[curKey]] <- 1
    }
  }
}

## read everything
pValEntry <- "weight"
goResults <- list()
for (curSet in names(dataDirs)) {
  temp <- f.read.files.from.directory.noStringsAsFactors(dataDirs[[curSet]], filePattern[[curSet]], ".txt", "", fileHeader[[curSet]], "\t", useHeader = TRUE)
  for (curVar in names(temp)) {
    colnames(temp[[curVar]]) <- gsub(gsubPattern[[curSet]], "", colnames(temp[[curVar]]))
    names(colnames(temp[[curVar]])) <- NULL
  }
  temp[[pValEntry]][temp[[pValEntry]] == "< 1e-30"] <- 1e-30
  if (!is.numeric(temp[[pValEntry]])) {
    before <- as.numeric(temp[[pValEntry]][,1])
    temp[[pValEntry]] <- matrix(as.numeric(temp[[pValEntry]]), ncol = ncol(temp$Term), dimnames = dimnames(temp$Term))
    after <- as.numeric(temp[[pValEntry]][,1])
    if (sum(before == after) != length(before)) {
      f.print.message("this should not happen")
    }
  }
  #temp[[pValEntry]] <- as.numeric(temp[[pValEntry]])
  #temp[[pValEntry]][is.na(temp[[pValEntry]])] <- 1e-30
  if (nrow(temp$Term) > 0) {
    goResults[[curSet]] <- temp
  }
}

selTerms <- list()
for (curSet in names(goResults)) {
  selTerms[[curSet]] <- f.plot.topGO.enrichment(goResults[[curSet]], setOrder[[curSet]], rDir, paste0("topGO_enrichment_", curSet, ".svg"), minSize = minSizeList[[curSet]], maxP = maxPlist[[curSet]], maxTerms = 100, pValEntry = pValEntry)
}

for (curSet in names(selTerms)) {
  commonCols <- intersect(setOrder[[curSet]], colnames(goResults[[curSet]]$Significant))
  if (length(commonCols) < length(setOrder[[curSet]])) {
    f.print.message("less common columns than expected:")
    #cat(paste(setdiff(setOrder[[curSet]], colnames(goResults[[curSet]]$Significant)), collapse = '\n'), '\n')
  }
  temp <- goResults[[curSet]]$Significant[selTerms[[curSet]],commonCols]
  tempSig <- goResults[[curSet]][[pValEntry]][selTerms[[curSet]],commonCols]
  curRange <- range(temp[tempSig < maxPlist[[curSet]]], na.rm = TRUE)
  f.print.message(curSet, curRange[1], curRange[2])
}


intervals <- c(-30, -20, seq(-10, 0, by = 2))
colorGrad <- f.yellowredblack(length(intervals))
f.open.figure(rDir, "topGO_pValueLegend.svg", TRUE, height = 3, width = 4)
plot(1:length(colorGrad), rep(1, length(colorGrad)), col = colorGrad, type = "h")
f.close.figure()


