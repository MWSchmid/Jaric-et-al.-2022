#!/usr/bin/env Rscript

rm(list=ls())
fileType <- "chipseeker"
outfilePrefix <- "/home/marc/tempIJA/homer_KEGG"
allFiles <- file.path("/home/marc/tempIJA", c("T1.txtOrCsv", "T2.txtOrCsv"))

fileType <- "topgo"
outfilePrefix <- "/home/marc/tempIJA/homer_BP_minNodeSize_5"
allFiles <- file.path("/home/marc/tempIJA", c("T1.txtOrCsv", "T2.txtOrCsv"))

rm(list=ls())

## arguments from commandline
myarg <- commandArgs()
argPos <- grep("--args", myarg, fixed = TRUE)
fileType <- myarg[argPos+1] # chipseeker, topgo
outfilePrefix <- myarg[argPos+2]
allFiles <- myarg[(argPos+3):length(myarg)]
pValThresh <- 0.05

keepOrder <- FALSE
if ("--keepOrder" %in% myarg) {
  keepOrder <- TRUE
  allFiles <- allFiles[1:(which(allFiles == "--keepOrder")-1)]
  orderFile <- myarg[which(myarg == "--keepOrder")+1]
  orderToUse <- scan(orderFile, what = "character")
}

# rest
suppressPackageStartupMessages({
  options(java.parameters = "-Xmx12g")
  library("gtools")
  library("gplots")
  library("ggplot2")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapperAds.R")
})

## a function for printing
f.print.message <- function(...) { cat("=== ", format(Sys.time(), "%Y %b %d %X"), paste("===", ...,"\n")) }

## get set names
setNames <- gsub("\\.txt|\\.csv|\\.txtOrCsv", "", basename(allFiles))
names(allFiles) <- setNames

## check type of data
allResults <- list()
if (fileType == "chipseeker") {
  pValEntry <- "pvalue"
  symbolEntry <- "symbol"
  termEntry <- "Description"
  for (curSet in names(allFiles)) {
    temp <- read.csv(allFiles[curSet], header = TRUE, stringsAsFactors = FALSE, row.names = 1)
    temp$ratioTestSet <- unlist(lapply(strsplit(temp$GeneRatio, '/', fixed = TRUE), function(x) as.numeric(x[1])/as.numeric(x[2])))
    temp$ratioBackground <- unlist(lapply(strsplit(temp$BgRatio, '/', fixed = TRUE), function(x) as.numeric(x[1])/as.numeric(x[2])))
    temp$enrichment <- temp$ratioTestSet/temp$ratioBackground
    temp$log2enrichment <- log2(temp$enrichment)
    temp$isSignificant <- temp[[pValEntry]] < pValThresh
    allResults[[curSet]] <- temp
  }
} else if (fileType == "topgo") {
  pValEntry <- "weight"
  symbolEntry <- "genesInTermSymbol"
  termEntry <- "Term"
  for (curSet in names(allFiles)) {
    temp <- read.table(allFiles[curSet], header = TRUE, stringsAsFactors = FALSE, sep = '\t', row.names = 1)
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
    temp$enrichment <- temp$Significant/temp$Expected
    temp$log2enrichment <- log2(temp$enrichment)
    temp$isSignificant <- temp[[pValEntry]] < pValThresh
    allResults[[curSet]] <- temp
  }
} else {
  f.print.message("!!!!! file type not supported, valid are chipseeker or topgo!!!!!\n")
  #quit("no", 0)
}

## merge different samples into one matrix
allTerms <- Reduce(union, lapply(allResults, rownames))
termToDescription <- rep("unknown", length(allTerms)); names(termToDescription) <- allTerms
enrichmentMatrix <- matrix(NA, ncol = length(setNames), nrow = length(allTerms), dimnames = list(allTerms, setNames))
enrichmentMatrixOnlySig <- matrix(NA, ncol = length(setNames), nrow = length(allTerms), dimnames = list(allTerms, setNames))
for (curSet in names(allResults)) {
  curTerms <- rownames(allResults[[curSet]])
  enrichmentMatrix[curTerms, curSet] <- allResults[[curSet]]$enrichment
  enrichmentMatrixOnlySig[curTerms, curSet] <- allResults[[curSet]]$enrichment
  enrichmentMatrixOnlySig[curTerms[!(allResults[[curSet]]$isSignificant)], curSet] <- NA
  termToDescription[curTerms] <- allResults[[curSet]][curTerms, termEntry]
}
if (keepOrder) {
  presentInOrder <- sum(orderToUse %in% rownames(enrichmentMatrix))
  if (presentInOrder < (length(orderToUse)-1)) {
    f.print.message("missing terms:")
    orderToUse[!(orderToUse %in% rownames(enrichmentMatrix))]
  }
  orderToUse <- orderToUse[orderToUse %in% rownames(enrichmentMatrix)]
  enrichmentMatrix <- enrichmentMatrix[orderToUse,]
  enrichmentMatrixOnlySig <- enrichmentMatrixOnlySig[orderToUse,]
} else {
  noNAcount <- sort(rowSums(!is.na(enrichmentMatrixOnlySig)), decreasing = TRUE)
  onTop <- names(noNAcount)[noNAcount == ncol(enrichmentMatrix)]
  absSumTop <- sort(rowSums(abs(enrichmentMatrixOnlySig[onTop,]), na.rm = TRUE), decreasing = TRUE)
  above <- names(absSumTop)
  others <- setdiff(rownames(enrichmentMatrix), onTop)
  absSum <- sort(rowSums(abs(enrichmentMatrixOnlySig[others,]), na.rm = TRUE), decreasing = TRUE)
  below <- names(absSum)
  if (ncol(enrichmentMatrix) == 2) {
    belowFirst <- sort(enrichmentMatrixOnlySig[below,1][!is.na(enrichmentMatrixOnlySig[below,1])], decreasing = TRUE)
    belowSecond <- sort(enrichmentMatrixOnlySig[below,2][!is.na(enrichmentMatrixOnlySig[below,2])], decreasing = TRUE)
    below <- c(names(belowFirst), names(belowSecond))
  }
  enrichmentMatrix <- enrichmentMatrix[c(above, below),]
  enrichmentMatrixOnlySig <- enrichmentMatrixOnlySig[c(above, below),]
}
out <- data.frame(enrichmentMatrix, description = termToDescription[rownames(enrichmentMatrix)])
write.csv(out, paste0(outfilePrefix, ".enrichMatrix.csv"))
out <- data.frame(enrichmentMatrixOnlySig, description = termToDescription[rownames(enrichmentMatrixOnlySig)])
write.csv(out, paste0(outfilePrefix, ".enrichMatrixOnlySignificant.csv"))

# we neet to revert the order for the plot
enrichmentMatrix <- enrichmentMatrix[rev(rownames(enrichmentMatrix)),]
enrichmentMatrixOnlySig <- enrichmentMatrixOnlySig[rev(rownames(enrichmentMatrixOnlySig)),]

## scatter plot with two samples or enrichment matrix with two or more samples
if ((ncol(enrichmentMatrix) == 2) & (sum(is.na(enrichmentMatrix)) == 0)) {
  equalLimits <- range(pretty(range(enrichmentMatrix)))
  forColor <- rep("black", nrow(enrichmentMatrix)); names(forColor) <- rownames(enrichmentMatrix)
  sigInX <- rownames(enrichmentMatrixOnlySig)[!is.na(enrichmentMatrixOnlySig[,1])]
  sigInY <- rownames(enrichmentMatrixOnlySig)[!is.na(enrichmentMatrixOnlySig[,2])]
  forColor[sigInX] <- "darkorchid3"
  forColor[sigInY] <- "darkgoldenrod3"
  forColor[intersect(sigInX, sigInY)] <- "black"
  pdf(paste0(outfilePrefix, ".scatterplot.pdf"), height = 5, width = 5)
  plot(enrichmentMatrix[,1], enrichmentMatrix[,2], xlim = equalLimits, ylim = equalLimits, pch = 16, las = 1, xlab = paste0("Enrichment in ", colnames(enrichmentMatrix)[1]), ylab = paste0("Enrichment in ", colnames(enrichmentMatrix)[2]), col = forColor)
  abline(a=0, b=1, lty = 3)
  invisible(dev.off())
}

## heatmap with stars if significant
toPlot <- log2(enrichmentMatrix)
if (sum(!is.na(toPlot)) == 0) {
  maxVal <- 0.1
} else {
  maxVal <- max(abs(toPlot[!is.na(toPlot)]))+0.2*max(abs(toPlot[!is.na(toPlot)]))
}
colorGradient <- colorRampPalette(c("orangered3", "white", "slateblue3"))(31)
xChars <- colnames(toPlot)
yChars <- termToDescription[rownames(toPlot)]
xPos <- 1:length(xChars)
yPos <- 1:length(yChars)
pdf(paste0(outfilePrefix, ".heatmap.pdf"), height = 4+nrow(enrichmentMatrix)/6, width = 8)
par(oma = c(1,26,1,1))
image(xPos, yPos, t(toPlot), xlab = "", ylab = "", main = "", yaxt = "n", xaxt = "n", zlim = c(-maxVal, maxVal), col = colorGradient)
axis(1, at = xPos, labels = xChars, outer = FALSE, las = 1)
axis(2, at = yPos, labels = yChars, outer = FALSE, las = 1)
forText <- round(enrichmentMatrix, 2)
for (curCol in colnames(forText)) {
  forText[,curCol][!is.na(enrichmentMatrixOnlySig[,curCol])] <- paste0(forText[,curCol][!is.na(enrichmentMatrixOnlySig[,curCol])], " *")
}
text(rep(xPos, each = length(yPos)), yPos, forText, cex = 1)
invisible(dev.off())

## heatmap with significant only
toPlot <- log2(enrichmentMatrixOnlySig)
pdf(paste0(outfilePrefix, ".heatmapOnlySignificant.pdf"), height = 4+nrow(enrichmentMatrix)/6, width = 8)
par(oma = c(1,26,1,1))
image(xPos, yPos, t(toPlot), xlab = "", ylab = "", main = "", yaxt = "n", xaxt = "n", zlim = c(-maxVal, maxVal), col = colorGradient)
axis(1, at = xPos, labels = xChars, outer = FALSE, las = 1)
axis(2, at = yPos, labels = yChars, outer = FALSE, las = 1)
forText <- round(enrichmentMatrixOnlySig, 2)
text(rep(xPos, each = length(yPos)), yPos, forText, cex = 1)
invisible(dev.off())










quit("no", 0)





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


