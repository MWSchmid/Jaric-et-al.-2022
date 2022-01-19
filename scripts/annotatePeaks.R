#!/usr/bin/env Rscript

rm(list=ls())
baseDir <- "/media/mwschmid/myData/MWSchmid/Ivana_ATAC"
scriptDir <- file.path(baseDir, "scripts")
infileName <- file.path(baseDir, "GitIgnore_peaks_merged", "union.3.BP.bed")
outfileName <- file.path(baseDir, "GitIgnore_peaks_merged", "union.3.BP.csAnnot")
rm(list=ls())

## arguments from commandline
myarg <- commandArgs()
argPos <- grep("--args", myarg, fixed = TRUE)
infileName <- myarg[argPos+1]
outfileName <- myarg[argPos+2]

skipExisting <- FALSE
if ("--skipExisting" %in% myarg) {
  skipExisting <- TRUE
}

# rest
suppressPackageStartupMessages({
  options(java.parameters = "-Xmx12g")
  library("XLConnect")
  library("biomaRt")
  library("ChIPseeker")
  library("GenomicRanges")
  library("GenomicFeatures")
  library("limma")
  library("DESeq2")
  library("gtools")
  library("org.Mm.eg.db")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapperAds.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper_overviewBandwidthFix.R")
})

## a function for printing
f.print.message <- function(...) { cat("=== ", format(Sys.time(), "%Y %b %d %X"), paste("===", ...,"\n")) }

#########################################################################################
# connect to mart first
ensembl <- f.connect.to.ensembl(martDataSet = "mmusculus_gene_ensembl", martName = "ENSEMBL_MART_ENSEMBL")
#txdb <- f.get.txdb(martDataSet = "mmusculus_gene_ensembl", martName = "ENSEMBL_MART_ENSEMBL") # doesn't work as function
if (skipExisting & file.exists(outfileName)) {
  out <- read.csv(outfileName, header = TRUE, stringsAsFactors = FALSE, row.names = 1)
} else {
  foundError <- TRUE
  numberOfTries <- 0
  while (foundError) {
    txdb <- tryCatch(makeTxDbFromBiomart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "www.ensembl.org"),
                     error = function(e) {NA}, finally = cat("### txdb: trying again\n"))
    if (!is.na(txdb)) {
      foundError <- FALSE
    } else {
      numberOfTries <- numberOfTries + 1
    }
    if (numberOfTries > 500) {
      f.print.message("ENSEMBL: tried more than 500 times to reach it, something is wrong.")
      break
    }
  }
  
  #########################################################################################
  # some annotation simplification and ordering
  annotationTranslator <- c("3' UTR" = "3pUTR", "5' UTR" = "5pUTR", "Distal Intergenic" = "distalIntergenic", 
                            "Downstream (<1kb)" = "DS1", "Downstream (1-2kb)" = "DS2", "Downstream (2-3kb)" = "distalIntergenic",
                            "Promoter (<=1kb)" = "prom_1", "Promoter (1-2kb)" = "prom_2", "Exon" = "exon", "Intron" = "intron",
                            "FirstIntron" = "firstIntron", "FirstExon" = "firstExon")
  annotationOrder <- c("distalIntergenic", "USTSS2", "USTSS1", "DSTSS1", "DSTSS2", "5pUTR", "firstExon", "firstIntron", "exon", "intron", "3pUTR", "DS1", "DS2")
  
  #########################################################################################
  # load data and check for multiple entrie
  cseData <- read.table(infileName, sep = '\t', header = FALSE, col.names = c("chrom", "start", "end", "name", "score", "strand"))
  cseData$chrom <- gsub("^chr", "", cseData$chrom)
  if (sum(table(paste0(cseData$chrom, '_', cseData$start, '_', cseData$end))>1) > 0) {
    f.print.message("ERROR: multiple peaks at the same location are not allowed!")
    quit("no", 0)
  }
  
  #########################################################################################
  # annotate peaks (simple table)
  f.print.message("identified peaks, doing annotation.")
  cseDataGranges <- makeGRangesFromDataFrame(cseData[,c("chrom", "start", "end")], seqnames.field="chrom", start.field="start", end.field=c("end"), ignore.strand=TRUE, seqinfo=NULL, starts.in.df.are.0based=FALSE) # MACS output is not 0-based
  cseDataAnnotated <- as.data.frame(annotatePeak(cseDataGranges, tssRegion = c(-5e2, 2e3), TxDb = txdb, level = "gene"), stringsAsFactors = FALSE)
  rownames(cseDataAnnotated) <- with(cseDataAnnotated, paste(seqnames, start, end, sep = '_'))
  out <- cseData; rownames(out) <- with(cseData, paste(chrom, start, end, sep = '_'))
  out$geneID <- cseDataAnnotated[rownames(out), "geneId"]
  out$anno <- cseDataAnnotated[rownames(out), "annotation"]
  temp <- out$anno
  temp[grep("intron 1 of", temp)] <- "FirstIntron"
  temp[grep("exon 1 of", temp)] <- "FirstExon"
  temp[grep("^Intron", temp)] <- "Intron"
  temp[grep("^Exon", temp)] <- "Exon"
  out$simpleAnno <- annotationTranslator[temp]
  out$distToTSS <- cseDataAnnotated[rownames(out), "distanceToTSS"]
  out$simpleAnno[(out$distToTSS > -2000) & out$distToTSS < -1000] <- "USTSS2"
  out$simpleAnno[(out$distToTSS > -1000) & out$distToTSS <= 0] <- "USTSS1"
  out$simpleAnno[(out$distToTSS < 2000) & out$distToTSS > 1000] <- "DSTSS2"
  out$simpleAnno[(out$distToTSS < 1000) & out$distToTSS > 0] <- "DSTSS1"
  write.csv(out, outfileName)
}

#########################################################################################
# add the gene description
if (!(skipExisting & file.exists(paste0(outfileName, ".withDescription")))) {
  uniqueGenes <- unique(out$geneID)
  geneDescr <- unique(getBM(attributes=c("ensembl_gene_id", "external_gene_name", "strand", "description"), filters="ensembl_gene_id", values=uniqueGenes, mart=ensembl))
  rownames(geneDescr) <- geneDescr$ensembl_gene_id
  strandTranslator <- c("1" = "+", "-1" = "-")
  out$geneName <- geneDescr[out$geneID, "external_gene_name"]
  out$strand <- strandTranslator[as.character(geneDescr[out$geneID, "strand"])]
  out$desc <- geneDescr[out$geneID, "description"]
  write.csv(out, paste0(outfileName, ".withDescription"))
}





