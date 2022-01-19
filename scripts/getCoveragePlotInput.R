#!/usr/bin/env Rscript

## arguments from commandline
myarg <- commandArgs()
argPos <- grep("--args", myarg, fixed = TRUE)

# packages
suppressPackageStartupMessages({
  library("biomaRt")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapperAds.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper_overviewBandwidthFix.R")
})

# some directories
GOIfile <- myarg[argPos+1]
outfileName <- myarg[argPos+2]
onlyLongest <- TRUE
if ("--allTranscripts" %in% myarg){
  onlyLongest <- FALSE
}

#GOIfile <- "/media/mwschmid/myData/MWSchmid/MarkusStoffel_ATAC/report/felipesGOIs.csv"
#outfileName <- "/media/mwschmid/myData/MWSchmid/MarkusStoffel_ATAC/GitIgnore_covPlots/forCoveragePlot.txt"

# connect to ensembl
ensembl <- f.connect.to.ensembl(martDataSet = "mmusculus_gene_ensembl", martName = "ENSEMBL_MART_ENSEMBL")

########################################################################
# load data
goiTab <- read.csv(GOIfile, stringsAsFactors = FALSE, header = TRUE)
rownames(goiTab) <- goiTab$ensemblID

#hurz <- listAttributes(ensembl)
#hurz[grep("chrom.*start", hurz$name),]

forBedUns <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "chromosome_name", "strand", "start_position", "end_position",
                                  "exon_chrom_start", "exon_chrom_end", "transcript_start", "transcript_end",
                                  "gene_biotype"),
                   filters = c("ensembl_gene_id", "biotype"), values = list(rownames(goiTab), "protein_coding"), mart = ensembl)
forBedUns <- subset(forBedUns, (!is.na(exon_chrom_start)) & (!is.na(exon_chrom_end)))
forBedUns$size <- forBedUns$exon_chrom_end - forBedUns$exon_chrom_start
forBed <- forBedUns[order(forBedUns$ensembl_gene_id, forBedUns$ensembl_transcript_id, forBedUns$exon_chrom_start),]

########################################################################
# get the representative/longest transcript
if (onlyLongest) {
  aggrList <- list(chrom = forBed$chromosome_name, strand = forBed$strand, gene = forBed$ensembl_gene_id, trans = forBed$ensembl_transcript_id)
  minSite <- aggregate(forBed$exon_chrom_start, by = aggrList, function(x) min(x, na.rm = TRUE))
  maxSite <- aggregate(forBed$exon_chrom_end, by = aggrList, function(x) max(x, na.rm = TRUE))
  colnames(minSite)[5] <- "start"
  minSite$end <- maxSite$x
  
  f.get.longest <- function(x) {
    if (is.vector(x) | nrow(x) == 1) {
      return(x)
    }
    return(x[which.max(x$size),])
  }
  minSite$size <- minSite$end-minSite$start
  minSiteLongestOnly <- lapply(split(minSite, minSite$gene), f.get.longest)
  minSiteLongestOnlyOut <- do.call("rbind", minSiteLongestOnly)
  selectedTranscripts <- minSiteLongestOnlyOut$trans
  forBed <- subset(forBed, ensembl_transcript_id %in% selectedTranscripts)
  f.get.out.format <- function(x) {
    out <- c(x[1,c(1:4,9,10)], # gene, transcript, chromosome, strand, transcript start, transcript end
             paste0(sort(x[,7], decreasing = FALSE), collapse = ','), # genomicStarts
             paste0(sort(x[,8], decreasing = FALSE), collapse = ',') # genomicEnds
    ) 
    return(as.character(out))
  }
} else {
  f.get.out.format <- function(x) {
    out <- c(x[1,c(1:4,5,6)], # gene, transcript, chromosome, strand, gene start, gene end
             paste0(sort(x[,7], decreasing = FALSE), collapse = ','), # genomicStarts
             paste0(sort(x[,8], decreasing = FALSE), collapse = ',') # genomicEnds
    ) 
    return(as.character(out))
  }
}

########################################################################
# filter and format output

out <- lapply(split(forBed, forBed$ensembl_transcript_id), f.get.out.format)
out <- do.call("rbind", out)
colnames(out) <- c("gene", "trans", "chrom", "strand", "start", "end", "genomicStart", "genomicEnd")
out[,"strand"][out[,"strand"] == "1"] <- "+"
out[,"strand"][out[,"strand"] == "-1"] <- "-"
out <- as.data.frame(out, stringsAsFactors = FALSE)
out$start <- as.numeric(out$start)
out <- out[with(out, order(chrom, start)),]
out$gene <- goiTab[out$gene, "symbol"]
write.table(out, outfileName, sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)





