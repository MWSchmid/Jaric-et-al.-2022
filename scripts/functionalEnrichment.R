#!/usr/bin/env Rscript

rm(list=ls())
baseDir <- "/media/mwschmid/myData/MWSchmid/Ivana_ATAC"
infileName <- file.path(baseDir, "GitIgnore_results", "union.3.BP.bothTimePoints.Rdata")
chipAnnoInfile <- "/home/marc/tempIJA/union.3.BP.csAnnot.csv"
homeAnnoInfile <- "/home/marc/tempIJA/union.3.BP.homer.csv"
outfilePrefix <- file.path(baseDir, "GitIgnore_results", "enrichment", "union.3.BP.bothTimePoints")
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
  library("org.Mm.eg.db")
  library("topGO")
  library("clusterProfiler")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapperAds.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper_overviewBandwidthFix.R")
})

ensembl <- f.connect.to.ensembl(martDataSet = "mmusculus_gene_ensembl", martName = "ENSEMBL_MART_ENSEMBL")

rDir <- dirname(outfilePrefix)
outfilePrefixNoPath <- basename(outfilePrefix)

##############################################################################
### load file
load(infileName) # sampleTab, deResults, (deResultsNoBatch), myData, myNormData, (design, myCont)
csAnno <- read.csv(chipAnnoInfile, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
hoAnno <- read.csv(homeAnnoInfile, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
hoAnno$Chr <- gsub("^chr", "", hoAnno$Chr)
rownames(hoAnno) <- paste0(hoAnno$Chr, '_', hoAnno$Start-1, '_', hoAnno$End)

##############################################################################
### collect significant peaks
contrastList <- list(
  oto_T1 = grep("^oto_T1", names(deResults), value = TRUE),
  oto_T2 = grep("^oto_T2", names(deResults), value = TRUE),
  otm_T1 = grep("^otm_T1", names(deResults), value = TRUE),
  otm_T2 = grep("^otm_T2", names(deResults), value = TRUE),
  otm_T1_Hannover1 = grep("^otm_T1_Hannover1", names(deResults), value = TRUE),
  otm_T2_Hannover1 = grep("^otm_T2_Hannover1", names(deResults), value = TRUE),
  otm_T1_Hannover2 = grep("^otm_T1_Hannover2", names(deResults), value = TRUE),
  otm_T2_Hannover2 = grep("^otm_T2_Hannover2", names(deResults), value = TRUE),
  otm_T1_Munster = grep("^otm_T1_Munster", names(deResults), value = TRUE),
  otm_T2_Munster = grep("^otm_T2_Munster", names(deResults), value = TRUE),
  otm_T1_Zurich = grep("^otm_T1_Zurich", names(deResults), value = TRUE),
  otm_T2_Zurich = grep("^otm_T2_Zurich", names(deResults), value = TRUE),
  otm_T1_Bern = grep("^otm_T1_Bern", names(deResults), value = TRUE),
  otm_T2_Bern = grep("^otm_T2_Bern", names(deResults), value = TRUE),
  time = grep("^time", names(deResults), value = TRUE)#,
  #lrt = grep("^LRT", names(deResults), value = TRUE)
)

allPeaks <- list()
for (curSet in names(contrastList)) {
  if (length(contrastList[[curSet]]) > 0) {
    tempPeakList <- list(
      sigPeaks_LFC_0.0 = Reduce(union, lapply(deResults[contrastList[[curSet]]], function(x) x$get_significant_entries(0.01, 0.01, 0.0)$any)),
      sigPeaks_LFC_0.5 = Reduce(union, lapply(deResults[contrastList[[curSet]]], function(x) x$get_significant_entries(0.01, 0.01, 0.5)$any)),
      sigPeaks_LFC_1.0 = Reduce(union, lapply(deResults[contrastList[[curSet]]], function(x) x$get_significant_entries(0.01, 0.01, 1.0)$any)),
      sigPeaks_LFC_1.5 = Reduce(union, lapply(deResults[contrastList[[curSet]]], function(x) x$get_significant_entries(0.01, 0.01, 1.5)$any)),
      sigPeaks_LFC_2.0 = Reduce(union, lapply(deResults[contrastList[[curSet]]], function(x) x$get_significant_entries(0.01, 0.01, 2.0)$any))
    )
    tempBedList <- lapply(tempPeakList, function(x) as.data.frame(do.call("rbind", strsplit(x, '_')), stringsAsFactors = FALSE))
    for (curPeaks in names(tempBedList)) {
      if (nrow(tempBedList[[curPeaks]]) == 0) { next }
      colnames(tempBedList[[curPeaks]]) <- c("chrom", "start", "end")
      tempBedList[[curPeaks]]$name <- paste0(curSet, '_', curPeaks, '_', 1:nrow(tempBedList[[curPeaks]]))
      forRowNames <- with(tempBedList[[curPeaks]], paste0(gsub("chr", "", chrom), '_', start, '_', end))
      tempBedList[[curPeaks]]$chipSeekerAnno <- csAnno[forRowNames, "simpleAnno"]
      tempBedList[[curPeaks]]$chipSeekerDistToTSS <- csAnno[forRowNames, "distToTSS"]
      tempBedList[[curPeaks]]$chipSeekerGene <- csAnno[forRowNames, "geneID"]
      tempBedList[[curPeaks]]$homerGene <- hoAnno[forRowNames, "GeneID"]
      tempBedList[[curPeaks]]$homerDistToTSS <- hoAnno[forRowNames, "DistanceToTSS"]
    }
    allPeaks[[curSet]] <- tempBedList
  }
}

##############################################################################
### extract gene lists
mappedGenes <- mappedkeys(org.Mm.egPATH2EG)
keggToEnsembl <- as.list(org.Mm.egPATH2EG[mappedGenes])
geneUniverses <- list(chipSeeker = unique(csAnno$geneID[nchar(csAnno$geneID)>0]), homer = unique(hoAnno$GeneID[nchar(hoAnno$GeneID)>0]))
geneUniversesKEGG <- list()
entrezToEnsemblByAnno <- list()
entrezToSymbolByAnno <- list()
ensemblToEntrezByAnno <- list()
for (annotationType in names(geneUniverses)) {
  temp <- geneUniverses[[annotationType]]
  geneDescription <- unique(getBM(attributes=c("ensembl_gene_id", "mgi_symbol", "entrezgene_id", "description"), filters='ensembl_gene_id', values=temp, mart=ensembl))
  geneDescription <- subset(geneDescription, !is.na(entrezgene_id))
  forMappingsEntrez <- unique(geneDescription$entrezgene_id)
  entrezToEnsembl <- sapply(forMappingsEntrez, function(x) paste0(unique(geneDescription$ensembl_gene_id[geneDescription$entrezgene_id == x]), collapse = '|')); names(entrezToEnsembl) <- as.character(forMappingsEntrez)
  entrezToSymbol <- sapply(forMappingsEntrez, function(x) paste0(unique(geneDescription$mgi_symbol[geneDescription$entrezgene_id == x]), collapse = '|')); names(entrezToSymbol) <- as.character(forMappingsEntrez)
  forMappingsEnsembl <- unique(geneDescription$ensembl_gene_id)
  ensemblToEntrez <- sapply(forMappingsEnsembl, function(x) as.character(unique(geneDescription$entrezgene_id[geneDescription$ensembl_gene_id == x]))); names(ensemblToEntrez) <- as.character(forMappingsEnsembl)
  geneUniversesKEGG[[annotationType]] <- as.character(forMappingsEntrez)
  entrezToEnsemblByAnno[[annotationType]] <- entrezToEnsembl
  entrezToSymbolByAnno[[annotationType]] <- entrezToSymbol
  ensemblToEntrezByAnno[[annotationType]] <- ensemblToEntrez
}

#curSet <- names(allPeaks)[1]
for (curSet in names(allPeaks)) {
  tempBedList <- allPeaks[[curSet]]
  #curPeaks <- names(tempBedList)[1]
  for (curPeaks in names(tempBedList)) {
    subData <- tempBedList[[curPeaks]]
    if (nrow(subData) == 0) { next }
    #annotationType <- "homer"
    for (annotationType in c("homer", "chipSeeker")) {
      geneTestSet <- subData[[paste0(annotationType, "Gene")]][abs(subData[[paste0(annotationType, "DistToTSS")]]) <= 2000]
      geneTestSet <- unique(geneTestSet[nchar(geneTestSet) > 0])
      geneTestSet <- geneTestSet[!is.na(geneTestSet)]
      if (length(geneTestSet) < 20) { next }
      enrichPrefix <- paste0(outfilePrefixNoPath, ".", curSet, ".", curPeaks, ".", annotationType)
      # GO term enrichment
      if (file.exists(file.path(rDir, paste0(enrichPrefix, "_BP_minNodeSize_5.txt")))|file.exists(file.path(rDir, paste0(enrichPrefix, "_BP_minNodeSize_5.txt.gz")))) {
        f.print.message("skipping topGO run (BP)")
      } else {
        topGOresults <- f.topGO.ensembl(geneUniverses[[annotationType]], geneTestSet, rDir, enrichPrefix, "BP", "org.Mm.eg.db", mart = ensembl, addGenes = TRUE, minNodeSize = 1, symbolColumn = "mgi_symbol", useLegacyEntrezgene = FALSE)
        topGOresults <- f.topGO.ensembl(geneUniverses[[annotationType]], geneTestSet, rDir, enrichPrefix, "BP", "org.Mm.eg.db", mart = ensembl, addGenes = TRUE, minNodeSize = 5, symbolColumn = "mgi_symbol", useLegacyEntrezgene = FALSE)
      }
      if (file.exists(file.path(rDir, paste0(enrichPrefix, "_MF_minNodeSize_5.txt")))|file.exists(file.path(rDir, paste0(enrichPrefix, "_MF_minNodeSize_5.txt.gz")))) {
        f.print.message("skipping topGO run (MF)")
      } else {
        topGOresults <- f.topGO.ensembl(geneUniverses[[annotationType]], geneTestSet, rDir, enrichPrefix, "MF", "org.Mm.eg.db", mart = ensembl, addGenes = TRUE, minNodeSize = 1, symbolColumn = "mgi_symbol", useLegacyEntrezgene = FALSE)
        topGOresults <- f.topGO.ensembl(geneUniverses[[annotationType]], geneTestSet, rDir, enrichPrefix, "MF", "org.Mm.eg.db", mart = ensembl, addGenes = TRUE, minNodeSize = 5, symbolColumn = "mgi_symbol", useLegacyEntrezgene = FALSE)
      }
      if (file.exists(file.path(rDir, paste0(enrichPrefix, "_CC_minNodeSize_5.txt")))|file.exists(file.path(rDir, paste0(enrichPrefix, "_CC_minNodeSize_5.txt.gz")))) {
        f.print.message("skipping topGO run (CC)")
      } else {
        topGOresults <- f.topGO.ensembl(geneUniverses[[annotationType]], geneTestSet, rDir, enrichPrefix, "CC", "org.Mm.eg.db", mart = ensembl, addGenes = TRUE, minNodeSize = 1, symbolColumn = "mgi_symbol", useLegacyEntrezgene = FALSE)
        topGOresults <- f.topGO.ensembl(geneUniverses[[annotationType]], geneTestSet, rDir, enrichPrefix, "CC", "org.Mm.eg.db", mart = ensembl, addGenes = TRUE, minNodeSize = 5, symbolColumn = "mgi_symbol", useLegacyEntrezgene = FALSE)
      }
      # KEGG
      if (file.exists(file.path(rDir, paste0(enrichPrefix, "_KEGG.csv")))) {
        f.print.message("skipping KEGG run")
        next
      }
      toTestKEGG <- unlist(ensemblToEntrez[geneTestSet])
      toTestKEGG <- as.character(toTestKEGG[!is.na(toTestKEGG)])
      allPossibleGenesKEGG <- geneUniversesKEGG[[annotationType]]
      entrezToEnsembl <- entrezToEnsemblByAnno[[annotationType]]
      entrezToSymbol <- entrezToSymbolByAnno[[annotationType]]
      ensemblToEntrez <- ensemblToEntrezByAnno[[annotationType]]
      out <- enrichKEGG(toTestKEGG, organism = "mmu", universe = allPossibleGenesKEGG, pAdjustMethod = "none", pvalueCutoff = 0.05, qvalueCutoff = 1) # , pvalueCutoff = 0.1, qvalueCutoff = 0.2
      out <- out[,c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID")]
      out$geneID <- gsub("/", ";", out$geneID, fixed = TRUE) # no clue why, but Rpl10 goes missing at the enrichKEGG step
      if (is.null(out)) { next }
      if (length(out) == 0) { next }
      if (nrow(out) == 0) { next }
      # add gene IDs
      out$symbol <- NA
      out$ensembl <- NA
      for (i in 1:nrow(out)) {
        temp <- unlist(strsplit(out$geneID[i], ";", fixed = TRUE))
        out$symbol[i] <- paste0(entrezToSymbol[temp], collapse = ';')
        out$ensembl[i] <- paste0(entrezToEnsembl[temp], collapse = ';')
      }
      write.csv(out, file.path(rDir, paste0(enrichPrefix, "_KEGG.csv")))
      out <- enrichKEGG(toTestKEGG, organism = "mmu", universe = allPossibleGenesKEGG, pAdjustMethod = "none", pvalueCutoff = 1, qvalueCutoff = 1) # , pvalueCutoff = 0.1, qvalueCutoff = 0.2
      out <- out[,c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID")]
      out$geneID <- gsub("/", ";", out$geneID, fixed = TRUE) # no clue why, but Rpl10 goes missing at the enrichKEGG step
      if (is.null(out)) { next }
      if (length(out) == 0) { next }
      if (nrow(out) == 0) { next }
      # add gene IDs
      out$symbol <- NA
      out$ensembl <- NA
      for (i in 1:nrow(out)) {
        temp <- unlist(strsplit(out$geneID[i], ";", fixed = TRUE))
        out$symbol[i] <- paste0(entrezToSymbol[temp], collapse = ';')
        out$ensembl[i] <- paste0(entrezToEnsembl[temp], collapse = ';')
      }
      write.csv(out, file.path(rDir, paste0("ALL_", enrichPrefix, "_KEGG.csv")))
    }
  }
}

