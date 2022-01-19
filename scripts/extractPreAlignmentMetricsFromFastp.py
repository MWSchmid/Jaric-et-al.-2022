#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
usage:
python extractPreAlignmentMetricsFromFastp.py fastqResults.json

numbers are written to stdout
"""

import json
import sys

# check input
try:
    fastpFileName = sys.argv[1]
except:
    print __doc__
    sys.exit(1)

with open(fastpFileName, "rb") as infile:
    fastqResults = json.load(infile)
    # before filtering
    numSeqs_BF = int(fastqResults["summary"]["before_filtering"]["total_reads"])/2 # this is otherwise both ends summed up, not done like this in the test samples
    seqLenForward_BF = float(fastqResults["summary"]["before_filtering"]["read1_mean_length"])
    seqLenReverse_BF = float(fastqResults["summary"]["before_filtering"]["read2_mean_length"])
    percGC_BF = float(fastqResults["summary"]["before_filtering"]["gc_content"])*100
    percQ20_BF = float(fastqResults["summary"]["before_filtering"]["q20_rate"])*100
    percQ30_BF = float(fastqResults["summary"]["before_filtering"]["q30_rate"])*100
    # after filtering
    numSeqs_AF = int(fastqResults["summary"]["after_filtering"]["total_reads"])/2 # this is otherwise both ends summed up, not done like this in the test samples
    seqLenForward_AF = float(fastqResults["summary"]["after_filtering"]["read1_mean_length"])
    seqLenReverse_AF = float(fastqResults["summary"]["after_filtering"]["read2_mean_length"])
    percGC_AF = float(fastqResults["summary"]["after_filtering"]["gc_content"])*100
    percQ20_AF = float(fastqResults["summary"]["after_filtering"]["q20_rate"])*100
    percQ30_AF = float(fastqResults["summary"]["after_filtering"]["q30_rate"])*100
    # common
    insertSize = int(fastqResults["insert_size"]["peak"])
    percDup = float(fastqResults["duplication"]["rate"])*100
    percAdapterReads = int(fastqResults["adapter_cutting"]["adapter_trimmed_reads"])/float(numSeqs_BF*2)*100 # here it's again both ends summed up
    numPassedFilter = int(fastqResults["filtering_result"]["passed_filter_reads"])
    percLowQualityReads = int(fastqResults["filtering_result"]["low_quality_reads"])/float(numSeqs_BF*2)*100 # here it's again both ends summed up
    percTooManyNReads = int(fastqResults["filtering_result"]["too_many_N_reads"])/float(numSeqs_BF*2)*100 # here it's again both ends summed up

headerList = [
    "numberOfReadsBeforeFilter",
    "averageLengthForwardBeforeFilter",
    "averageLengthReverseBeforeFilter",
    "percentGCBeforeFilter",
    "percentQ20BasesBeforeFilter",
    "percentQ30BasesBeforeFilter",
    "numberOfReadsAfterFilter",
    "averageLengthForwardAfterFilter",
    "averageLengthReverseAfterFilter",
    "percentGCAfterFilter",
    "percentQ20BasesAfterFilter",
    "percentQ30BasesAfterFilter",
    "insertSize",
    "percentDuplication",
    "percentReadsWithAdapters",
    "numberOfReadsPassingFilter",
    "percentLowQualityReads",
    "percentHighNBaseContent"]
    
valueList = [
    numSeqs_BF,
    seqLenForward_BF,
    seqLenReverse_BF,
    percGC_BF,
    percQ20_BF,
    percQ30_BF,
    numSeqs_AF,
    seqLenForward_AF,
    seqLenReverse_AF,
    percGC_AF,
    percQ20_AF,
    percQ30_AF,
    insertSize,
    percDup,
    percAdapterReads,
    numPassedFilter,
    percLowQualityReads,
    percTooManyNReads]

outHeader = ','.join([str(x) for x in headerList])
sys.stdout.write(outHeader + ';')
outString = ','.join([str(x) for x in valueList])
sys.stdout.write(outString)
