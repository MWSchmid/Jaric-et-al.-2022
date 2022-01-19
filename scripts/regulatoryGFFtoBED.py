# -*- coding: utf-8 -*-
"""
usage:
python regulatoryGFFtoBED.py infileName outDir
"""

import sys
import gzip

try:
    infileName = sys.argv[1]
    outDir = sys.argv[2]
except:
    print >> sys.stderr, __doc__
    sys.exit(1)

def myopen(fileName, mode="r"):
    """open either a regular or a compressed file"""
    if fileName.endswith(".gz"):
        return gzip.open(fileName, mode=mode)
    else:
        return open(fileName, mode=mode)

outFileHandler = {}
featureTypeCounter = {}
featureTypeSizes = {}
with myopen(infileName, "r") as infile:
    for line in infile:
        (chrom, source, feature, start, end, score, strand, phase, rest) = line.strip().split('\t')
        restDict = dict([x.split('=') for x in rest.split(';')])
        curID = "NA"
        curType = "NA"
        for key, value in restDict.items():
            if key == "ID":
                curID = value
            #if key == "feature_type" or key == "stable_id":
            if key == "feature_type" or key == "motif_feature_type":
                curType = value.translate(None, ' ')
            if key == "binding_matrix":
                curType = value.split('_')[0]
        if curID == "NA" and curType == "NA":
            print >> sys.stderr, line
            sys.exit(1)
        try:
            featureTypeCounter[curType] += 1
        except KeyError:
            featureTypeCounter[curType] = 1
        if curID == "NA":
            curID = '_'.join([curType, str(featureTypeCounter[curType])])
        outString = '\t'.join([chrom, start, end, curID, score, strand])
        featureSize = abs(int(end)-int(start))
        try:
            print >> outFileHandler[curType], outString
            featureTypeSizes[curType].append(featureSize)
        except KeyError:
            outFileHandler[curType] = open(outDir+'/'+''.join([curType, ".bed"]), "wb")
            print >> outFileHandler[curType], outString
            featureTypeSizes[curType] = [featureSize]

# close outfiles and close them as well
for curType, outfile in outFileHandler.items():
    outfile.close()
    maxSize = max(featureTypeSizes[curType])
    minSize = min(featureTypeSizes[curType])
    meanSize = sum(featureTypeSizes[curType])/float(len(featureTypeSizes[curType]))
    below100 = 0
    for elem in featureTypeSizes[curType]:
        if elem < 100:
            below100 += 1
    print >> sys.stderr, "%s: %d\t%d\t%4.1f\t%d\t%3.1f" % (curType, featureTypeCounter[curType], minSize, meanSize, maxSize, (below100/float(featureTypeCounter[curType]))*100)

