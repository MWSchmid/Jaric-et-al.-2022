import sys
import pysam
from itertools import chain

import pyBigWig
import numpy as np

### usage information
usage = """

python %s folderWithBigWigs txtFile

txtFile:

* tab separated, two columns:
* sample ID (must match the bigwig file name up to the extension)
* group

""" % sys.argv[0]

class bigWigConnection(object):
    """A bigWig handler"""
    def __init__(self, path):
        self.path = path
    
    def getChromSizes(self):
        infile = pyBigWig.open(self.path)
        temp = infile.chroms()
        out = dict([[key, int(value)] for key, value in temp.items()])
        infile.close()
        return out
    
    def smoothCoverageSave(self, covArray, winLen=9, winType="flat"):
        """Calculate smoothened coverage, see http://www.scipy.org/Cookbook/SignalSmooth.
        Note that winLen should be an odd integer. winType == flat meand moving average."""
        x = covArray
        if x.size < winLen:
            raise ValueError("Input vector needs to be bigger than window size.")
        if x.ndim != 1:
            raise ValueError("Smooth only accepts 1 dimension arrays.")
        if window not in ["flat", "hanning", "hamming", "bartlett", "blackman"]:
            raise ValueError("Valid winTypes are flat, hanning, hamming, bartlett, and blackman.")
        s = np.r_[2 * x[0] - x[winLen:1:-1], x, 2 * x[-1] - x[-1:-winLen:-1]]
        if window == "flat":
            w = np.ones(winLen, 'd')
        else:
            w = eval("np." + window + "(winLen)")
        y = np.convolve(w / w.sum(), s, mode="same")
        return y[winLen-1:-winLen+1]
    
    def smoothCoverageFlat(self, covArray, winLen=9):
        """Calculate smoothened coverage using moving average.
        Note that winLen should be an odd integer."""
        if (winLen % 2) == 0:
            winLen += 1
        x = covArray
        s = np.r_[2 * x[0] - x[winLen:1:-1], x, 2 * x[-1] - x[-1:-winLen:-1]]
        w = np.ones(winLen, 'd')
        y = np.convolve(w / w.sum(), s, mode="same")
        return y[winLen-1:-winLen+1]
    
    def smoothCoverageSpec(self, covArray, winLen=9, winType="blackman"):
        """Calculate smoothened coverage with a specified window type.
        Note that winLen should be an odd integer."""
        if (winLen % 2) == 0:
            winLen += 1
        x = covArray
        s = np.r_[2 * x[0] - x[winLen:1:-1], x, 2 * x[-1] - x[-1:-winLen:-1]]
        w = eval("np." + winType + "(winLen)")
        y = np.convolve(w / w.sum(), s, mode="same")
        return y[winLen-1:-winLen+1]
    
    def getRawCoverage(self, chrom, start, end):
        """Retrieve an array with the genome coverage."""
        infile = pyBigWig.open(self.path)
        temp = infile.values(chrom, start, end)
        out = np.array(temp)
        infile.close()
        return out
    
    def getSmoothCoverage(self, chrom, start, end, winLen=9, winType="flat"):
        """Retrieve an array with the smoothened genome coverage.
        CHECK/TODO: one might have to extend the coverage array first"""
        raw = self.getRawCoverage(chrom, start, end)
        if winType == "flat":
            out = self.smoothCoverageFlat(raw, winLen)
        else:
            out = self.smoothCoverageSpec(raw, winLen, winType)
        return out

try:
    bigWigFolder = sys.argv[1]
    sampleTabFile = sys.argv[2]
except:
    sys.exit(usage)

# load groupings
sampleGroups = {}
with open(sampleTabFile, 'r') as infile:
    for line in infile:
        sampleID, group = line.strip().split('\t')
        filename = bigWigFolder + '/' + sampleID + ".bw"
        try:
            sampleGroups[group].append(filename )
        except KeyError:
            sampleGroups[group] = [filename]

# average them
for group, fileList in sampleGroups.items():
    outfileName = bigWigFolder + '/' + group + '.bw'
    outfile = pyBigWig.open(outfileName, "w")
    infile = bigWigConnection(fileList[0])
    chromSizes = infile.getChromSizes()
    headerList = [(chrom, size) for chrom, size in chromSizes.items()]
    outfile.addHeader(headerList)
    numFiles = len(fileList)
    for chrom, size in headerList:
        size = int(size)
        chromCov = np.zeros(size)
        for curFile in fileList:
            infile = bigWigConnection(curFile)
            tempCov = infile.getRawCoverage(chrom, 0, size)
            chromCov += tempCov
        aveCov = chromCov/float(numFiles)
        #ignore = len(aveCov)%10
        #binnedCov = np.mean(aveCov[:-ignore].reshape(-1, 10), axis=1)
        #if ignore:
        #    toAdd = np.mean(aveCov[-ignore:])
        #    binnedCov = np.append(binnedCov, toAdd)
        #outfile.addEntries(chrom, 0, values=binnedCov, span=10, step=10)
        outfile.addEntries(chrom, 0, values=aveCov, span=1, step=1)
    outfile.close()


#bigWigFolder="/home/marc/tempJV"
#group = "custom"
#fileList = ["/home/marc/tempJV/A13.bw", "/home/marc/tempJV/A14.bw", "/home/marc/tempJV/A15.bw", "/home/marc/tempJV/A17.bw", "/home/marc/tempJV/A18.bw"]

#outfileName = bigWigFolder + '/' + group + '.bw'
#outfile = pyBigWig.open(outfileName, "w")
#infile = bigWigConnection(fileList[0])
#chromSizes = infile.getChromSizes()
#headerList = [(chrom, size) for chrom, size in chromSizes.items()]
#outfile.addHeader(headerList)
#numFiles = len(fileList)
#for chrom, size in headerList:
#    size = int(size)
#    chromCov = np.zeros(size)
#    for curFile in fileList:
#        infile = bigWigConnection(curFile)
#        tempCov = infile.getRawCoverage(chrom, 0, size)
#        chromCov += tempCov
#    aveCov = chromCov/float(numFiles)
#    outfile.addEntries(chrom, 0, values=aveCov, span=1, step=1)
#outfile.close()







