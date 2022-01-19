import sys
import pysam
from itertools import chain

import pyBigWig
import numpy as np

### usage information
usage = """

python %s flanksize annotationInfile bamInfile/bwInfile

if the infile is a bam, the coverage will be scaled to RPM

annotationInfile:

gene, trans, chrom, strand, start, end, ...

""" % sys.argv[0]

class geneWithCov:
    
    def __init__(self, gene, transcript, chrom, strand, start, end):
        self.gene = gene
        self.transcript = transcript
        self.chrom = chrom
        self.strand = strand
        self.start = int(start)
        self.end = int(end)
        self.cov = [0.0]*(self.end-self.start)
        self.numReadsFound = float(0)
    
    def __str__(self):
        if self.strand == "-":
            self.cov = self.cov[::-1]
        cov = ';'.join([str(x) for x in self.cov])
        out = '\t'.join([self.gene, self.transcript, str(self.numReadsFound), cov])
        return out
    
    def addCoverage(self, bamHandler):
        for al in infile.fetch(self.chrom, self.start, self.end):
            tags = dict(al.tags)
            try:
                hit = float(tags["XW"]) # count weights
            except KeyError:
                hit = 1
            self.numReadsFound += hit
            hasOutOfTranscript = False
            #print >> sys.stderr, '\t'.join([str(al.pos), str(hit)])
            if hit > 0:
                cs = al.pos
                for cigar_entry in al.cigar:
                    cigar_op = cigar_entry[0]
                    cigar_len = cigar_entry[1]
                    if cigar_op == 0:
                        for i in xrange(cs, (cs + cigar_len), 1):
                            try:
                                self.cov[i-self.start] += hit
                            except IndexError:
                                hasOutOfTranscript = True
                        cs += cigar_len
                    else:
                        cs += cigar_len
            #if hasOutOfTranscript:
                #print >> sys.stderr, '\t'.join([self.gene, str(self.start), str(self.end), str(al.pos), str(cs), str(al.pos-self.start), str(self.end-cs)])
    
    def smoothCoveragePostLoading(self, scaleTo, winLen=9, winType="flat"):
        dummy = bigWigConnection("")
        raw = np.array(self.cov)
        raw = raw/float(scaleTo)*1e6
        if winType == "flat":
            out = dummy.smoothCoverageFlat(raw, winLen)
        else:
            out = dummy.smoothCoverageSpec(raw, winLen, winType)
        self.cov = out
        
    def addRawCov(self, bwConnection):
        """Add coverage given a bigWigConnection."""
        self.cov = bwConnection.getRawCoverage(self.chrom, self.start, self.end)
        return None

    def addSmoothCov(self, bwConnection, winLen=147, winType="flat"):
        """Add smoothened coverage given a bigWigConnection."""
        self.cov = bwConnection.getSmoothCoverage(self.chrom, self.start, self.end, winLen, winType)
        return None

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
    flankSize = int(sys.argv[1])
    annotationInfile = sys.argv[2]
    bigWigInfile = sys.argv[3]
except:
    sys.exit(usage)

useRawCoverage = False
if "--useRawCoverage" in sys.argv:
    useRawCoverage = True

# load genes
alreadySeen = set([])
genes = []
with open(annotationInfile, 'rb') as infile:
    for line in infile:
        fields = line.strip().split('\t')
        gene, transcript, chrom, strand, start, end = fields[0:6]
        if gene in alreadySeen:
            continue
        genes.append(geneWithCov(gene, transcript, chrom, strand, int(start)-flankSize, int(end)+flankSize))
        alreadySeen.add(gene)

if bigWigInfile.endswith(".bw"):
    BWconn = bigWigConnection(bigWigInfile)
    for gene in genes:
        if useRawCoverage:
            gene.addRawCov(BWconn)
        else:
            gene.addSmoothCov(BWconn, winLen=147, winType="flat")
        print gene
elif bigWigInfile.endswith(".bam"):
    # get total number of hits first
    totalHits = 0
    infile = pysam.Samfile(bigWigInfile, "rb")
    for al in infile.fetch():
        tags = dict(al.tags)
        try:
            totalHits += float(tags["XW"]) # count weights
        except KeyError:
            totalHits += 1
    infile.close()
    print >> sys.stderr, "scaling to %10.2f hits" % (totalHits)
    # now the covarage
    infile = pysam.Samfile(bigWigInfile, "rb")
    for gene in genes:
        gene.addCoverage(infile)
        if not useRawCoverage:
            gene.smoothCoveragePostLoading(scaleTo=totalHits, winLen=300, winType="flat")
        print gene
    infile.close()
else:
    sys.exit(usage)

