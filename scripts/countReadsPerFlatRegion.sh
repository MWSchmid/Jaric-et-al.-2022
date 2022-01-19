#!/bin/bash
#
# countReadsPerFlatRegion.sh -- count the number of reads within simple genomic regions. 
# This is not meant for RNA-Seq read counting, as there is no exact matching and hierarchy
# taken care of. It's useful for ChIP-Seq data (count reads per peak/region).
#
# Authors: Marc W. Schmid <contact@mwschmid.ch>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#

me=$(basename "$0")

## defaults
inputDir=""
outputDir=""
prefix=""
inputFile=""
bamFiles=""
threads=1
memory=4
maxMem=4

## Exit status codes (following <sysexits.h>)
EX_OK=0			# successful termination
EX_USAGE=64		# command line usage error
EX_DATAERR=65		# data format error
EX_NOINPUT=66		# cannot open input
EX_NOUSER=67		# addressee unknown
EX_NOHOST=68		# host name unknown
EX_UNAVAILABLE=69	# service unavailable
EX_SOFTWARE=70		# internal software error
EX_OSERR=71		# system error (e.g., can't fork)
EX_OSFILE=72		# critical OS file missing
EX_CANTCREAT=73		# can't create (user) output file
EX_IOERR=74		# input/output error
EX_TEMPFAIL=75		# temp failure; user is invited to retry
EX_PROTOCOL=76		# remote error in protocol
EX_NOPERM=77		# permission denied
EX_CONFIG=78		# configuration error

## helper functions

function die () {
    rc="$1"
    shift
    (echo -n "$me: ERROR: ";
        if [ $# -gt 0 ]; then echo "$@"; else cat; fi) 1>&2
    exit $rc
}

## usage info

usage () {
    cat <<__EOF__
Usage:
  $me [options] [INDIR] [OUTDIR] [OUTPREFIX] [REGIONFILE] [BAMFILE_A] [BAMFILE_B] ...

This script counts the number of reads per genomic region. DO NOT USE IT FOR RNA-SEQ!

Arguments:
INDIR: Directory with the input file.
OUTDIR: Directory in which all output will be store.
OUTPREFIX: Prefix for output. The output file will be named <OUTPREFIX>.counts.txt
REGIONFILE: Name of the file with the genomic regions (bed/saf/gtf).
BAMFILE_X: Bam file(s) for which the reads shall be counted. Full paths must be given!

Details on the REGIONFILE:
Reads are counted with featureCounts which takes either GTF or SAF. SAF is a tab-delimited
format with five columns - here an example:
GeneID	Chr	Start	End	Strand
497097	chr1	3204563	3207049	-
497097	chr1	3411783	3411982	-
497098	chr1	3660633	3661579	-
BED files will be converted to SAF before counting. GTF works, but only entries with
a feature type of "simpleRegion" are considered. I recommend using BED/SAF for this script.

Details on counting:
# 1) Multireads are counted fractionally (1/n with NH-tag:n). # does not work at the moment - it says --fraction requires -M, even though both is specified...
2) Reads are counted if the overlap by at least 30 bases.
3) Only primary (forward read) alignments are counted.
4) If flagged, duplicated reads are ignored.
5) Alignments can be counted multiple times if regions overlap or are very close to each other.
Options:
  -v            Enable verbose logging (no effect)
  -h            Print this help text
  -t		Number of available threads
  -m            Amount of memory to be allocated (per core, in GB)
__EOF__
}

warn () {
  (echo -n "$me: WARNING: ";
      if [ $# -gt 0 ]; then echo "$@"; else cat; fi) 1>&2
}

have_command () {
  type "$1" >/dev/null 2>/dev/null
}

require_command () {
  if ! have_command "$1"; then
    die $EX_UNAVAILABLE "Could not find required command '$1' in system PATH. Aborting."
  fi
}

is_absolute_path () {
    expr match "$1" '/' >/dev/null 2>/dev/null
}

input_exists () {
  echo -n "Checking input file ${1}..."
  if [ -e $1 ]; then
    echo "[ok]"
  else
    echo "[failed]"
    die $EX_NOINPUT "Could not find input file: ${1}."
  fi
}

output_exists () {
  echo -n "Checking output file ${1}... "
  if [ -e $1 ]; then
    echo "[ok]"
  else
    echo "[failed]"
    die $EX_OSFILE "Could not find output file: ${1}."
  fi
}

remove_if_present () {
  if [ -e $1 ]; then
    rm $1
  fi
}

## parse command-line

short_opts='hvt:m:'
long_opts='help,verbose,threads,memory'

getopt -T > /dev/null
rc=$?
if [ "$rc" -eq 4 ]; then
    # GNU getopt
    args=$(getopt --name "$me" --shell sh -l "$long_opts" -o "$short_opts" -- "$@")
    if [ $? -ne 0 ]; then
        die $EX_USAGE "Type '$me --help' to get usage information."
    fi
    # use 'eval' to remove getopt quoting
    eval set -- $args
else
    # old-style getopt, use compatibility syntax
    args=$(getopt "$short_opts" "$@")
    if [ $? -ne 0 ]; then
        die $EX_USAGE "Type '$me --help' to get usage information."
    fi
    set -- $args
fi

while [ $# -gt 0 ]; do
    case "$1" in
	--threads|-t)  shift; threads=$1 ;;
	--memory|-m)   shift; memory=$1 ;;	
        --verbose|-v)  verbose='--verbose' ;;
        --help|-h)     usage; exit 0 ;;
        --)            shift; break ;;
    esac
    shift
done

maxMem=$(($threads * $memory))

## sanity checks

# the required arguments must be present 
if [ $# -lt 5 ]; then
    die $EX_USAGE "Missing required arguments. Type '$me --help' to get usage help."
fi

inputDir=$1
shift
outputDir=$1
shift
prefix=$1
shift
inputFile=$1
shift
bamFiles=("$@")

## main
echo "=== ${me}: Starting at `date '+%Y-%m-%d %H:%M:%S'`"

require_command featureCounts

# checking input
input_exists "${inputDir}/${inputFile}"
for curBam in "${bamFiles[@]}"; do
input_exists $curBam
done

# change BED to SAF
if [[ $inputFile =~ \.bed$ ]]; then
annotationFile="${outputDir}/${inputFile//.bed}.saf"
echo "GeneID\tChr\tStart\tEnd\tStrand" > ${annotationFile}
awk -v OFS='\t' '{if (($1!="#") && ($1!="track") && ($1!="browser")) {print $4,$1,$2,$3,"."}}' ${inputDir}/${inputFile} >> ${annotationFile}
else
annotationFile="${inputDir}/${inputFile}"
fi

# check if GTF or SAF
if [[ $annotationFile =~ \.saf$ ]]; then
formSpecOptions="-F SAF"
else
formSpecOptions="-F GTF -t simpleRegion"
fi

# run featureCounts
bamFiles=$(printf "%s " "${bamFiles[@]}")
bamFiles=${bamFiles::-1}
command="featureCounts -T ${threads} -O --primary --ignoreDup --minOverlap 30 ${formSpecOptions}\
	 -a ${annotationFile} -o ${outputDir}/${prefix}.counts.txt ${bamFiles}"
echo "=== ${me}: Running: ${command}"
eval $command
rc=$?
echo "=== ${me}: featureCounts ended with exit code $rc"

# checking output
output_exists "${outputDir}/${prefix}.counts.txt"
remove_if_present "${outputDir}/${inputFile//.bed}.saf"

## All done.
echo "=== ${me}: Script done at `date '+%Y-%m-%d %H:%M:%S'`."
exit $rc
