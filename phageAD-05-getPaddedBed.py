#!/usr/bin/env python 

import sys

# arguments: inputFileName, upstream5primeLength, downstream5primeLEngth, mode = ['forward', 'reverse']

# input file should be a bedgraph file

# EXAMPLE RUNLINE
## python cz-8-getPaddedSeq.py 8_coverage/170804_4_S4.bedgraph 25 50 forward > example.padded.bed

# this script doesn't deal with cicular genome 20171214

# BUG FOUND: 20180327 the BED file should be 0-based, while SAM file pos is 1-based and genomecov bedgraph is 1-based 

# BUF FIX VERSION: ~/data/cz-Feb2018/cz-04-getPaddedBed.py

# BUG FOUND: 20180421, the BED interval cacluation was shifted

# BUG FOUND: 20191028, bedgraph file lines can have intervals > 1 base

# INPUT is a bedgraph file
# BAM file coordinates are 1-based 
# the coordinates from in the bedgraph files are already converted to 0-based, i.e. the input bedgraph 
# files follows the bedgraph standards
# INPUT as a bedgraph file
# coordiantes are 0-based, half-opened
# http://genome.ucsc.edu/goldenPath/help/bedgraph.html
# OUTPUT as a BED file:
# coordinates are 0-based, half-opened interval, i.e. [chromStart, chromEnd)

# Overall, the coordinate systems used in the order of the pipeline are:
# .bam (including all BAM files in 04-align/): POS at 4th column, 1-based
# .bedgraph (in 05-pad/, as *.forward.bedgraph and *.reverse.bedgraph): 0-based, half-opened
# .bed (in 05-pad/, as *.forward.padded.bed and *.reverse.padded.bed, and the overall *.padded.bed): 0-based, half-opened

inputFileName = sys.argv[1]
upLength = int(sys.argv[2])
downLength = int(sys.argv[3])
mode = sys.argv[4]

def padLength(pos, mode, up, down):
    # half open result 
    if (mode == 'forward'):
        # upstream
        left = pos - up
        # downstream
        right = pos + down
    elif (mode == 'reverse'):
        # upstream
        right = pos + up + 1
        # downstream
        left = pos - down + 1
    else:
        print('ERROR: invalid mode')
        sys.exit()
    return (left, right)

def formatBedLine(paddedRange, fieldList, mode):
    if (mode == 'forward'):
        strand = '+'
    elif (mode == 'reverse'):
        strand = '-'
    else:
        print('ERROR: invalid mode')
        sys.exit()
    # paddedRange is a tuple of 2 elements, left and right position
    # BED: chrom, chromStart, chromEnd, name, score, strand
    name = fieldList[0] + ':' + str(paddedRange[0]) + '-' + str(paddedRange[1])
    score = str(fieldList[3])
    printList = [chrom, str(paddedRange[0]), str(paddedRange[1]), name, score, strand]
    printStr = '\t'.join(printList)
    return(printStr)

with open(inputFileName) as inputFile:
    for inputLine in inputFile:
        # print(inputLine)
        fieldList = inputLine.split()
        if (fieldList[3] != "0"):
            # print(fieldList[1])
            chrom = fieldList[0]
            # BED file is 0-based, when SAM file `pos` and bedtools genomecov results are 1-based
            chromStart = int(fieldList[1])
            chromEnd = int(fieldList[2])
            for pos5prime in range(chromStart, chromEnd):
                # pos5prime = chromStart - 1 
                paddedRange = padLength(pos5prime, up = upLength, down = downLength, mode = mode)
                outLine = formatBedLine(paddedRange, fieldList = fieldList, mode = mode)
                print(outLine)

