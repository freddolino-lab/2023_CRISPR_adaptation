#!/usr/bin/env python3

"""
The script is to take a BED file with negative positions (chromStart and chromEnd) and position exceed the linear representation of a circular genome with single chromosome, and with strand information, and output the corresponding sequences in fasta format, with sequence id containing chromosome position information, strand, and BED score value

USAGE:
python3 apoCas9-05-getFastaFromBed.py <BED_FILENAME> <REF_FASTA_FILENAME>

OUTPUT:
write to standard output.
A FASTA file with entry names formatted as ">chrom:chromStart-chromEnd(strand)[score]" and entry sequences in one lines
Notice that :
1. chromStart, chromEnd values are exactly same as the input BED file (0-based, half-opened). Therefore the coordinates are 0-based, chromStart is and chromeEnd is not included in the region of interest
2. the entry sequences in the output FASTA format are not wrapped to fixed length, i.e. the full squences would be written as one line. This is because the original purpose of this script is for CRISPR Adaptation project extracting sequences from regions of fixed length of 75

UPDATES:
20200426: support multiple entries in reference FASTA file (FASTAREF), check Python3 compatibility
"""

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def parseGenomeFasta(genomeFastaFilename):
    """
    parse FASTA file taking a file name (string) and return a dictionary with chromosome names as keys and chromosome sequences as values
    """
    # chromName = genomeFastaLines[0].rstrip()
    # chromSeq = ''.join(genomeFastaLines[1:]).replace('\n', '')
    # return(chromName, chromSeq)
    genomeSeqDict = dict()
    for record in SeqIO.parse(genomeFastaFilename, 'fasta'):
        genomeSeqDict[record.id] = str(record.seq)
    # print(genomeSeqDict.keys())
    # print([len(genomeSeqDict['NC_017501.1']), len(genomeSeqDict['pYZEJS040'])])
    return(genomeSeqDict)

def reverseComplement(seqString):
    """
    get reverse complement of a DNA sequence string
    """
    seq = Seq(seqString, IUPAC.ambiguous_dna)
    seqRc = seq.reverse_complement()
    return(str(seqRc))

def getSeqRobust(regionStart, regionEnd, strand, chromosomeSeq):
    # BED file positions are 0-based, half-open
    chromosomeSize = len(chromosomeSeq)
    # regionStart < 0
    if (regionStart < 0):
        leftSeq = chromosomeSeq[(chromosomeSize + regionStart):(chromosomeSize + 1)]
        regionStart = 0
    else:
        leftSeq = ''
    # regionEnd >= chromosomeSize
    if (regionEnd >= chromosomeSize):
        rightSeq = chromosomeSeq[0:(regionEnd - chromosomeSize)]
        regionEnd = chromosomeSize
    else:
        rightSeq = ''
    # literal sequence from genomic region
    chromosomeRegionSeq = leftSeq + chromosomeSeq[regionStart:regionEnd] + rightSeq
    if (strand == '+'):
        fastaSeq = chromosomeRegionSeq
    else:
        fastaSeq = reverseComplement(chromosomeRegionSeq)
    return(fastaSeq)

def parseBedToFasta(bedLine, genomeSeq):
    """
    bedLine: a string of an unprocessed line in a BED files
    genomeSeq: a dictionary with chromosome names as keys and chromosome sequences as values
    """
    bedList = bedLine.rstrip().split('\t')
    chrom = bedList[0]
    chromStart = int(bedList[1])
    chromEnd = int(bedList[2])
    score = bedList[4]
    strand = bedList[5]
    fastaName = bedList[3] + '(' + strand + ')[' + score + ']'
    fastaSeq = getSeqRobust(chromStart, chromEnd, strand, genomeSeq[chrom])
    return(fastaName, fastaSeq)

if __name__=='__main__':
    # read inputs
    bedFilename = sys.argv[1]
    genomeFastaFilename = sys.argv[2]
    # parse genomeFastaFile
    genomeSeqDict = parseGenomeFasta(genomeFastaFilename)
    # print(genomeSeqDict.keys())
    # # parse BED file
    with open(bedFilename) as bedFile:
        for line in bedFile:
            # parse BED to Fasta
            fastaName, fastaSeq = parseBedToFasta(line, genomeSeqDict)
            # format Fasta for output
            print('>' + fastaName + '\n' + fastaSeq)
