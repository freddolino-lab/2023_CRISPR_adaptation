#!/usr/bin/env python3

"""
take a Fasta file (name) as input and write the reverse complement of sequences to a new Fasta file
USAGE:
python3 cz-07-reverseComplementFasta.py <input.fasta> <output.fasta>
EXAMPLE:
(py3.4.3) diaorch@serenity:~/data/cz-Feb2018$ python3 cz-07-reverseComplementFasta.py 04-cov/czFeb2018_167-9_S25.padded.fasta 07-GATT/czFeb2018_167-9_S25.padded.rc.fasta
"""

import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def convertReverseComplementRecord(iRecord):
    """
    take a SeqIO record as input, return the SeqIO record with 1)sequnce name with (rc) attached to the end; and 2) sequence as the reverse complement of input sequence
    """
    rcRecord = SeqRecord(iRecord.seq.reverse_complement(), id = iRecord.id + '(rc)', name = '', description = '')
    return (rcRecord)

def convertReverseComplementIter(iSeqIter):
    """
    take an iterator of input Fasta seq records as input, return an iterator of output Fasta seq records with 1) (rc) attached to end of sequence name; and 2) reverse complement of input sequences as output sequences
    """
    rcSeqIter = [convertReverseComplementRecord(record) for record in iSeqIter]
    return(rcSeqIter)

def main(iName, oName):
    # adopted from example from:
    # http://biopython.org/wiki/SeqIO
    # 'generator expression'
    inputSeqIterator = SeqIO.parse(iName, 'fasta')
    outputSeqIterator = convertReverseComplementIter(inputSeqIterator)
    SeqIO.write(outputSeqIterator, oName, 'fasta')


if __name__ == '__main__':
    inputFilename = sys.argv[1]
    outputFilename = sys.argv[2]
    main(iName = inputFilename, oName = outputFilename)
