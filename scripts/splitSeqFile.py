#!/usr/bin/env python

from Bio import SeqIO
import sys, os

def usage():
    print """
    
    splitSeqFile.py <infile> <inSeqFormat> <outSeqFormat> <splitSize>
    
    """

def printSeq(seq):
    # all possible adjustments to sequence before printing
    # eg Prodigal puts a mark at the end of translated sequences
    seq = seq.replace("*", "")
    return seq

if len(sys.argv) != 5:
    usage()
    sys.exit("Insufficient arguments provided.")

acceptedFormats = ['fasta', 'fastq']

infile = sys.argv[1]
intype = sys.argv[2]
outtype = sys.argv[3]
splitSize = int(sys.argv[4])

if intype not in acceptedFormats or outtype not in acceptedFormats:
    usage()
    sys.ext("Infile or outfile format provided is not legal. Exit.") 

inhandle = SeqIO.parse(infile, intype)

FILE, EXT = os.path.splitext(infile)
seqCount = 0
fileCount = 0
outFile = "%s.split%d.%s%s" % (FILE, splitSize, str(fileCount).zfill(5), EXT)
print outFile
outHandle = open(outFile, "w")
for seq in inhandle:
    if outtype == "fasta":
        outHandle.write(">%s\n%s\n" % (seq.id, printSeq(str(seq.seq))))
    elif outtype == "fastq":
        outHandle.write("@%s\n%s\n+\n%s\n" % (seq.id, printSeq(str(seq.seq)), ''.join([chr(i+33) for i in seq.letter_annotations['phred_quality']]))) 
    seqCount += 1
    if seqCount%splitSize == 0:
        outHandle.close()
        fileCount += 1
        outFile = "%s.split%d.%s%s" % (FILE, splitSize, str(fileCount).zfill(5), EXT)
        print outFile
        outHandle = open(outFile, "w")

outHandle.close()
