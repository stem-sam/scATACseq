#Python script to count duplicate reads for each cell
#in a single cell ATAC-seq experiment.

import argparse
import os
import subprocess
import pysam

def submitter(commander):
        """Submits commands directly to the command line and waits for the process to finish."""
        submiting = subprocess.Popen(commander,shell=True)
        submiting.wait()

parser = argparse.ArgumentParser(description='A program to count duplicate reads for each cell for scATAC-seq analysis.')
parser.add_argument('-B','--inbam', help='Input bam file (not deduplicated)',dest='inbam',required=True)
parser.add_argument('-D','--indedup', help='Input deduplicated bam file',dest='indedup',required=True)
parser.add_argument('-I','--index',help='Index table',dest='indexfile',required=True)
parser.add_argument('-O','--output',help='Output file name',dest='outfile',required=True)
args = parser.parse_args()

print "Building index table..."

indexdic = {}
indices = open(args.indexfile,'r')
for line in indices:
	liner = line.strip().split()
	indexdic[liner[0]] = [0,0]

print "Counting mapped reads..."

bamfile = pysam.Samfile(args.inbam,'rb')
for read in bamfile.fetch():
	try:
		cellname = read.qname.split(':')[0]
		indexdic[cellname][0] += 1
	except KeyError:
		continue

print "Counting deduplicated reads..."

dedupfile = pysam.Samfile(args.indedup,'rb')
for read in dedupfile.fetch():
        try:
                cellname = read.qname.split(':')[0]
                indexdic[cellname][1] += 1
        except KeyError:
                continue

print "Writing results..."

outter = open(args.outfile,'w')
print >> outter, "Barcode\tMapped\tDeduplicated"
for barcode in indexdic.keys():
	if indexdic[barcode][0] > 0:
		print >> outter, barcode + "\t" + str(indexdic[barcode][0]) + "\t" + str(indexdic[barcode][1])

outter.close()
