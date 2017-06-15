import sys
import pysam

#print len(sys.argv)
if len(sys.argv) != 6:
	sys.exit('Usage: python sc_atac_individual_dhs_readcounter.py [Input Bam file] [Input Index table or "NoTags"] [DHS BED] [Output file] [Include sites with no reads? (True/False)]')

inbam = sys.argv[1]
indextable = sys.argv[2]
type1bed = sys.argv[3]
outfile = sys.argv[4]
includezeroes = sys.argv[5]
outbinary = outfile.replace(".txt",".binary.txt")

if indextable != 'NoTags':
	descer = open(indextable,'r')
	cells = [x.strip().split()[0] for x in descer.readlines() if '@' not in x]
	descer.close()
else:
	cells = ['GM12878']

cellsdic = {}
for x,cell in enumerate(cells):
	cellsdic[cell] = x

def lister(bedfile):
	currfile = open(bedfile,'r')
	currrecout = [line.strip().split()[0:3] for line in currfile]
	currfile.close()
	return currrecout

#print "Counting celltype-specific reads..."
print "Building DHS map..."
rec1list = lister(type1bed)
bamfile = pysam.Samfile(inbam,'rb')

def counter(bedtuple,outsfile,binariesout,first=False):
	#tempindex = sorted(cells)
	templen = len(cells)
	if first:
		print >> outsfile, "chr\tstart\tend\tannot\t" + "\t".join(cells)
		print >> binariesout, "chr\tstart\tend\tannot\t" + "\t".join(cells)
	#recct = 0
	if includezeroes:
		for rec in bedtuple:
			#recct += 1
			recname = rec[0] + "_" + rec[1] + "_" + rec[2]
			currcounts = [0]*templen
			currbinary = ["0"]*templen
			reads = bamfile.fetch(rec[0], int(rec[1]), int(rec[2]))
			for read in reads:
				readname = read.qname.split(':')[0]
				#if 'CTF' in readname or 'AMBIG' in readname:
				#	continue
				#if readname in cells:
				#	currcounts[cells.index(readname)] += 1
				#	currbinary[cells.index(readname)] = 1
				try:
					currcounts[cellsdic[readname]] += 1
					currbinary[cellsdic[readname]] = "1"
				except KeyError:
					pass
				#if indextable == 'NoTags':
				#	currcounts[0] += 1
			#if includezeroes:
			print >> outsfile, "\t".join(rec[0:3]) + "\t" + recname + "\t" + "\t".join([str(x) for x in currcounts])
			#print >> binariesout, "\t".join(rec[0:3]) + "\t" + recname + "\t" + "\t".join([str(x) for x in currbinary])
			print >> binariesout, "\t".join(rec[0:3]) + "\t" + recname + "\t" + "\t".join(currbinary)

	else:
		for rec in bedtuple:
			#recct += 1
			recname = rec[0] + "_" + rec[1] + "_" + rec[2]
			currcounts = [0]*templen
			currbinary = [0]*templen
			reads = bamfile.fetch(rec[0], int(rec[1]), int(rec[2]))
			for read in reads:
				readname = read.qname.split(':')[0]
				#if 'CTF' in readname or 'AMBIG' in readname:
				#	continue
				#if readname in cells:
				#	currcounts[cells.index(readname)] += 1
				#	currbinary[cells.index(readname)] = 1
				try:
					currcounts[cellsdic[readname]] += 1
					currbinary[cellsdic[readname]] = "1"
				except KeyError:
					pass
			if sum(currbinary) > 0:
			#if sum(currcounts) > 9:
				print >> outsfile, "\t".join(rec[0:3]) + "\t" + recname + "\t" + "\t".join([str(x) for x in currcounts])
				#print >> binariesout, "\t".join(rec[0:3]) + "\t" + recname + "\t" + "\t".join([str(x) for x in currbinary])
				print >> binariesout, "\t".join(rec[0:3]) + "\t" + recname + "\t" + "\t".join(currbinary)


dhsmat = open(outfile,'w')
dhsbinary = open(outbinary,'w')
print "Counting DHS reads..."
counter(rec1list,dhsmat,dhsbinary,first=True)
dhsmat.close()
