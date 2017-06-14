import sys
import pysam

#print len(sys.argv)
if len(sys.argv) != 5:
	sys.exit('Usage: python sc_atac_readcounter.py [Input Sam file] [Input Index table (NA if no table)] [Cell Type 1 BED] [Cell Type 2 BED] [Output file]')

inputbam = sys.argv[1]
indextable = sys.argv[2]
type1bed = sys.argv[3]
outfile = sys.argv[4]


type1desc = type1bed.split(".merged.whitelist.bed")[0].split(".")[-1]

totalct = {}
celltype1ct = {}
descriptor = {}
descdic = {}
if indextable != 'NA':
	descer = open(indextable,'r')
	for line in descer:
		liner = line.strip().split()
		descdic[liner[0]] = liner[1]

print "Counting total reads..."
bamfile = pysam.Samfile(inputbam,'rb')
for read in  bamfile.fetch():
	tagger = read.qname.split(':')[0]
	if 'CTF' in tagger or 'AMBIG' in tagger:
		continue
	try:
		totalct[tagger] += 1
	except KeyError:
		totalct[tagger] = 1
		celltype1ct[tagger] = 0
		try:
			descriptor[tagger] = descdic[tagger]
		except KeyError:
			descriptor[tagger] = 'bkgd'

bamfile.close()

def lister(bedfile):
        currfile = open(bedfile,'r')
        currrecout = [line.strip().split()[0:3] for line in currfile]
        currfile.close()
        return currrecout

def tuper(bedfile):
	currfile = open(bedfile,'r')
	currrecout = ()
	for line in currfile:
        	liner = line.strip().split()
        	currrec = (liner[0],liner[1],liner[2])
        	currrecout = currrecout + (currrec,)
	currfile.close()
	return currrecout

print "Counting celltype-specific reads..."
print "Building " + type1desc + " map..."
rec1tup = lister(type1bed)
bamfile = pysam.Samfile(inputbam,'rb')

def counter(bedtuple,countdic,labeler):
	for rec in bedtuple:
		reads = bamfile.fetch(rec[0], int(rec[1]), int(rec[2]))
		for read in reads:
			readname = read.qname.split(':')[0]
			if 'CTF' in readname or 'AMBIG' in readname:
				continue
			countdic[readname] += 1

print "Counting " + type1desc + " reads..."
counter(rec1tup,celltype1ct,type1desc)

outter = open(outfile,'w')
print >> outter, 'Tag\tTotal\t' + type1desc
for tag in sorted(totalct.keys()):
	print >> outter, tag + '\t' + descriptor[tag] + '\t' + str(totalct[tag]) + "\t" + str(celltype1ct[tag])

outter.close()
