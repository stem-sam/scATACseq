#!/usr/bin/env python
import argparse
#import os
import subprocess
#import glob
import sys
sys.path.append('/net/shendure/vol1/home/cusanovi/bin/Levenshtein/')
import Levenshtein
import gzip
#import io
import cStringIO
io_method = cStringIO.StringIO

parser = argparse.ArgumentParser(description='A program to fix erroneous barcodes in scATAC data.')
parser.add_argument('-1','--Read1', help='Input fastq file 1',dest='input1',required=True)
parser.add_argument('-2','--Read2', help='Input fastq file 2',dest='input2',required=True)
parser.add_argument('-O1','--output1', help='Output fastq file 1',dest='output1',required=True)
parser.add_argument('-O2','--output2', help='Output fastq file 2',dest='output2',required=True)
parser.add_argument('-L','--log', help='Output log file',dest='logfile',required=True)
parser.add_argument('-X','--nextseq',help='NextSeq run indicator',dest='nextseq',action="store_true")
parser.add_argument('-Z','--gzip',help='Gzip indicator',dest='gzip',action="store_true")
args = parser.parse_args()

def submitter(commander):
        """Submits commands directly to the command line and waits for the process to finish."""
        submiting = subprocess.Popen(commander,shell=True)
        submiting.wait()

def editcheck(barc,reflist):
	try:
		reflist[barc]
		eddist = '0'
	except KeyError:
		winner = '_CTF' + '_'*(len(barc)-4)
		winner_ed = 10
		runnerup_ed = 10
		for barcode in reflist.keys():
			curred = Levenshtein.distance(barc,barcode)
			if curred <= winner_ed:
				runnerup_ed = winner_ed
				winner = barcode
				winner_ed = curred
		if winner_ed > 3:
			winner = '_CTF' + '_'*(len(barc)-4)
		if runnerup_ed - winner_ed < 2:
			winner = '_AMBIG' + '_'*(len(barc)-6)
		barc = winner
		eddist = str(winner_ed)
	return(barc,eddist)

complements = {'A':'T', 'T':'A', 'G':'C', 'C':'G','N':'N'}
def reverse_complement(x):
	xrev = x[::-1]
	xrevcomp = ''.join([complements[z] for z in xrev])
	return xrevcomp

nex_i7 = {"ATTACTCG":"","TCCGGAGA":"","CGCTCATT":"","GAGATTCC":"","ATTCAGAA":"","GAATTCGT":"","CTGAAGCT":"","TAATGCGC":"","CGGCTATG":"","TCCGCGAA":"","TCTCGCGC":"","AGCGATAG":""}
pcr_i7 = {"TCGGATTCGG":"","GCGGCTGCGG":"","AGATTACGTT":"","CTAACTAGGT":"","CATAGCGACC":"","CCGCTAAGAG":"","ATGGAACGAA":"","GCGTTCCGTT":"","GGTTATCGAA":"","GCATCGTATG":"","AATACGATAA":"","TTCCGTCGAC":"","TCCGGCTTAT":"","ACCAGGCGCA":"","AGAGGAGAAT":"","GTACTCCTAT":"","GCTAACGGAT":"","AGTTGAATCA":"","TGATTAGGTA":"","TCGTAGCATC":"","TCTTGAGGTT":"","AGGTCAGCTT":"","TATTAGACTT":"","CTCAATTAGT":"","TCGCCGCCGG":"","CCGTATGATT":"","AACGCGCAGA":"","CTCGTCGTAG":"","CTAATTGCGA":"","CGCGGCCATA":"","AATATTACTT":"","ATTGGCAGAT":"","ATGGCGCCTG":"","ATAAGGACTC":"","TAGTAAGCCG":"","ATTATGCAAG":"","TTGGCAAGCC":"","TTGATTGGCG":"","GCATATGAGC":"","GAACTCGACT":"","CTAGCCAGCC":"","TGCGACCTCT":"","ATTCTTAGCT":"","TTGATACGAT":"","TATAATAGTT":"","TTGCCGTAGG":"","AGACCATATC":"","TTGGTAAGGA":"","CAGCTAGCGG":"","CTAAGCCTTG":"","CGTTACCGCT":"","GACTGGACCA":"","GCAAGACCGT":"","TCAATCTCCT":"","ATACCTCGAC":"","TAGAGGCGTT":"","TAGGTAACTT":"","TTCGAATATT":"","TGGACGACTA":"","GTAGGCTGCA":"","GTAGGATAAG":"","CGTCGAGCGC":"","ACTATTCATT":"","TTGCTTAGAT":"","CGAATGGAGC":"","CTATATAGCC":"","CTACTAATAA":"","TGGTTGCCGT":"","TCCTCTGCCG":"","GATTCTTGAA":"","GTAGCAGCTA":"","CCTCAGCTCC":"","AAGTAGCTCA":"","TATTGCTGGA":"","CCAGATACGG":"","AACGAATTCG":"","CGCTTATCGT":"","AAGTACGCGA":"","GATCTTCGCA":"","TCTTAGCCTG":"","TTATTGAGGC":"","TTGCGAGCAT":"","GCTTGAAGAG":"","AGTCCGCTGC":"","TAAGTCCTGA":"","AGTTCTCATG":"","CAGACTAAGG":"","TCTATCGCTG":"","GCGCTATGGT":"","CATTATTATT":"","AGCCGTAGTT":"","TGATATTGCG":"","ACGGCGTTAA":"","GGCTTACTCC":"","GCGCGTTCAT":"","GAGCGCGATG":""}
pcr_i5 = {"CTCCATCGAG":"","TTGGTAGTCG":"","GGCCGTCAAC":"","CCTAGACGAG":"","TCGTTAGAGC":"","CGTTCTATCA":"","CGGAATCTAA":"","ATGACTGATC":"","TCAATATCGA":"","GTAGACCTGG":"","TTATGACCAA":"","TTGGTCCGTT":"","GGTACGTTAA":"","CAATGAGTCC":"","GATGCAGTTC":"","CCATCGTTCC":"","TTGAGAGAGT":"","ACTGAGCGAC":"","TGAGGAATCA":"","CCTCCGACGG":"","CATTGACGCT":"","TCGTCCTTCG":"","TGATACTCAA":"","TTCTACCTCA":"","TCGTCGGAAC":"","ATCGAGATGA":"","TAGACTAGTC":"","GTCGAAGCAG":"","AGGCGCTAGG":"","AGATGCAACT":"","AAGCCTACGA":"","GTAGGCAATT":"","GGAGGCGGCG":"","CCAGTACTTG":"","GGTCTCGCCG":"","GGCGGAGGTC":"","TAGTTCTAGA":"","TTGGAGTTAG":"","AGATCTTGGT":"","GTAATGATCG":"","CAGAGAGGTC":"","TTAATTAGCC":"","CTCTAACTCG":"","TACGATCATC":"","AGGCGAGAGC":"","TCAAGATAGT":"","TAATTGACCT":"","CAGCCGGCTT":"","AGAACCGGAG":"","GAGATGCATG":"","GATTACCGGA":"","TCGTAACGGT":"","TGGCGACGGA":"","AGTCATAGCC":"","GTCAAGTCCA":"","ATTCGGAAGT":"","GTCGGTAGTT":"","AGGACGGACG":"","CTCCTGGACC":"","TAGCCTCGTT":"","GGTTGAACGT":"","AGGTCCTCGT":"","GGAAGTTATA":"","TGGTAATCCT":"","AAGCTAGGTT":"","TCCGCGGACT":"","TGCGGATAGT":"","TGGCAGCTCG":"","TGCTACGGTC":"","GCGCAATGAC":"","CTTAATCTTG":"","GGAGTTGCGT":"","ACTCGTATCA":"","GGTAATAATG":"","TCCTTATAGA":"","CCGACTCCAA":"","GCCAAGCTTG":"","CATATCCTAT":"","ACCTACGCCA":"","GGAATTCAGT":"","TGGCGTAGAA":"","ATTGCGGCCA":"","TTCAGCTTGG":"","CCATCTGGCA":"","CTTATAAGTT":"","GATTAGATGA":"","TATAGGATCT":"","AGCTTATAGG":"","GTCTGCAATC":"","CGCCTCTTAT":"","GTTGGATCTT":"","GCGATTGCAG":"","TGCCAGTTGC":"","CTTAGGTATC":"","GAGACCTACC":"","ATTGACCGAG":""}
nex_i5 = {"TATAGCCT":"","ATAGAGGC":"","CCTATCCT":"","GGCTCTGA":"","AGGCGAAG":"","TAATCTTA":"","CAGGACGT":"","GTACTGAC":""}

FULL_BARC = {}
for i in nex_i7.keys():
	for j in pcr_i7.keys():
		for k in pcr_i5.keys():
			for l in nex_i5.keys():
				currbarc = i + j + k + l
				FULL_BARC[currbarc] = ""

totreads = 0
exactmatch = 0
editmatch = 0
failed = 0
readcount = 0
prevouts = 0

#infasta1 = open(args.input1,'r')
#infasta2 = open(args.input2,'r')
outfasta1 = open(args.output1,'w')
outfasta2 = open(args.output2,'w')
#if1 = gzip.open(args.input1,'rb')
#infasta1 = io.BufferedReader(if1)
#if2 = gzip.open(args.input2,'rb')
#infasta2 = io.BufferedReader(if2)
p1 = subprocess.Popen(["zcat", args.input1], stdout = subprocess.PIPE)
infasta1 = io_method(p1.communicate()[0])
p2 = subprocess.Popen(["zcat", args.input2], stdout = subprocess.PIPE) 
infasta2 = io_method(p2.communicate()[0])
#infasta1 = gzip.open(args.input1,'rb')
#infasta2 = gzip.open(args.input2,'rb')
#outfasta1 = gzip.open(args.output1,'wb')
#outfasta2 = gzip.open(args.output2,'wb')
#p3 = subprocess.Popen("gzip -c > " + args.output1, shell=True, stdin=subprocess.PIPE)
#p4 = subprocess.Popen("gzip -c > " + args.output2, shell=True, stdin=subprocess.PIPE)
reads1 = []
reads2 = []
for line in infasta1:
	seqbarc = line.strip('@').split(':')[0]
	barcodes = line.strip().split()[1].split(":")[3]
	seqbarc = "".join(barcodes.split('+'))
	if args.nextseq:
		#barcodes = line.strip().split()[1].split(":")[3]
		#seqbarc = "".join(barcodes.split('+'))
		seqbarca = seqbarc[0:18]
		seqbarcb = reverse_complement(seqbarc[18:36])
		seqbarc = seqbarca + seqbarcb
	read1 = infasta1.next().strip()
	plus = infasta1.next()
	qual1 = infasta1.next().strip()
	dumpbarc = infasta2.readline()
        read2 = infasta2.readline().strip()
        plus = infasta2.readline()
        qual2 = infasta2.readline().strip()
	try:
		FULL_BARC[seqbarc]
		#print >> outfasta1, '@' + seqbarc + ':' + str(totreads) + '#0000/1'
		#print >> outfasta1, read1 + '\n+\n' + qual1
		#print >> outfasta2, '@' + seqbarc + ':' + str(totreads) + '#0000/2'
		#print >> outfasta2, read2 + '\n+\n' + qual2
                #outfasta1.write('@' + seqbarc + ':' + str(totreads) + '#0000/1\n' + read1 + '\n+\n' + qual1 + '\n')
                #outfasta2.write('@' + seqbarc + ':' + str(totreads) + '#0000/2\n' + read2 + '\n+\n' + qual2 + '\n')
		reads1.append('@' + seqbarc + ':' + str(totreads) + '#0000/1\n' + read1 + '\n+\n' + qual1 + '\n')
		reads2.append('@' + seqbarc + ':' + str(totreads) + '#0000/2\n' + read2 + '\n+\n' + qual2 + '\n')
		readcount += 1
		totreads += 1
		exactmatch += 1
	except KeyError:
		b1 = seqbarc[0:8]
		b2 = seqbarc[8:18]
		b3 = seqbarc[18:28]
		b4 = seqbarc[28:36]
		b1ed = editcheck(b1,nex_i7)
		b2ed = editcheck(b2,pcr_i7)
		b3ed = editcheck(b3,pcr_i5)
		b4ed = editcheck(b4,nex_i5)
		seqbarc = b1ed[0] + b2ed[0] + b3ed[0] + b4ed[0]
		seqed = b1ed[1] + b2ed[1] + b3ed[1] + b4ed[1]
		#print >> outfasta1, '@' + seqbarc + ':' + str(totreads) + '#' + seqed + '/1'
                #print >> outfasta1, read1 + '\n' + '+\n' + qual1
                #print >> outfasta2, '@' + seqbarc + ':' + str(totreads) + '#' + seqed + '/2'
                #print >> outfasta2, read2 + '\n' + '+\n' + qual2
                #outfasta1.write('@' + seqbarc + ':' + str(totreads) + '#0000/1\n' + read1 + '\n+\n' + qual1 + '\n')
                #outfasta2.write('@' + seqbarc + ':' + str(totreads) + '#0000/2\n' + read2 + '\n+\n' + qual2 + '\n')
		reads1.append('@' + seqbarc + ':' + str(totreads) + '#0000/1\n' + read1 + '\n+\n' + qual1 + '\n')
                reads2.append('@' + seqbarc + ':' + str(totreads) + '#0000/2\n' + read2 + '\n+\n' + qual2 + '\n')
                totreads += 1
		readcount += 1
		if 'CTF' in seqbarc or 'AMBIG' in seqbarc:
                	failed += 1
		else:
			editmatch += 1
	if readcount == 250000:
		#if prevouts == 0:
		#	p3 = subprocess.Popen("gzip -c > " + args.output1, shell=True, stdin=subprocess.PIPE)
		#	p4 = subprocess.Popen("gzip -c > " + args.output2, shell=True, stdin=subprocess.PIPE)	
		#p3.stdin.write("\n".join(reads1) + "\n")
		#p4.stdin.write("\n".join(reads2) + "\n") 			
		outfasta1.writelines(reads1)
		outfasta2.writelines(reads2)
		readcount = 0
		reads1 = []
		reads2 = []

if readcount > 0:
	#p3.stdin.write("\n".join(reads1) + "\n")
        #p4.stdin.write("\n".join(reads2) + "\n")
	outfasta1.writelines(reads1)
	outfasta2.writelines(reads2)

#infasta1.close()
#infasta2.close()
#if1.close()
#if2.close()
outfasta1.close()
outfasta2.close()
#p3.communicate()
#p4.communicate()


logout = open(args.logfile,'w')
print >> logout, 'total=' + str(totreads) + '\texact=' + str(round(float(exactmatch)/totreads,2)) + '\tby_ed=' + str(round(float(editmatch)/totreads,2)) + '\tfail=' + str(round(float(failed)/totreads,2))
logout.close()

if args.gzip:
	zipper = 'gzip ' + args.output1 + '; gzip ' + args.output2 + ';'
	submitter(zipper)

