#Python pipeline to convert NextSeq BCL files to fastq,
#fix errors in barcodes and trim adapters from fastqs.

#Current issues:
# - Depends on perl scripts (yuck!) in Darren's
#   home directory to run properly.

# - Does not currently allow you to modify
#   default parameters of the perl scripts.

# - Hardcodes trimmomatic as located in Darren's
#   home directory.

import argparse
import os
import subprocess
import glob

parser = argparse.ArgumentParser(description='A program to convert NextSeq BCL files to fastq files for scATAC-seq analysis.')
parser.add_argument('-R','--rundir', help='Run directory containing BCL files',dest='rundir',required=True)
parser.add_argument('-O','--outdir', help='Output directory',dest='outdir',required=True)
parser.add_argument('-P','--prefix',help='Output file prefix',dest='prefix',required=True)
args = parser.parse_args()

def submitter(commander):
        """Submits commands directly to the command line and waits for the process to finish."""
        submiting = subprocess.Popen(commander,shell=True)
        submiting.wait()

if args.outdir[-1] != '/':
	args.outdir = args.outdir + '/'

try:
	os.mkdir(args.outdir)
except OSError:
	print 'Outdir already exists...'

print "Making fastq files..."

print "Converting BCL files..."
fastqer = 'module load modules modules-init modules-gs bcl2fastq/2.16 fastqc/0.10.1; bcl2fastq --runfolder-dir ' + args.rundir + ' -o ' + args.outdir + '/fastq/ --no-lane-splitting'
#submitter(fastqer)

print "Cleaning fastq files..."
fastqdic = {'@':'winner'}
for i,read in enumerate(['R1','R2']):
	curfiles = glob.glob(args.outdir + 'fastq/Undetermined_S0_*' + read + '*')
	if len(curfiles) > 1:
		cater = 'zcat ' + args.outdir + 'fastq/Undetermined_S0_L00*_' + read + '_001.fastq.gz > ' + args.outdir + args.prefix + '.' + str(i+1) + '.temp.fq'
#		submitter(cater)
	else:
		cper = 'cp ' + args.outdir + 'fastq/Undetermined_S0_' + read + '_001.fastq.gz ' + args.outdir + args.prefix + '.' + str(i+1) + '.temp.fq.gz; gunzip ' + args.outdir + args.prefix + '.' + str(i+1) + '.temp.fq.gz'
#		submitter(cper)
#	infile = open(args.outdir + args.prefix + '.' + str(i+1) + '.temp.fq','r')
#	outfile = open(args.outdir + args.prefix + '.' + str(i+1) + '.fq','w')
	z = 1
#	for line in infile:
#		try:
#			fastqdic[line[0]]
#                	barcodeclean = line.strip().split()[1].split(':')[3].replace('+','')
#                	print >> outfile, '@' + barcodeclean + ':' + str(z)
#			z += 1
#		except KeyError:
#                	print >> outfile, line.strip()
#	infile.close()
#	outfile.close()
	cleaner = 'rm ' + args.outdir + args.prefix + '.' + str(i+1) + '.temp.fq; gzip ' + args.outdir + args.prefix + '.' + str(i+1) + '.fq'
#	submitter(cleaner)

print "Fixing barcodes..."

#splitter = 'perl /net/shendure/vol1/home/cusanovi/scatac/NCP_fastq_10bpbarcode_split.pl -1 ' + args.outdir + args.prefix + '.1.fq.gz -2 ' + args.outdir + args.prefix + '.2.fq.gz -O ' + args.outdir + args.prefix + ' -X'
#submitter(splitter)

splitter = 'python /net/shendure/vol1/home/cusanovi/scatac/NCP_fastq_10bpbarcode_split.py -1 ' + args.outdir + args.prefix + '.1.fq -2 ' + args.outdir + args.prefix + '.2.fq -O1 ' + args.outdir + args.prefix + '.split.1.fq -O2 ' + args.outdir + args.prefix + '.split.2.fq -L ' + args.outdir + args.prefix + '.split.log -X -Z'
submitter(splitter)

print "Trimming adapters..."

trimmer = 'java -Xmx1G -jar /net/shendure/vol1/home/cusanovi/bin/Trimmomatic-0.32/trimmomatic-0.32.jar PE '+ args.outdir + args.prefix + '.split.1.fq.gz ' + args.outdir + args.prefix + '.split.2.fq.gz ' + args.outdir + args.prefix + '.split.1.trimmed.paired.fastq.gz ' + args.outdir + args.prefix + '.split.1.trimmed.unpaired.fastq.gz ' + args.outdir + args.prefix + '.split.2.trimmed.paired.fastq.gz ' + args.outdir + args.prefix + '.split.2.trimmed.unpaired.fastq.gz ILLUMINACLIP:/net/shendure/vol1/home/cusanovi/bin/Trimmomatic-0.32/adapters/NexteraPE-PE.fa:2:30:10:1:true MINLEN:20'
submitter(trimmer)

print "Cleaning up..."

cleaner = 'rm ' + args.outdir + args.prefix + '.split.1.fq.gz; rm ' + args.outdir + args.prefix + '.split.2.fq.gz; rm ' + args.outdir + args.prefix + '.split.1.trimmed.unpaired.fastq.gz; rm ' + args.outdir + args.prefix + '.split.2.trimmed.unpaired.fastq.gz;'
#submitter(cleaner)
