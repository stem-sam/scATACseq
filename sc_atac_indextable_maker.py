#Script to generate a table of indexes and sample names
#For each sample name, input the nextera i7 index numbers,
#pcr i7 index numbers, pcr i5 index numbers and nextera i5
#index numbers in the format nex_i7:pcr_i7:pcr_i5:nex_i5.
#To assign multiple indexes, separate numbers between colons
#with commas (e.g. '1,3,5').  To input a range, separate
#start and end with '-' (e.g. '1-6').
import sys

indices = sys.argv[1]
name = sys.argv[2]
indexsplit = indices.strip().split(':')

def indexsplitter(indexrange):
	if len(indexrange) < 3:
		indexout = [int(indexrange)-1]
	if "-" in indexrange:
		indexout = range(int(indexrange.split("-")[0])-1,int(indexrange.split("-")[1]))
	if "," in indexrange:
		indexout = [int(x)-1 for x in indexrange.split(",")]
	return indexout

nexi7_indices = indexsplitter(indexsplit[0])
pcri7_indices = indexsplitter(indexsplit[1])
pcri5_indices = indexsplitter(indexsplit[2])
nexi5_indices = indexsplitter(indexsplit[3])

nex_i7 = ["ATTACTCG","TCCGGAGA","CGCTCATT","GAGATTCC","ATTCAGAA","GAATTCGT","CTGAAGCT","TAATGCGC","CGGCTATG","TCCGCGAA","TCTCGCGC","AGCGATAG"]
#Illumina P7 barcodes
#pcr_i7 = ["TAAGGCGA","CGTACTAG","AGGCAGAA","TCCTGAGC","GGACTCCT","TAGGCATG","CTCTCTAC","CAGAGAGG","GCTACGCT","CGAGGCTG","AAGAGGCA","GTAGAGGA","TGCTTCAG","CCTAAGAC","CGATCAGT","TGCAGCTA","TCGACGTC","AGCAGGAG","GATGCCGT","TGAGGACT","GAGGAGAG","CCTGCAGA","AGCCGAGT","TCGTCCGA"]
#Martin P7 barcodes
pcr_i7 = ["TCGGATTCGG","GCGGCTGCGG","AGATTACGTT","CTAACTAGGT","CATAGCGACC","CCGCTAAGAG","ATGGAACGAA","GCGTTCCGTT","GGTTATCGAA","GCATCGTATG","AATACGATAA","TTCCGTCGAC","TCCGGCTTAT","ACCAGGCGCA","AGAGGAGAAT","GTACTCCTAT","GCTAACGGAT","AGTTGAATCA","TGATTAGGTA","TCGTAGCATC","TCTTGAGGTT","AGGTCAGCTT","TATTAGACTT","CTCAATTAGT","TCGCCGCCGG","CCGTATGATT","AACGCGCAGA","CTCGTCGTAG","CTAATTGCGA","CGCGGCCATA","AATATTACTT","ATTGGCAGAT","ATGGCGCCTG","ATAAGGACTC","TAGTAAGCCG","ATTATGCAAG","TTGGCAAGCC","TTGATTGGCG","GCATATGAGC","GAACTCGACT","CTAGCCAGCC","TGCGACCTCT","ATTCTTAGCT","TTGATACGAT","TATAATAGTT","TTGCCGTAGG","AGACCATATC","TTGGTAAGGA","CAGCTAGCGG","CTAAGCCTTG","CGTTACCGCT","GACTGGACCA","GCAAGACCGT","TCAATCTCCT","ATACCTCGAC","TAGAGGCGTT","TAGGTAACTT","TTCGAATATT","TGGACGACTA","GTAGGCTGCA","GTAGGATAAG","CGTCGAGCGC","ACTATTCATT","TTGCTTAGAT","CGAATGGAGC","CTATATAGCC","CTACTAATAA","TGGTTGCCGT","TCCTCTGCCG","GATTCTTGAA","GTAGCAGCTA","CCTCAGCTCC","AAGTAGCTCA","TATTGCTGGA","CCAGATACGG","AACGAATTCG","CGCTTATCGT","AAGTACGCGA","GATCTTCGCA","TCTTAGCCTG","TTATTGAGGC","TTGCGAGCAT","GCTTGAAGAG","AGTCCGCTGC","TAAGTCCTGA","AGTTCTCATG","CAGACTAAGG","TCTATCGCTG","GCGCTATGGT","CATTATTATT","AGCCGTAGTT","TGATATTGCG","ACGGCGTTAA","GGCTTACTCC","GCGCGTTCAT","GAGCGCGATG"]
#Illumina P5 barcodes
#pcr_i5 = ["TAGATCGC","CTCTCTAT","TATCCTCT","AGAGTAGA","GTAAGGAG","ACTGCATA","AAGGAGTA","CTAAGCCT","GAGCTACT","ACTAATGA","TGAGGCAT","AGCTCTGA","GCTCCTAC","TCCGTAGG","AGTACGAC","ATAGTCTT"]
#Martin P5 barcodes
pcr_i5 = ["CTCCATCGAG","TTGGTAGTCG","GGCCGTCAAC","CCTAGACGAG","TCGTTAGAGC","CGTTCTATCA","CGGAATCTAA","ATGACTGATC","TCAATATCGA","GTAGACCTGG","TTATGACCAA","TTGGTCCGTT","GGTACGTTAA","CAATGAGTCC","GATGCAGTTC","CCATCGTTCC","TTGAGAGAGT","ACTGAGCGAC","TGAGGAATCA","CCTCCGACGG","CATTGACGCT","TCGTCCTTCG","TGATACTCAA","TTCTACCTCA","TCGTCGGAAC","ATCGAGATGA","TAGACTAGTC","GTCGAAGCAG","AGGCGCTAGG","AGATGCAACT","AAGCCTACGA","GTAGGCAATT","GGAGGCGGCG","CCAGTACTTG","GGTCTCGCCG","GGCGGAGGTC","TAGTTCTAGA","TTGGAGTTAG","AGATCTTGGT","GTAATGATCG","CAGAGAGGTC","TTAATTAGCC","CTCTAACTCG","TACGATCATC","AGGCGAGAGC","TCAAGATAGT","TAATTGACCT","CAGCCGGCTT","AGAACCGGAG","GAGATGCATG","GATTACCGGA","TCGTAACGGT","TGGCGACGGA","AGTCATAGCC","GTCAAGTCCA","ATTCGGAAGT","GTCGGTAGTT","AGGACGGACG","CTCCTGGACC","TAGCCTCGTT","GGTTGAACGT","AGGTCCTCGT","GGAAGTTATA","TGGTAATCCT","AAGCTAGGTT","TCCGCGGACT","TGCGGATAGT","TGGCAGCTCG","TGCTACGGTC","GCGCAATGAC","CTTAATCTTG","GGAGTTGCGT","ACTCGTATCA","GGTAATAATG","TCCTTATAGA","CCGACTCCAA","GCCAAGCTTG","CATATCCTAT","ACCTACGCCA","GGAATTCAGT","TGGCGTAGAA","ATTGCGGCCA","TTCAGCTTGG","CCATCTGGCA","CTTATAAGTT","GATTAGATGA","TATAGGATCT","AGCTTATAGG","GTCTGCAATC","CGCCTCTTAT","GTTGGATCTT","GCGATTGCAG","TGCCAGTTGC","CTTAGGTATC","GAGACCTACC","ATTGACCGAG"]
nex_i5 = ["TATAGCCT","ATAGAGGC","CCTATCCT","GGCTCTGA","AGGCGAAG","TAATCTTA","CAGGACGT","GTACTGAC"]

for nexi7_id in nexi7_indices:
	for pcri7_id in pcri7_indices:
		for pcri5_id in pcri5_indices:
			for nexi5_id in nexi5_indices:
				print nex_i7[nexi7_id] + pcr_i7[pcri7_id] + pcr_i5[pcri5_id] + nex_i5[nexi5_id] + '\t' + name

