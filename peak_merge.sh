#!/usr/bin/env bash
zcat /net/shendure/vol10/projects/human_EBs/nobackup/NextSeq_hEB_July_Sept2016_merged/sort/macs/*.[0-9]_peaks.narrowPeak.gz | awk 'BEGIN{OFS="\t"}{print $1,$2+$10,$2+$10+1}' | bedtools slop -i stdin -g /net/shendure/vol10/projects/scATAC/nobackup/genomes/hg19/chromInfo.txt -b 75 | sort -k1,1 -k2,2n | bedtools merge -i stdin | bedtools intersect -v -a stdin -b /net/shendure/vol10/projects/scATAC/nobackup/genomes/annotations/hg19/wgEncodeDacMapabilityConsensusExcludable.bed > /net/shendure/vol10/projects/human_EBs/nobackup/NextSeq_hEB_July_Sept2016_merged/sort/macs/master_combined_summits.whitelist.bed;

bedtools window -w 2500 -a /net/shendure/vol10/projects/human_EBs/nobackup/NextSeq_hEB_July_Sept2016_merged/sort/macs/master_combined_summits.whitelist.bed -b /net/shendure/vol10/projects/scATAC/nobackup/genomes/annotations/hg19/gencode.v19.annotation.transcripts.tss.bed > /net/shendure/vol10/projects/human_EBs/nobackup/NextSeq_hEB_July_Sept2016_merged/sort/macs/master_combined_summits.within.2.5kb.of.tss.whitelist.bed;

zcat /net/shendure/vol10/projects/human_EBs/nobackup/NextSeq_hEB_July_Sept2016_merged/sort/macs/*.[0-9]_peaks.narrowPeak.gz | sort -k1,1 -k2,2n | bedtools merge -i stdin | bedtools intersect -v -a stdin -b /net/shendure/vol10/projects/scATAC/nobackup/genomes/annotations/hg19/wgEncodeDacMapabilityConsensusExcludable.bed > /net/shendure/vol10/projects/human_EBs/nobackup/NextSeq_hEB_July_Sept2016_merged/sort/macs/master_combined_peaks.whitelist.bed;

bedtools window -w 2500 -a /net/shendure/vol10/projects/human_EBs/nobackup/NextSeq_hEB_July_Sept2016_merged/sort/macs/master_combined_peaks.whitelist.bed -b /net/shendure/vol10/projects/scATAC/nobackup/genomes/annotations/hg19/gencode.v19.annotation.transcripts.tss.bed > /net/shendure/vol10/projects/human_EBs/nobackup/NextSeq_hEB_July_Sept2016_merged/sort/macs/master_combined_peaks.within.2.5kb.of.tss.whitelist.bed;




