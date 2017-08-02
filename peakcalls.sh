module load cython/0.18.0
#TAGS=`find /net/shendure/vol10/projects/mouse_atlas/nobackup/tissues/ -name "*.[0-9].tagAlign.gz" -exec echo {} \;`
TAGS=`find /net/shendure/vol10/projects/human_EBs/nobackup/NextSeq_hEB_July_Sept2016_merged/sort/bam -name "*.[0-9].bam" -exec echo {} \;`
#TAGS=`find $1 -maxdepth 1 -name "*.true.nodups.tagAlign.gz" -exec echo {} \;`
#TAGS=`find /net/shendure/vol10/projects/mouse_atlas/nobackup/tissues/ -name "*.true.nodups.tagAlign.gz" -exec echo {} \;`
for TAG in $TAGS;
do
	echo $TAG
	TAGSHORT=`echo $TAG | sed 's/[.]bam//g' | sed 's/bam/macs/g'`
	#echo $BAMSHORT
	if [ ! -f ${TAGSHORT}_peaks.narrowPeak.gz ]; then
		macs2 callpeak -t $TAG -f BAM -g hs --nomodel --keep-dup all --extsize 200 --shift -100 --call-summits -n $TAGSHORT
		sort -k 8gr,8gr ${TAGSHORT}_peaks.narrowPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | gzip -c > ${TAGSHORT}_peaks.narrowPeak.gz
	else
		echo "${TAGSHORT}_peaks.narrowPeak.gz already exists!"
	fi
done
