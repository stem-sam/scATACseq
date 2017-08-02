TAGS=`find /net/shendure/vol10/projects/human_EBs/nobackup/NextSeq_hEB_July_Sept2016_merged/matrix/ -name "*.true.nodups_peaks_merged_clusters.sigopen*.great_enrichments.txt" -exec echo {} \;`
for TAG in $TAGS;
do
	echo $TAG
	TAGSHORT=`echo $TAG | sed 's/great_enrichments[.]txt/great_enrichments.clean.txt/g'`
	#echo $BAMSHORT
	if [ ! -f ${TAGSHORT}.z ]; then
		sed -i 's/&amp;#xef;/i/g' $TAG;
		sed -i "s/'/prime/g" $TAG;
		Rscript /net/shendure/vol10/projects/mouse_development/scripts/sc_atac_great_cleaner.R $TAG $TAGSHORT;
	else
		echo "${TAGSHORT} already exists!"
	fi
done
