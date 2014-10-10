awk -F" " '{ 
		gsub(",$", "", $12)
		split($12, arr,",")
		for (f = 1; f <= length(arr); f ++)
		{
			printf("%s\t%d\t%d\t%s\t%s\n",$1,$4,$5,$7,arr[f])
		}
	}' gencode.v19.TSS.notlow.gff | sort -k5,5 > gencode.v19.TSS.notlow.transcriptIds.sortedByTranscriptId.bed
awk -F" " '{  if ($1!~/#/ && ($3 == "transcript") ) { print $12"\t"$10"\t"$16"\t"$18"\t"$24"\t"$26 } }' gencode.v19.annotation.gtf | sed 's/;//g' | sed 's/\"//g' | sort | uniq > transcriptId_gene_mapping_sortedByTranscriptId
join gencode.v19.TSS.notlow.transcriptIds.sortedByTranscriptId.bed transcriptId_gene_mapping_sortedByTranscriptId -1 5 -2 1  --check-order |  awk -F" " '{ print $2"\t"$3"\t"$4"\t"$1"\t.\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10 }' | sort -k1,1 -k2,2n > gencode.v19.TSS.notlow.withTranscriptDetails.bed
# gencode.v19.TSS.notlow.withTranscriptDetails.bed is in bed6+ format. Added fields are geneId,status, geneName, transcriptName, level
bedtools window -w 500000 -a allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.bed -b gencode.v19.TSS.notlow.withTranscriptDetails.bed  > allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.with_nearby_genes.bed
split allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.with_nearby_genes.bed -a 3 -d -l 100000 splits/allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.with_nearby_genes.bed
for splitFile in splits/allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.with_nearby_genes.bed[0-9][0-9][0-9]
do
	echo """
	awk '{ printf(\"%s\\t%d\\t%d\\tline_%d\\n\",\$1,(\$2+\$3)/2,(\$2+\$3)/2,NR)}' $splitFile | bedtools slop -b 500 -i stdin -g $HOME/TR/hg19.genome > $splitFile.distal_master_peaks.with_nearby_genes.1000.bed
	""" >  splits/targetGenes.part1.distalPeaks.`basename $splitFile`.bsub
	echo """
	awk -F\"\\t\" '{ print \$9\"\\t\"\$10\"\\t\"\$11\"\\t\"NR}'  $splitFile | bedtools slop -i stdin -b 500 -g $HOME/TR/hg19.genome > $splitFile.nearby_genes_plusminus500bp.bed """ > splits/targetGenes.part2.distalPeaks.`basename $splitFile`.bsub
done
for splitFile in splits/*bed[0-9][0-9][0-9].distal_master_peaks.with_nearby_genes.1000.bed
do
	for cellLine in `cut -f1 cellLinesWithH3K27acAndDnase.txt`
	do
		signal_file=`\ls $HOME/nearline/decorate_dnase/allchipseq_signal/wgEncodeBroadHistone${cellLine}H3k27ac*.bw`
		echo """
			bigWigAverageOverBed $signal_file $splitFile stdout | awk '{ print \$NF}' > splits/`basename $splitFile`.H3K27ac_signal.targetgenes.${cellLine}.txt 
		""" > splits/targetGenes.genes.`basename $splitFile`.${cellLine}.bsub
	done
done
for splitFile in splits/*bed[0-9][0-9][0-9].nearby_genes_plusminus500bp.bed
do
	for cellLine in `cut -f1 cellLinesWithH3K27acAndDnase.txt`
        do
        	signal_file=`\ls $HOME/nearline/decorate_dnase/allchipseq_signal/wgEncodeBroadHistone${cellLine}H3k27ac*.bw`
        	echo """
        	bigWigAverageOverBed $signal_file $splitFile stdout | awk '{ print \$NF}' > splits/`basename $splitFile`.H3K27ac_signal.masterPeaks.${cellLine}.txt
        	""" > splits/targetGenes.distalpeaks.`basename $splitFile`.${cellLine}.bsub
        done
done
for splitFile in `\ls -lv splits/*bed[0-9][0-9][0-9] | awk '{ print $NF}'`
do
	subscript=`echo $splitFile | awk -F"bed" '{ print $2}'`
	echo $subscript
	echo """Rscript correlateH3K27acProfiles_splitFiles.R $subscript """ > splits/targetGeneCorrelations.$subscript.bsub
done
for splitFile in `\ls -lv splits/*bed[0-9][0-9][0-9] | awk '{ print $NF}'`
do
        subscript=`echo $splitFile | awk -F"bed" '{ print $2 }'`
	echo $subscript
	paste $splitFile splits/output_correlations.$subscript > splits/targetGeneCorrelations.$subscript.bed
	awk '{
	split($NF,arr,"_")
        if (arr[1] != "NA" && arr[2] != "NA" && arr[1] > 0.75 && arr[2] < 0.01)
        {
                print
        }
      }' splits/targetGeneCorrelations.$subscript.bed > splits/targetGeneCorrelations.greater_than_0.75_and_p_less_than_0point01.$subscript.bed
done

for splitFile in `\ls -lv splits/targetGeneCorrelations.greater_than_0.75_and_p_less_than_0point01.*.bed  | awk '{ print $NF}'`
do
	echo $splitFile
	cat $splitFile >> targetGeneCorrelations.greater_than_0.75_and_p_less_than_0point01.bed
done
