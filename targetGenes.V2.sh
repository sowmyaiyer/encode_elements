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
bedtools window -w 50000 -a allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.bed -b gencode.v19.TSS.notlow.withTranscriptDetails.bed  > allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.with_nearby_genes.bed

awk '{ printf("%s\t%d\t%d\tline_%d\n",$1,($2+$3)/2,($2+$3)/2,NR)}' $HOME/nearline/decorate_dnase/allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.with_nearby_genes.bed | bedtools slop -b 500 -i stdin -g $HOME/TR/hg19.genome > $HOME/nearline/decorate_dnase/distal_master_peaks.with_nearby_genes.1000.bed

awk -F"\t" '{ print $9"\t"$10"\t"$11"\t"NR}' allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.with_nearby_genes.bed | bedtools slop -i stdin -b 500 -g $HOME/TR/hg19.genome > nearby_genes_plusminus500bp.bed
for cellLine in `cut -f1 cellLinesWithH3K27acAndDnase.txt`
do
	signal_file=`\ls $HOME/nearline/decorate_dnase/allchipseq_signal/wgEncodeBroadHistone${cellLine}H3k27ac*.bw`
	echo """
	bigWigAverageOverBed $signal_file $HOME/nearline/decorate_dnase/nearby_genes_plusminus500bp.bed stdout | awk '{ print \$NF}' > targetGenes/H3K27ac_signal.targetgenes.${cellLine}.txt 
	""" > targetGenes.genes.${cellLine}.bsub
	echo """
	bigWigAverageOverBed $signal_file $HOME/nearline/decorate_dnase/distal_master_peaks.with_nearby_genes.1000.bed stdout | awk '{ print \$NF}' > targetGenes/H3K27ac_signal.masterPeaks.${cellLine}.txt 
	""" > targetGenes.distalpeaks.${cellLine}.bsub
done
Rscript correlateH3K27acProfiles.R
split H3K27ac_signal_master_peaks_and_target_genes.txt -a 3 -d -l 100000 splits/H3K27ac_signal_master_peaks_and_target_genes
for f in $HOME/nearline/decorate_dnase/splits/H3K27ac_signal_master_peaks_and_target_genes*
do
        echo """Rscript correlateH3K27acProfiles_part2.R $f $f.targetGeneCorrelations.txt""" > correlateH3K27acProfiles_part2.`basename $f`.bsub
done
if [[ -f targetGeneCorrelations.txt ]]; then
        rm targetGeneCorrelations.txt
fi
for splitFile in `\ls -lv splits/*targetGeneCorrelations.txt | awk '{ print $NF}'`
do
        cat $splitFile >> targetGeneCorrelations.txt
done
paste allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.with_nearby_genes.bed targetGeneCorrelations.txt > targetGeneCorrelations.bed
awk '{ 
	split($NF, arr,"_") 
	if (arr[4] != "NA" && arr[5] != "NA" && arr[4] > 0.7 && arr[5] < 0.01)
	{
		print 
	}
}' targetGeneCorrelations.bed > targetGeneCorrelations.greater_than_0.7_and_p_less_than_0point01.bed
