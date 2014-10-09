for distance in {"distal","proximal"}
do
# Collapse peaks with identical coordinates but differing signal values, pick max value to shade box
sort -k1,1 -k2,2n -k4,4 -k5,5 allEncodeDnasePeaks.H3K27ac_only.sorted.merged.${distance}_intersect_allchipseq_peaks.bed | bedtools groupby -i stdin -g 1,2,3,4,5 -c 6 -o max > allEncodeDnasePeaks.H3K27ac_only.sorted.merged.${distance}_intersect_allchipseq_peaks.duplicate_peaks_collapsed.bed
rm tempChipSeq/allEncodeDnasePeaks.H3K27ac_only.sorted.merged.${distance}_intersect_allchipseq_peaks.duplicate_peaks_collapsed.*.bed
awk '{ print $0 >> "tempChipSeq/allEncodeDnasePeaks.H3K27ac_only.sorted.merged.${distance}_intersect_allchipseq_peaks.duplicate_peaks_collapsed."$5".bed"}' allEncodeDnasePeaks.H3K27ac_only.sorted.merged.${distance}_intersect_allchipseq_peaks.duplicate_peaks_collapsed.bed
rm $HOME/public_html/encode_elements/allEncodeDnasePeaks.H3K27ac_only.sorted.merged.${distance}_intersect_allchipseq_peaks.merged_by_tf.bed

# Merge peaks of same TF across different cell types
for f in $HOME/public_html/encode_elements/tempChipSeq/allEncodeDnasePeaks.H3K27ac_only.sorted.merged.${distance}_intersect_allchipseq_peaks.duplicate_peaks_collapsed.*.bed
do
	echo $f
        awk '{ print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$1"_"$2"_"$3 }' $f > $f.temp.txt
        sort -k1,1 -k2,2n  $f.temp.txt | bedtools merge -i stdin -c 4,6,7,5 -o distinct,collapse,collapse,distinct > $f.collapsed.txt
        awk -f pickChIPSeqInterval.awk $f.collapsed.txt > $f.maxScoringPeaks.txt
        cat $f.maxScoringPeaks.txt >> $HOME/public_html/encode_elements/allEncodeDnasePeaks.H3K27ac_only.sorted.merged.${distance}_intersect_allchipseq_peaks.merged_by_tf.bed
done
sort -k1,1 -k2,2n $HOME/public_html/encode_elements/allEncodeDnasePeaks.H3K27ac_only.sorted.merged.${distance}_intersect_allchipseq_peaks.merged_by_tf.bed > $HOME/public_html/encode_elements/allEncodeDnasePeaks.H3K27ac_only.sorted.merged.${distance}_intersect_allchipseq_peaks.merged_by_tf.sorted.bed
# END Merge peaks of same TF across different cell types

# Rescale scores to 0-1000 scale
maxChipSeqScore=`awk -f getMaxScore.awk allEncodeDnasePeaks.H3K27ac_only.sorted.merged.${distance}_intersect_allchipseq_peaks.merged_by_tf.sorted.bed`
minChipSeqScore=`awk -f getMinScore.awk allEncodeDnasePeaks.H3K27ac_only.sorted.merged.${distance}_intersect_allchipseq_peaks.merged_by_tf.sorted.bed`
echo $maxChipSeqScore $minChipSeqScore
awk -vmin=$minChipSeqScore -vmax=$maxChipSeqScore '{printf("%s\t%d\t%d\t%s\t%.0f\t%s\n",$1,$2,$3,$4,(($5-min)/(max-min))*1000, $6) }' allEncodeDnasePeaks.H3K27ac_only.sorted.merged.${distance}_intersect_allchipseq_peaks.merged_by_tf.sorted.bed > allEncodeDnasePeaks.H3K27ac_only.sorted.merged.${distance}_intersect_allchipseq_peaks.merged_by_tf.sorted.normalized.bed
# END Rescale scores to 0-1000 scale

# Convert to bigBed format for display. -as option details additional columns to be displayed on the detail page
~/bedToBigBed -tab -type=bed5+1 -as=display/tf.track.as allEncodeDnasePeaks.H3K27ac_only.sorted.merged.${distance}_intersect_allchipseq_peaks.merged_by_tf.sorted.normalized.bed ~/hg19.genome display/allEncodeDnasePeaks.H3K27ac_only.sorted.merged.${distance}_intersect_allchipseq_peaks.merged_by_tf.sorted.normalized.bb
done

# Requirement to collapse TF peaks in line with DNase peaks - new ChIP-seq file from cluster
max=`awk -f getMaxScore.awk allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal_intersect_allchipseq_peaks.tf_collapsed.bed`
min=`awk -f getMinScore.awk allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal_intersect_allchipseq_peaks.tf_collapsed.bed`
awk -vmin=$min -vmax=$max '{
		split($4,tfs,",")
		name=""
		if (length(tfs) > 5)
		{
			name="multiple"
		} else {
			name=$4
		}
		printf("%s\t%d\t%d\t%s\t%.0f\t%s\n",$1,$2,$3,name,(($5-min)/(max-min))*1000, $6) }
		' allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal_intersect_allchipseq_peaks.tf_collapsed.bed > allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal_intersect_allchipseq_peaks.tf_collapsed.normalized.bed
~/bedToBigBed -tab -type=bed5+1 -as=display/tf.track.2.as allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal_intersect_allchipseq_peaks.tf_collapsed.normalized.bed ~/hg19.genome display/allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal_intersect_allchipseq_peaks.tf_collapsed.normalized.bb

max=`awk -f getMaxScore.awk allEncodeDnasePeaks.H3K27ac_only.sorted.merged.proximal_intersect_allchipseq_peaks.tf_collapsed.bed`
min=`awk -f getMinScore.awk allEncodeDnasePeaks.H3K27ac_only.sorted.merged.proximal_intersect_allchipseq_peaks.tf_collapsed.bed`
awk -vmin=$min -vmax=$max '{
                split($4,tfs,",")
                name=""
                if (length(tfs) > 5)
                {
                        name="multiple"
                } else {
                        name=$4
                }
                printf("%s\t%d\t%d\t%s\t%.0f\t%s\n",$1,$2,$3,name,(($5-min)/(max-min))*1000, $6) }
                ' allEncodeDnasePeaks.H3K27ac_only.sorted.merged.proximal_intersect_allchipseq_peaks.tf_collapsed.bed > allEncodeDnasePeaks.H3K27ac_only.sorted.merged.proximal_intersect_allchipseq_peaks.tf_collapsed.normalized.bed

~/bedToBigBed -tab -type=bed5+1 -as=display/tf.track.2.as allEncodeDnasePeaks.H3K27ac_only.sorted.merged.proximal_intersect_allchipseq_peaks.tf_collapsed.normalized.bed ~/hg19.genome display/allEncodeDnasePeaks.H3K27ac_only.sorted.merged.proximal_intersect_allchipseq_peaks.tf_collapsed.normalized.bb

# END Requirement to collapse TF peaks in line with DNase peaks
