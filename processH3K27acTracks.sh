paste allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.bed processH3K27ac/H3K27ac_percentiles.*.numbers.txt > allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.H3K27ac_percentiles.numbers.txt
awk -f processH3K27ac_for_display.awk allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.H3K27ac_percentiles.numbers.txt > allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.H3K27ac_percentiles.for_display.normalized.bed
$HOME/bedToBigBed -tab -type=bed5+1 -as=display/h3k27ac.as allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.H3K27ac_percentiles.for_display.normalized.bed $HOME/hg19.genome display/allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.H3K27ac_percentiles.for_display.normalized.bb

# Requirement to group multiple H3K27ac for the same distal peaks
sort -k1,1 -k2,2n allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.H3K27ac_percentiles.for_display.normalized.bed | bedtools groupby -i stdin -g 1,2,3 -c 4,6 -o collapse,collapse  > allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.H3K27ac_collapsed.temp
awk '{
	detail=""
	split($4,cellLines,",")
	split($5,percentiles,",")
	max=-1
	argmax=0
	for (p=1; p <= length(percentiles); p ++)
	{
		if (percentiles[p] > max)
		{
			max=percentiles[p]
			argmax=p
		}	
		detail=sprintf("%s %s(%s)",detail,cellLines[p],percentiles[p])
	}
	printf("%s\t%s\t%s\t%s\t%d\t.\t0\t0\t232,149,16\t%s\t%s(%s)\n",$1,$2,$3,cellLines[argmax],((max-95)/(100-95)*1000),detail,cellLines[argmax],max)
	
	
}' allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.H3K27ac_collapsed.temp > allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.maxH3K27ac.formatted.bed

# Requirement to go back to multiple lines but with color. Arrrggghhh !!!

awk '{printf("%s\t%s\t%s\t%s\t%d\t.\t0\t0\t232,149,16\t%s\n",$1,$2,$3,$4,$5,$6) }' allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.H3K27ac_percentiles.for_display.normalized.bed > allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.H3K27ac_percentiles.for_display.normalized.color.bed

# END Arrrggghhh !!!

$HOME/bedToBigBed -tab -type=bed9+2 -as=sandbox/h3k27ac_withDetails.as allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.maxH3K27ac.formatted.bed $HOME/hg19.genome sandbox/allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.maxH3K27ac.bb
# END Requirement to group multiple H3K27ac for the same distal peaks

# Requiement for name to read "multiple" if there are more than 1 cell line with > 95% H3K27ac
awk -F"\t" '{ 
	split($10,arr," "); 
	if (length(arr) > 1) {
		name="multiple"
	} else {
		name=$4
	}
	printf("%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,name,$5,$6,$7,$8,$9,$10,$11)
}' allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.maxH3K27ac.formatted.bed > allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.maxH3K27ac.formatted.withNameChange.bed
$HOME/bedToBigBed -tab -type=bed9+2 -as=sandbox/h3k27ac_withDetails.as allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.maxH3K27ac.formatted.withNameChange.bed $HOME/hg19.genome sandbox/allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.maxH3K27ac.formatted.withNameChange.bb

# END Requiement for name to read "multiple" if there are more than 1 cell line with > 95% H3K27ac


# Prepare for clustering
rm allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.H3K27ac_percentiles.numbers.for_clustering.txt
head -1 allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.H3K27ac_percentiles.numbers.txt | awk -F"\t" '{
	for (f=9; f <= 22; f ++)
        {
                split($f,arr,"=")
                printf("%s\t",arr[1])
        }
        printf("\n")	
}' > allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.H3K27ac_percentiles.numbers.for_clustering.txt

awk -F"\t" '{
        for (f=9; f <= 22; f ++)
        {
                split($f,arr,"=")
                printf("%s\t",arr[2])
        }
        printf("\n")
}'  allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.H3K27ac_percentiles.numbers.txt >> allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.H3K27ac_percentiles.numbers.for_clustering.txt
