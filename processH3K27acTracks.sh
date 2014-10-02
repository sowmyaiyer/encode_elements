paste allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.bed processH3K27ac/H3K27ac_percentiles.*.numbers.txt > allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.H3K27ac_percentiles.numbers.txt
awk -f processH3K27ac_for_display.awk allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.H3K27ac_percentiles.numbers.txt > allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.H3K27ac_percentiles.for_display.normalized.bed
$HOME/bedToBigBed -tab -type=bed5+1 -as=display/h3k27ac.as allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.H3K27ac_percentiles.for_display.normalized.bed $HOME/hg19.genome display/allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.H3K27ac_percentiles.for_display.normalized.bb


# Prepare for clustering
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
