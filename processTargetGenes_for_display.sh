awk -f targetGenes.awk targetGeneCorrelations.greater_than_0.7_and_p_less_than_0point01.bed > sandbox/targetGenes.bed
sort -k1,1 -k2,2n sandbox/targetGenes.bed > sandbox/targetGenes.sorted.bed
$HOME/bedToBigBed -tab -type=bed12+6 -as=sandbox/targetGenes.as sandbox/targetGenes.sorted.bed ~/hg19.genome sandbox/targetGenes.bb

# Requirement to collapse multiple lines to one with highest correlation
# head targetGeneCorrelations.greater_than_0.7_and_p_less_than_0point01.bed
#chr1	87640	87790	Helas3	10	Helas3	10	1	chr1	52473	52473	ENST00000606857.1	.	+	ENSG00000268020.2	KNOWN	OR4G4P	OR4G4P-001	2	Gm12878,H1hesc,Helas3,Hepg2,Hmec,Hsmm,Hsmmt,Huvec,K562,Nha,Nhdfad,Nhek,Nhlf,Osteobl_0.15996,0.39951,0.30516,0,0,0.10347,0,0.11137,0.06225,0,0.228,0,0,0.07799_0,0.12723,0.12758,0,0,0,0,0.11514,0,0,0.044,0,0,0_0.78_0.000945395382869261
# chr1	125340	125490	Helas3	9	Helas3	9	1	chr1	129173	129173	ENST00000471248.1	.	-	ENSG00000238009.2	NOVEL	RP11-34P13.7	RP11-34P13.7-002	2	Gm12878,H1hesc,Helas3,Hepg2,Hmec,Hsmm,Hsmmt,Huvec,K562,Nha,Nhdfad,Nhek,Nhlf,Osteobl_0,0,0,0,0,0,0,0,0.1411,0,0,0,0,0_0,0,0.00304,0,0,0,0,0,0.2822,0,0,0,0,0_0.73_0.00281344827578829
awk '{ split($NF,arr,"_"); print $0"\t"arr[4]"\t"$9"_"$10"_"$11"_"$12"_"$14"_"$15"_"$16"_"$17"_"$18"_"$19 }' targetGeneCorrelations.greater_than_0.7_and_p_less_than_0point01.bed | sort -k1,1 -k2,2n -k15,15 |  bedtools groupby -g 1,2,3,15 -c 21,22 -o collapse,collapse > targetGeneCorrelations.greater_than_0.7_and_p_less_than_0point01.collapsed.bed
 # head targetGeneCorrelations.greater_than_0.7_and_p_less_than_0point01.collapsed.bed
#chr1	87640	87790	ENSG00000268020.2	0.78	chr1_52473_52473_ENST00000606857.1_+_ENSG00000268020.2_KNOWN_OR4G4P_OR4G4P-001_2
#chr1	125340	125490	ENSG00000185097.2	0.73	chr1_622053_622053_ENST00000332831.2_-_ENSG00000185097.2_KNOWN_OR4F16_OR4F16-001_2
#chr1	125340	125490	ENSG00000228463.4	0.73	chr1_259121_259121_ENST00000335577.4_-_ENSG00000228463.4_NOVEL_AP006222.2_AP006222.2-004_2
#chr1	125340	125490	ENSG00000238009.2	0.73,0.73	chr1_129173_129173_ENST00000471248.1_-_ENSG00000238009.2_NOVEL_RP11-34P13.7_RP11-34P13.7-002_2,chr1_129217_129217_ENST00000477740.1_-_ENSG00000238009.2_NOVEL_RP11-34P13.7_RP11-34P13.7-003_2

awk '{
	transcripts_all=""
	split($5, corr,",")
	split($6,transcripts,",")
	max=-1
	argmax=0
	for (f=1; f <= length(corr); f ++)
	{
		if (corr[f] > max)
		{
			max=corr[f]
			argmax=f
		}
		split(transcripts[f],thisTranscript,"_")
		transcripts_all=sprintf("%s %s",transcripts_all,thisTranscript[4])
	}
	split(transcripts[argmax],t,"_")
	printf("%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,t[8],max,t[1],t[2],t[3],t[4],t[5],t[6],t[7],t[9],transcripts_all)
}' targetGeneCorrelations.greater_than_0.7_and_p_less_than_0point01.collapsed.bed > targetGeneCorrelations.greater_than_0.7_and_p_less_than_0point01.maxCorrTranscript.bed
awk -f targetGenes.maxTranscript.awk targetGeneCorrelations.greater_than_0.7_and_p_less_than_0point01.maxCorrTranscript.bed > targetGeneCorrelations.greater_than_0.7_and_p_less_than_0point01.maxCorrTranscript.display.bed
sort -k1,1 -k2,2n targetGeneCorrelations.greater_than_0.7_and_p_less_than_0point01.maxCorrTranscript.display.bed > targetGeneCorrelations.greater_than_0.7_and_p_less_than_0point01.maxCorrTranscript.display.sorted.bed
~/bedToBigBed -tab -type=bed12+7 -as=sandbox/targetGenes.as targetGeneCorrelations.greater_than_0.7_and_p_less_than_0point01.maxCorrTranscript.display.sorted.bed ~/hg19.genome sandbox/targetGene.maxCorrelatedTranscript.bb
# END Requirement to collapse multiple lines to one with highest correlation
