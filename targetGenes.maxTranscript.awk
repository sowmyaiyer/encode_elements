#chr1	770925	771075	RP11-206L10.11	0.77	chr1	766985	766985	ENST00000415295.1	+	ENSG00000228794.4	NOVEL	RP11-206L10.11-007	 ENST00000415295.1
#chr1	773060	773210	RP11-206L10.8	0.75	chr1	745541	745541	ENST00000447500.1	-	ENSG00000230092.3	KNOWN	RP11-206L10.8-002	 ENST00000447500.1 ENST00000590817.1
{
	if ($7 < $2)
	{
		chr=$6
		start=$7
		end=$3
		name=$4
		strand="-"
		thickStart=start
		thickEnd=end
		itemRgb=0
		blockCount=2
		blockSize1=1
		blockSize2=$3-$2
		blockStart1=0
		blockStart2=$2-start
	} else {
		chr=$1
		start=$2
		end=$7
		name=$4
		strand="+"
		thickStart=start
		thickEnd=end
		itemRgb=0
		blockCount=2
		blockSize1=$3-$2
		blockSize2=1
		blockStart1=0
		blockStart2=$7-start-1
	}
printf("%s\t%d\t%d\t%s\t%.0f\t%s\t%d\t%d\t%d\t%d\t%d,%d\t%d,%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",chr,start,end,name,$5*1000,strand,thickStart,thickEnd,itemRgb,blockCount,blockSize1,blockSize2,blockStart1,blockStart2,$11,$4,$13,$9,$12,$10,$14,$5)
}
