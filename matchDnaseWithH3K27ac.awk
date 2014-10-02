{
	printf("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t",$1,$2,$3,$4,$5,$6,$7,$8)
	for (f=9; f <= 22; f ++)
	{
		split($f,arr,"=")
		if (arr[2] >= 0.95)
		{
			printf("%s,",arr[1])	
		}
	}
	printf("\n")
}
