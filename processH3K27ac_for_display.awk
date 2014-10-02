{
	for (f=9; f <= 22; f ++)
	{
		split($f,arr,"=")
		if (arr[2] >= 0.95)
		{
			printf("%s\t%d\t%d\t%s\t%.0f\t%.2f\n",$1,$2,$3,arr[1],((arr[2]-0.95)/(1.00-0.95))*1000,arr[2]*100)
		}
	}
}
