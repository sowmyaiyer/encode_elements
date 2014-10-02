{
	found=0
	for (f=9; f <= 22; f ++)
	{
		split($f,arr,"=")
		if (arr[2] >= 0.95)
		{
			found=1
		}
	}
	if (found == 0)
	{
		print $0 >> "noSignificantH3K27ac.bed"
	}
}
