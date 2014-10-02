#wget "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
#gunzip gencode.v19.annotation.gtf.gz
awk -F" " '{
	if ($3=="gene" && $16~/\"KNOWN\"/ && $26~/1|2/)
	{
		gsub("\"","",$18)
		gsub(";","",$18)
		gsub("\"","",$10)
		gsub(";","",$10)
		if ($7 == "+")
		{
			start=$4
				
		} else if ($7 == "-") {
			start=$5
		}
		printf("%s\t%d\t%d\t%s_%s\n",$1,start-500, start+500, $18,$10)
	}
}' gencode.v19.annotation.gtf > gencode.v19.genes.tss.plusminus500bp.bed


