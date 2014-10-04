correlation <- function(x)
{
	if ((sum(x[1:14]) == 0) & (sum(x[15:28]) == 0))
		return(NA)
	else
	{
		distal_profile <- paste(x[1:14], collapse=",")
		gene_profile <- paste(x[15:28], collapse=",")
		cor_results <- cor.test(as.numeric(x[1:14]), as.numeric(x[15:28]), method="spearman")
		cor <- round(cor_results$estimate, digits=2)
		cor_p <- cor_results$p.value
		return (paste(distal_profile, gene_profile, cor, cor_p, sep="_"))
	}
}
input_matrix_file <- commandArgs(TRUE)[1]
output_correlations_file <- commandArgs(TRUE)[2]
cellLines <- scan("cellLinesWithH3K27acAndDnase.txt", what="character")
test <- as.matrix(read.table(input_matrix_file, sep="\t"))
correlations <- apply(test, 1, FUN=correlation)
correlations_with_cellLineInfo <- paste(paste(cellLines, collapse=","),correlations,sep="_")
write(x=correlations_with_cellLineInfo, file=output_correlations_file, sep="\n")
