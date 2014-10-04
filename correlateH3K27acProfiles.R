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

cellLines <- scan("cellLinesWithH3K27acAndDnase.txt", what="character")
count <- 1
for (cellLine in cellLines)
{
	cat(cellLine,"\n")
	h3k27ac_master_peak <- as.numeric(scan(paste("/home/si14w/nearline/decorate_dnase/targetGenes/H3K27ac_signal.masterPeaks.",cellLine,".txt",sep="")))
	if (count == 1)
        {
                master_matrix <- matrix(nrow=length(h3k27ac_master_peak))
        }
	master_matrix <- cbind(master_matrix, h3k27ac_master_peak)
	count <- count + 1
}
for (cellLine in cellLines)
{
	cat(cellLine,"\n")
	h3k27ac_nearby_gene <- as.numeric(scan(paste("/home/si14w/nearline/decorate_dnase/targetGenes/H3K27ac_signal.targetgenes.",cellLine,".txt",sep="")))
	master_matrix <- cbind(master_matrix, h3k27ac_nearby_gene)
}

cat(nrow(master_matrix), ncol(master_matrix), "\n")
test <- master_matrix[,-1]
write.table(test, file="H3K27ac_signal_master_peaks_and_target_genes.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
