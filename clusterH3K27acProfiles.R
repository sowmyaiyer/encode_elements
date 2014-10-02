allH3K27ac <- as.matrix(read.table("allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.H3K27ac_percentiles.numbers.for_clustering.txt",sep="\t", header=TRUE, row.names=NULL))
#allH3K27ac <- as.data.frame(read.table("test.txt",sep="\t", header=TRUE, row.names=NULL))
cat("read file\n")
kmeans_results <- kmeans(allH3K27ac, centers=5)
m <- cbind(allH3K27ac, kmeans_results$cluster)

write.table(m, file="allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.clusters.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
