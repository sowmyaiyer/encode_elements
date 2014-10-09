normalizeto1000Scale <- function(x, max, min)
{
	x <- as.numeric(x)
	return (round((((x-min)/(max-min)) * 1000)))
}
getCleanCellLineNames <- function(x)
{
	cellLines <- unlist(strsplit(x, split=", "))
	cellLineCleanNames <- as.character(sapply(cellLines, FUN=getCleanCellName))
	return (paste(cellLineCleanNames, collapse=","))
}
getCleanCellName <- function(cellLine)
{
	return (cellLineMapping[cellLine,1])
}

#chr1	540640	540790	Hsmm	38	Hsmm,Nhdfad	38,29	2
#chr1	540880	541030	Nhdfad	26	Nhdfad	26	1
#chr1	559280	559430	Nhdfad	59.000000	Hsmm,Nhdfad	43,59.000000	2
cellLineMapping <- as.data.frame(read.table("cellLineMapping", row.names=1, sep="\t"))
df <- as.data.frame(read.table("allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.bed"))
#df <- as.data.frame(read.table("test.bed",sep="\t"))
cat("read file\n")
max_sig <- max(as.numeric(df[,5]))
min_sig <- min(as.numeric(df[,5]))
normalized_score_signal <- sapply(df[,5], FUN=normalizeto1000Scale, max=max_sig, min=min_sig)

cell_types <- strsplit(as.character(df[,6]), split=",")
unique_cell_types <- lapply(cell_types, FUN=unique)
unique_cell_types_formatted <- gsub(gsub(gsub(as.character(unique_cell_types), pattern="\"",replacement="",perl=TRUE), pattern="\\)", replacement="",perl=TRUE), pattern="c\\(", replacement="",perl=TRUE)

unique_cell_types_clean <- sapply(unique_cell_types_formatted, FUN=getCleanCellLineNames)
total_peaks_merged <- df[,8]
number_unique_cell_types_merged <- as.numeric(lapply(unique_cell_types, FUN=length))

dnase_signal_track <- data.frame(df[,1],df[,2],df[,3],paste("re",gsub(df[,1],pattern="chr",replacement=""),".",floor((df[,2]+df[,3])/2),sep=""),normalized_score_signal, ".", df[,2],df[,3], "94,207,74", unique_cell_types_clean,number_unique_cell_types_merged,total_peaks_merged,gsub(df[,7],pattern="\\.000000", replacement=""))

write.table(x=dnase_signal_track, file="allEncodeDnasePeaks.normalized.withPeakNames.and_details.bed", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
