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
df <- as.data.frame(read.table("allEncodeDnasePeaks.H3K27ac_only.sorted.merged.proximal.bed"))
#df <- as.data.frame(read.table("test.bed",sep="\t"))
cat("read file\n")
print(df)
max_sig <- max(as.numeric(df[,5]))
min_sig <- min(as.numeric(df[,5]))

#max_number <- max(as.numeric(df[,8]))
#min_number <- min(as.numeric(df[,8]))

#max_composite <- max(as.numeric(df[,5])*as.numeric(df[,8]))
#min_composite <- min(as.numeric(df[,5])*as.numeric(df[,8]))

normalized_score_signal <- sapply(df[,5], FUN=normalizeto1000Scale, max=max_sig, min=min_sig)
#normalized_score_number <- sapply(df[,8], FUN=normalizeto1000Scale, max=max_number, min=min_number)
#normalized_score_composite <- sapply(df[,5]*df[,8], FUN=normalizeto1000Scale, max=max_composite, min=min_composite)

##dnase_signal_track <- data.frame(df[,1],df[,2],df[,3],paste("masterPeak_",1:nrow(df),sep=""),normalized_score_signal)
#dnase_signal_track <- data.frame(df[,1],df[,2],df[,3],df[,8],normalized_score_signal)
cell_types <- strsplit(as.character(df[,6]), split=",")
unique_cell_types <- lapply(cell_types, FUN=unique)
unique_cell_types_formatted <- gsub(gsub(gsub(as.character(unique_cell_types), pattern="\"",replacement="",perl=TRUE), pattern="\\)", replacement="",perl=TRUE), pattern="c\\(", replacement="",perl=TRUE)
unique_cell_types_clean <- sapply(unique_cell_types_formatted, FUN=getCleanCellLineNames)
total_peaks_merged <- df[,8]
number_unique_cell_types_merged <- as.numeric(lapply(unique_cell_types, FUN=length))
print(number_unique_cell_types_merged)
# Added following like 09/25/2014
dnase_signal_track <- data.frame(df[,1],df[,2],df[,3],paste("re",gsub(df[,1],pattern="chr",replacement=""),".",floor((df[,2]+df[,3])/2),sep=""),normalized_score_signal, unique_cell_types_clean,number_unique_cell_types_merged,total_peaks_merged)
#dnase_number_track <- data.frame(df[,1],df[,2],df[,3],paste("masterPeak_",1:nrow(df),sep=""),normalized_score_number)
#dnase_composite_track <- data.frame(df[,1],df[,2],df[,3],paste("masterPeak_",1:nrow(df),sep=""),normalized_score_composite)

#write.table(x=dnase_signal_track, file="allEncodeDnasePeaks.H3K27ac_only.sorted.merged.proximal.signal.normalized.withPeakNames.bed", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(x=dnase_signal_track, file="allEncodeDnasePeaks.H3K27ac_only.sorted.merged.proximal.signal.normalized.withPeakNames.and_details.bed", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
#write.table(x=dnase_number_track, file="allEncodeDnasePeaks.H3K27ac_only.sorted.merged.proximal.number.bed", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
#write.table(x=dnase_composite_track, file="allEncodeDnasePeaks.H3K27ac_only.sorted.merged.proximal.signal_times_number.bed", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
