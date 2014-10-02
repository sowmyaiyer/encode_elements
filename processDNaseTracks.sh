Rscript processDNaseTrackForDisplay_distal.R
$HOME/bedToBigBed -tab -type=bed5+3 -as=display/dnase.signaltrack.as allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.signal.normalized.withPeakNames.and_details.bed $HOME/hg19.genome display/allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.signal.normalized.withPeakNames.and_details.bb

Rscript processDNaseTrackForDisplay_proximal.R
$HOME/bedToBigBed -tab -type=bed5+3 -as=display/dnase.signaltrack.as allEncodeDnasePeaks.H3K27ac_only.sorted.merged.proximal.signal.normalized.withPeakNames.and_details.bed $HOME/hg19.genome display/allEncodeDnasePeaks.H3K27ac_only.sorted.merged.proximal.signal.normalized.withPeakNames.and_details.bb
