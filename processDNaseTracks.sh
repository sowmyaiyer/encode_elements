Rscript processDNaseTrackForDisplay_distal.R
$HOME/bedToBigBed -tab -type=bed9+3 -as=display/dnase.track.as allEncodeDnasePeaks.H3K27ac_only.sorted.merged.distal.signal.normalized.withPeakNames.and_details.bed $HOME/hg19.genome display/dnase_track_distal.bb

Rscript processDNaseTrackForDisplay_proximal.R
$HOME/bedToBigBed -tab -type=bed9+3 -as=display/dnase.track.as allEncodeDnasePeaks.H3K27ac_only.sorted.merged.proximal.signal.normalized.withPeakNames.and_details.bed $HOME/hg19.genome display/dnase_track_proximal.bb
