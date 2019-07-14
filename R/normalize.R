
normalize_from_mri <- function(signal, arrival_frame, echotime) {
  signal <- seperate_tissue_signal(signal, arrival_frame)
  return(-1/echotime * log(signal$signal/signal$baseline))
}

normalize_from_ct <- function(signal, arrival_frame=11) {
  signal <- seperate_tissue_signal(signal, arrival_frame)
  return(signal$signal-signal$baseline)
}

seperate_tissue_signal <- function(signal, arrival_frame) {
  baseline <- mean(signal[1:(arrival_frame-1)])
  signal <- signal[arrival_frame:length(signal)]
  return(list(baseline=baseline, signal=signal))
}

