
vetperf.num <- function(aifnum, roinum, arrival_frame, delta_t, density) {
  # deconvolution
  f_s_res = get_flow_scaled_residue(aif = aifnum, roi = roinum, arrival_frame = arrival_frame, delta_t = delta_t)
  # get parameters
  return(get_perfusion_params(f_s_res, delta_t = delta_t, density = density))
}

vetperf.data <- function(data, arrival_frame, delta_t, density) {
  
}
