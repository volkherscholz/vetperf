
get_perfusion_params <- function(flow_scaled_residue, delta_t, density) {
  # calculate mtt
  mtt = delta_t * sum(flow_scaled_residue) / max(flow_scaled_residue)
  # calculate cbf
  cbf = 1/density * max(flow_scaled_residue)
  # calculate cbv
  cbv = 1/density * delta_t * sum(flow_scaled_residue)
  #
  return(list(mtt = mtt, cbf = cbf, cbv = cbv))
}
