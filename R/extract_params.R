#' Extract perfusion parameters
#' 
#' Low level function used to extract perfusion parameters from a flow-scaled residue function
#'
#' @param flow_scaled_residue Flow scaled residue function
#' @param delta_t Time between Frames
#' @param density Density of brain tissue
#'
#' @return List with perfusion parameters
#' @export
#' 
get_perfusion_params <- function(flow_scaled_residue, delta_t, density) {
  # calculate mtt
  mtt <- delta_t * sum(flow_scaled_residue) / max(flow_scaled_residue)
  # calculate cbf
  cbf <- 1/density * max(flow_scaled_residue)
  # calculate cbv
  cbv <- 1/density * delta_t * sum(flow_scaled_residue)
  #
  return(list(MTT = mtt, CBF = cbf, CBV = cbv))
}
