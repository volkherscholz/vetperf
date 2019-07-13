
vetperf.num <- function(aifnum, roinum, arrival_frame, delta_t, echo_time, density = 1.04, lambda = NULL, iterations = 10) {
  # deconvolution
  f_s_res_fit = perfresidue(
    aif = aifnum,
    roi = roinum,
    arrival_frame = arrival_frame,
    delta_t = delta_t,
    echo_time = echo_time,
    lambda = lambda,
    iterations = iterations)
  f_s_res = f_s_res_fit$residue
  # get parameters
  return(get_perfusion_params(f_s_res, delta_t = delta_t, density = density))
}

vetperf.data <- function(
  data,
  arrival_frame,
  delta_t,
  echo_time,
  animalcol = "animal",
  sidecol = "side",
  roicol = "roi",
  valcol = "value",
  baseroi = NULL,
  density = 1.04,
  lambda = NULL,
  iterations = 10) {
  # define internal function hardcoding the hyperparameters
  perfparams <- function(data) {
    aif = data[roi == 'aif', get(valcol)]
    params = data[
      get(roicol) != 'aif',
      vetperf.num(
        aif, get(valcol),
        arrival_frame = arrival_frame,
        delta_t = delta_t,
        echo_time = echo_time,
        density = density,
        lambda = lambda,
        iterations = iterations),
      by = get(roicol)]
    
    # check whether to normalize rois
    if (!is.null(baseroi)) {
      paramsnorm = params[get == baseroi, c("MTT", "CBF", "CBV")]
      params = params[, .SD / paramsnorm, by = get, .SDcols = c("MTT", "CBF", "CBV")]
      params = params[get != baseroi]
    }
    # rename roi column
    setnames(params, "get", roicol)
    return(params)
  }
  
  perfdata = data[, perfparams(.SD), by = list(get(animalcol), get(sidecol))]
  setnames(perfdata, c("get", "get.1"), c(animalcol, sidecol))
  return(perfdata)
}
