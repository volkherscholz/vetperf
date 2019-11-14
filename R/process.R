

#' Calculate perfusion parameters from numeric vectors
#' 
#' \code{vetperf.num} takes an arterial input function as well as a ROI measurement and
#' computes the perfusion parameters MTT, CBF and CBV
#' 
#' Internally, it performs a l1-l2 weighted regression employing cross-validation to obtain the best
#' fit for the response
#'
#' @param aifnum Numeric vector, the arterial input function
#' @param roinum Numeric vector, the ROI measurement
#' @param arrival_frame Integer, the frame number at which the bolus arrived
#' @param delta_t Float, the time between two consecutive images
#' @param echo_time Float, MRI echo time of the sequence
#' @param density Float, density of the brain (optional)
#' @param alpha Float, elastic-net hyperparameter (l1 vs l2 weight), optional, default 0.1
#' @param lambda Float, l1-l2 glmnet hyperparameter. Will be determined via cross-validation
#'   if not provided (recommend)
#' @param iterations Integer, the number of iterations used for cross-validation (optional, default 10)
#'
#' @return List with MTT, CBF, CBV parameters
#' @export 
#' 
vetperf.num <- function(aifnum, roinum, arrival_frame, delta_t, echo_time, density = 1.04, alpha = 0.1, lambda = NULL, iterations = 10) {
  # deconvolution
  f_s_res_fit <- perfresidue(
    aif = aifnum,
    roi = roinum,
    arrival_frame = arrival_frame,
    delta_t = delta_t,
    echo_time = echo_time,
    alpha = alpha,
    lambda = lambda,
    iterations = iterations)
  f_s_res <- f_s_res_fit$residue
  # get parameters
  perfparams <- get_perfusion_params(f_s_res, delta_t = delta_t, density = density)
  return(perfparams)
}

.datatable.aware=TRUE;
#' Calculate perfusion parameters
#'
#' \code{vetperf.data} takes a data.table object containing AIF and ROIs measurements
#' and replaces them with perfusion parameters
#'
#' High level function to obtain perfusion parameters from aif and roi measurements.
#' Given a data.table object holding measurements for one or more animals, one or more sides
#' and one or more ROIs, this functions extracts the associated perfusion parameters 
#' and optionally normalizes them with respect to a ROI of choice.
#' The data.table is assumed to have five columns:
#'   animal: factor indicating the animal
#'   side: factor indicating the side
#'   roi: factor indicating the ROI, must include a level corresponding to
#'        the arterial input function
#'   timestep: column indicating the frame at which the measurement was taken
#'   value: column holding the measurements
#'
#' @param data Data.table object containing aif and roi measurements
#' @param arrival_frame Integer, the frame number at which the bolus arrived
#' @param delta_t Float, the time between two consecutive images
#' @param echo_time Float, MRI echo time of the sequence
#' @param animalcol String, the name of the column holding the animal information
#' @param sidecol String, the name of the column indicating which side is measured
#' @param roicol String, the name of the column which indicates the ROI
#' @param aifroi String, the ROI which contains the AIF measurement (default "aif")
#' @param valcol String, the name of the column which holds the data
#' @param baseroi String, the name of the white matter roi which should be used for normalization 
#'   (if normalized values should be returned)
#' @param density Float, density of the brain (optional)
#' @param alpha Float, elastic-net hyperparameter (l1 vs l2 weight), optional, default 0.1
#' @param lambda Float, l1-l2 glmnet hyperparameter. Will be determined via cross-validation
#'   if not provided (recommend)
#' @param iterations Integer, the number of iterations used for cross-validation (optional, default 10)
#'
#' @return Data.table object holding perfusion parameters (MTT, CBF, CBV)
#' @export
#'
#' @examples
#' perfdata = vetperf.data(animals, arrival_frame = 9, delta_t = 1.6, echo_time = 0.03, baseroi = "cs", iterations = 1)
#' 
vetperf.data <- function(
  data,
  arrival_frame,
  delta_t,
  echo_time,
  animalcol = "animal",
  sidecol = "side",
  roicol = "roi",
  aifroi = "aif",
  valcol = "value",
  baseroi = NULL,
  density = 1.04,
  alpha = 0.1,
  lambda = NULL,
  iterations = 10) {
  # rename columns to defaults
  data.table::setnames(data, c(animalcol, sidecol, roicol, valcol), c("animal", "side", "roi", "value"))
  #.animal <- eval(substitute(animalcol))
  #.side <- eval(substitute(sidecol))
  #.roi <- eval(substitute(roicol))
  #.value <- eval(substitute(valcol))
  
  # define internal function hardcoding the hyperparameters
  perfparams <- function(data) {
    aif <- data[roi == aifroi, value]
    params = data[
      roi != aifroi,
      vetperf.num(
        aif, value,
        arrival_frame = arrival_frame,
        delta_t = delta_t,
        echo_time = echo_time,
        density = density,
        alpha = alpha,
        lambda = lambda,
        iterations = iterations),
      by = roi]
    
    # check whether to normalize rois
    if (!is.null(baseroi)) {
      paramsnorm <- params[roi == baseroi, c("MTT", "CBF", "CBV")]
      params <- params[, .SD / paramsnorm, by = roi, .SDcols = c("MTT", "CBF", "CBV")]
      params <-  params[roi != baseroi]
    }
    return(params)
  }
  
  perfdata <- data[, perfparams(.SD), by = list(animal, side)]
  # rename to original
  data.table::setnames(data, c("animal", "side", "roi", "value"), c(animalcol, sidecol, roicol, valcol))
  return(perfdata)
}
