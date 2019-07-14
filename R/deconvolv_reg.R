#' Extract a flow scaled residue function from AIF and ROI measurements
#' 
#' \code{perfresidue} is a low-level function which givem AIF and ROI measurements 
#' obtains an estimate for the flow-scaled residue function
#' 
#' We refer to the Review Article by Fieselmann et. al. (see \url{https://doi.org/10.1155/2011/467563})
#' for detailed definitions. This function implements an improved version of the algorithm described
#' in the previously mentioned paper, employing machine learning techniques to obtain
#' better estimates for the flow-scaled residue function.
#'
#' @param aif Numeric vector, the arterial input function
#' @param roi Numeric vector, the ROI measurement
#' @param arrival_frame Integer, the frame number at which the bolus arrived
#' @param delta_t Float, the time between two consecutive images
#' @param echo_time Float, MRI echo time of the sequence
#' @param alpha Float, l1 hyperparameter (l1 vs l2 weight), optional, default 0.1
#' @param lambda Float, l1-l2 glmnet hyperparameter. Will be determined via cross-validation
#'   if not provided (recommend)
#' @param iterations Integer, the number of iterations used for cross-validation (optional, default 10)
#'
#' @return List with the flow-scaled residue function, the fitted lambda, the regression matrix 
#'   as well as the signal used (rescaled ROI information)
#' @export
perfresidue <- function(aif, roi, arrival_frame, delta_t, echo_time, alpha=0.1, lambda = NULL, iterations = 10) {
  aif_norm <- normalize_from_mri(aif, arrival_frame = arrival_frame, echotime = echo_time)
  roi_norm <- normalize_from_mri(roi, arrival_frame = arrival_frame, echotime = echo_time)
  # get convolution matrix
  elem <- get_reg_elements(aif_norm, roi_norm, delta_t)
  # perform l1-l2 weighted fit
  if (is.null(lambda)) {
    lambda <- get_min_lambda(elem$A, elem$signal, iterations = iterations)
  }
  # return fit
  cvfit <- glmnet::glmnet(elem$A, elem$signal, alpha=0.1, intercept=FALSE, lower=0., lambda = lambda)
  return(list(residue = stats::coef(cvfit), lambda = lambda, A = elem$A, signal = elem$signal))
}


get_reg_elements <- function(aif_norm, roi_norm, delta_t) {
  t_sym <- stats::toeplitz(c(rep(0, length(aif_norm)-1), aif_norm))
  A <- t_sym[length(aif_norm):(2 * length(aif_norm)-1), 1:length(aif_norm)]
  # 
  return(list(A = delta_t * A, signal = roi_norm))
}


get_min_lambda <- function(A, signal, iterations) {
  # from
  # https://stats.stackexchange.com/questions/97777/variablity-in-cv-glmnet-results
  # do cross validation to get best lambda and average over many runs
  MSEs <- NULL
  for (i in 1:iterations){
    cv <- glmnet::cv.glmnet(A, signal, alpha=0.1, intercept=FALSE, lower=0., nfolds = 3)  
    MSEs <- cbind(MSEs, cv$cvm)
  }
  rownames(MSEs) <- cv$lambda
  lambda.min <- as.numeric(names(which.min(rowMeans(MSEs))))
  return(lambda.min)
}
