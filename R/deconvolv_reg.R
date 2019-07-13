
get_reg_elements <- function(aif_norm, roi_norm, delta_t) {
  t_sym = toeplitz(c(rep(0, length(aif_norm)-1), aif_norm))
  A = t_sym[length(aif_norm):(2 * length(aif_norm)-1), 1:length(aif_norm)]
  # 
  return(list(A = delta_t * A, signal = roi_norm))
}

perfresidue <- function(aif, roi, arrival_frame, delta_t, echo_time, alpha=0.1, lambda = NULL, iterations = 10) {
  aif_norm = normalize_from_mri(aif, arrival_frame = arrival_frame, echotime = echo_time)
  roi_norm = normalize_from_mri(roi, arrival_frame = arrival_frame, echotime = echo_time)
  # get convolution matrix
  elem = get_reg_elements(aif_norm, roi_norm, delta_t)
  # perform l1-l2 weighted fit
  if (is.null(lambda)) {
    lambda = get_min_lambda(elem$A, elem$signal, iterations = iterations)
  }
  # return fit
  cvfit = glmnet(elem$A, elem$signal, alpha=0.1, intercept=FALSE, lower=0., lambda = lambda)
  return(list(residue = coef(cvfit), lambda = lambda, A = elem$A, signal = elem$signal))
}

get_min_lambda <- function(A, signal, iterations) {
  # from
  # https://stats.stackexchange.com/questions/97777/variablity-in-cv-glmnet-results
  # do cross validation to get best lambda and average over many runs
  MSEs <- NULL
  for (i in 1:iterations){
    cv <- cv.glmnet(A, signal, alpha=0.1, intercept=FALSE, lower=0., nfolds = 3)  
    MSEs <- cbind(MSEs, cv$cvm)
  }
  rownames(MSEs) <- cv$lambda
  lambda.min <- as.numeric(names(which.min(rowMeans(MSEs))))
  return(lambda.min)
}
