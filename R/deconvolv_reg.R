
get_reg_elements <- function(aif_norm, roi_norm, delta_t) {
  t_sym = toeplitz(c(rep(0, length(aif_norm)-1), aif_norm))
  A = t_sym[length(aif_norm):(2 * length(aif_norm)-1), 1:length(aif_norm)]
  # 
  return(list(A = delta_t * A, signal = roi_norm))
}

get_flow_scaled_residue <- function(aif, roi, arrival_frame, delta_t, alpha=0.1) {
  aif_norm = normalize_from_mri(aif, arrival_frame = arrival_frame)
  roi_norm = normalize_from_mri(roi, arrival_frame = arrival_frame)
  # get convolution matrix
  elem = get_reg_elements(aif_norm, roi_norm, delta_t)
  # perform l1-l2 weighted fit
  cvfit = cv.glmnet(elem$A, elem$signal, alpha=0.1, intercept=FALSE, lower=0.)
  # return best fit
  return(coef(cvfit, s = "lambda.min"))
}

