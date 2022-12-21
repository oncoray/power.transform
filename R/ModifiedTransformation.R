.transform_shifted_optimisation <- function(parameters, shift_range, lambda_range, x, type){

  shift <- parameters[1]
  lambda <- parameters[2]

  # Set log-likelihood function.
  loglik_FUN <- switch(
    type,
    "box_cox"=..box_cox_loglik,
    "yeo_johnson"=..yeo_johnson_loglik
  )

  # Apply shift.
  x <- x - shift

  # Compute log-likelihood.
  llf <- suppressWarnings(do.call(
    loglik_FUN,
    args=list(
      "lambda"=lambda,
      "x"=x)))

  # Check boundary conditions and update llf, if necessary. This function steers
  # the optimiser away from the boundary.
  llf <- apply_boundary(
    llf=llf,
    shift=shift,
    shift_range=shift_range,
    lambda=lambda,
    lambda_range=lambda_range)

  return(llf)
}



.transformation_robust_shifted_optimisation <- function(parameters, shift_range, lambda_range, x, type){
  # This follows the algorithm from Raymaekers J, Rousseeuw PJ. Transforming
  # variables to central normality. Mach Learn. 2021.
  # doi:10.1007/s10994-021-05960-5. However, because we are optimising over
  # shift and lambda directly, we don't really need to select an initial
  # lambda-0 value using rectified optimisation and a first re-weighted
  # optimisation step. We directly use lambda to compute weighted
  # log-likelihood.

  shift <- parameters[1]
  lambda <- parameters[2]

  # Set transformation function.
  transform_FUN <- switch(
    type,
    "box_cox"=..box_cox_transform,
    "yeo_johnson"=..yeo_johnson_transform
  )

  # Set log-likelihood function.
  loglik_FUN <- switch(
    type,
    "box_cox"=..box_cox_loglik,
    "yeo_johnson"=..yeo_johnson_loglik
  )

  # Apply shift.
  x <- x - shift

  # Perform transformation for lambda.
  y <- suppressWarnings(do.call(
    transform_FUN,
    args=list(
      "lambda"=lambda,
      "x"=x)))

  # Compute M-estimates for locality and scale
  robust_estimates <- huber_estimate(y)

  # Check problematic values.
  if(!is.finite(robust_estimates$sigma)) return(NA_real_)
  if(robust_estimates$sigma == 0.0) return(NA_real_)

  # Compute weights.
  weights <- as.numeric(abs(y - robust_estimates$mu) / robust_estimates$sigma <= stats::qnorm(0.99))
  if(!all(is.finite(weights))) return(NA_real_)
  if(sum(weights) == 0.0) return(NA_real_)

  # Compute log-likelihood.
  llf <- suppressWarnings(do.call(
    loglik_FUN,
    args=list(
      "lambda"=lambda,
      "x"=x,
      "w"=weights)))

  # Check boundary conditions and update llf, if necessary. This function steers
  # the optimiser away from the boundary.
  llf <- apply_boundary(
    llf=llf,
    shift=shift,
    shift_range=shift_range,
    lambda=lambda,
    lambda_range=lambda_range)

  return(llf)
}
