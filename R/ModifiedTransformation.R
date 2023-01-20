.transform_shifted_optimisation <- function(parameters, x, type, ...){

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

  return(llf)
}



.transformation_robust_shifted_optimisation <- function(parameters, x, type, weight_method="tukey_window", ...){
  # This follows the algorithm from Raymaekers J, Rousseeuw PJ. Transforming
  # variables to central normality. Mach Learn. 2021.
  # doi:10.1007/s10994-021-05960-5. However, because we are optimising over
  # shift and lambda directly, we can't select an initial
  # lambda-0 value using reweighted optimisation and a first re-weighted
  # optimisation step. We directly use lambda to compute weighted
  # log-likelihood using Tukey's bisquare function as a weight.

  shift <- parameters[1]
  lambda <- parameters[2]

  # Set log-likelihood function.
  loglik_FUN <- switch(
    type,
    "box_cox"=..box_cox_loglik,
    "yeo_johnson"=..yeo_johnson_loglik
  )

  # Make sure that x is sorted.
  if(is.unsorted(x)) x <- sort(x)

  # Apply shift.
  x <- x - shift

  # Set weight function used to compute weights.
  weight_fun <- switch(
    weight_method,
    "tukey_window" = tukey_tapered_cosine_window,
    "trim_transformation" = transformed_step_weighting,
    "trim_residual" = residual_step_weighting,
    "tukey_biweight_residual" = residual_tukey_biweight,
    "huber_weight_residual" = residual_huber_weight)

  # Compute weights
  w <- do.call(
    weight_fun,
    args=c(
      list(
        "x" = x,
        "lambda" = lambda,
        "type" = type),
      list(...)))

  if(!all(is.finite(w))) return(NA_real_)
  if(sum(w) == 0.0) return(NA_real_)

  # Compute log-likelihood.
  llf <- suppressWarnings(do.call(
    loglik_FUN,
    args=list(
      "lambda"=lambda,
      "x"=x,
      "w"=w)))

  return(llf)
}
