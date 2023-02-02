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



.transformation_robust_shifted_optimisation <- function(parameters, x, type, weight_method="original_cosine", ...){
  # Compute weighted log-likelihood function.

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
    "original_step" = original_step,
    "original_triangle" = original_triangle,
    "original_cosine" = original_cosine,
    "transformed_step" = transformed_step,
    "transformed_triangle" = transformed_triangle,
    "transformed_cosine" = transformed_cosine,
    "residual_step" = residual_step,
    "residual_triangle" = residual_triangle,
    "residual_cosine" = residual_cosine)

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
