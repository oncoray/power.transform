.transformation_robust_optimisation <- function(x, type){
  # This follows the algorithm from Raymaekers J, Rousseeuw PJ. Transforming
  # variables to central normality. Mach Learn. 2021.
  # doi:10.1007/s10994-021-05960-5

  # Sort x.
  x <- sort(x)

  # Compute z-values according to the inverse cumulative density function.
  z_expected <- compute_expected_z(x=x)

  # Step 1: Compute initial estimate for lambda.
  # Standard method based on optimising log-likelihood of the normal
  # distribution.

  optimal_lambda <- suppressWarnings(
    stats::optimise(
      ..transformation_rectified_optimisation,
      interval=c(-4, 4),
      x=x,
      z=z_expected,
      type=type,
      maximum=FALSE))

  optimal_lambda <- ifelse(
    is.finite(optimal_lambda$objective),
    optimal_lambda$minimum,
    1.0)

  # Step 2: Compute lambda from reweighted maximum likelihood.
  optimal_lambda <- ..transformation_reweighted_optimisation(
    lambda_0=optimal_lambda,
    x=x,
    type=type,
    ii=1L)

  # Step 3: Compute lambda from reweighted maximum likelihood again.
  optimal_lambda <- ..transformation_reweighted_optimisation(
    lambda_0=optimal_lambda,
    x=x,
    type=type,
    ii=2L)

  if(!is.finite(optimal_lambda)) return(1.0)

  return(optimal_lambda)
}



..transformation_rectified_optimisation <- function(lambda, x, z, type){

  rectifier_FUN <- switch(
    type,
    "box_cox"=..box_cox_transform_rectified,
    "yeo_johnson"=..yeo_johnson_transform_rectified
  )

  # Compute values after power transformation with the appropriate lambda value.
  y <- do.call(
    rectifier_FUN,
    args=list(
      "lambda"=lambda,
      "x"=x))

  # Compute M-estimates for locality and scale
  robust_estimates <- huber_estimate(y)

  # Check problematic values.
  if(!is.finite(robust_estimates$sigma)) return(NA_real_)
  if(robust_estimates$sigma == 0.0) return(NA_real_)

  # Compute residuals.
  residual <- (y - robust_estimates$mu) / robust_estimates$sigma - z

  # Compute Tukey bisquare function to truncate weights of outlier residuals.
  truncated_weights <- numeric(length(residual)) + 1.0
  valid_residuals <- which(abs(residual) <= 0.5)

  if(length(valid_residuals) > 0){
    truncated_weights[valid_residuals] <- 1.0 - (1.0 - (residual[valid_residuals] / 0.5)^2)^3
  }

  return(sum(truncated_weights))
}



..transformation_reweighted_optimisation <- function(lambda_0, x, type, ii){

  if(!is.finite(lambda_0)) return(NA_real_)

  # Set transformation function.
  if(ii == 1){
    # For the initial step, use the rectified transformations
    transform_FUN <- switch(
      type,
      "box_cox"=..box_cox_transform_rectified,
      "yeo_johnson"=..yeo_johnson_transform_rectified
    )

  } else {
    transform_FUN <- switch(
      type,
      "box_cox"=..box_cox_transform,
      "yeo_johnson"=..yeo_johnson_transform
    )
  }

  # Set log-likelihood function.
  loglik_FUN <- switch(
    type,
    "box_cox"=..box_cox_loglik,
    "yeo_johnson"=..yeo_johnson_loglik
  )

  # Perform transformation for lambda_0
  y <- do.call(
    transform_FUN,
    args=list(
      "lambda"=lambda_0,
      "x"=x
    )
  )

  # Compute M-estimates for locality and scale
  robust_estimates <- huber_estimate(y)

  # Check problematic values.
  if(!is.finite(robust_estimates$sigma)) return(NA_real_)
  if(robust_estimates$sigma == 0.0) return(NA_real_)

  # Compute weights.
  weights <- as.numeric(abs(y - robust_estimates$mu) / robust_estimates$sigma <= stats::qnorm(0.99))
  if(sum(weights) == 0.0) return(NA_real_)

  # Compute optimal lambda.
  optimal_lambda <- suppressWarnings(
    stats::optimise(
      loglik_FUN,
      interval=c(-4.0, 4.0),
      x=x,
      w=weights,
      maximum=TRUE))

  lambda <- ifelse(
    is.finite(optimal_lambda$objective),
    optimal_lambda$maximum,
    1.0)

  return(lambda)
}
