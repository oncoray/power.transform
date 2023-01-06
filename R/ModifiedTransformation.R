.transform_shifted_optimisation <- function(parameters, shift_range, lambda_range, x, type, ...){

  shift <- parameters[1]
  lambda <- parameters[2]

  # Check boundary conditions.
  if(shift < shift_range[1] || shift > shift_range[2] || lambda < lambda_range[1] || lambda > lambda_range[2]) return(NA_real_)

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



.transformation_robust_shifted_optimisation <- function(parameters, shift_range, lambda_range, x, z, type, weight_method="tailed", k=0.5, central_weight=0.80, ...){
  # This follows the algorithm from Raymaekers J, Rousseeuw PJ. Transforming
  # variables to central normality. Mach Learn. 2021.
  # doi:10.1007/s10994-021-05960-5. However, because we are optimising over
  # shift and lambda directly, we can't select an initial
  # lambda-0 value using reweighted optimisation and a first re-weighted
  # optimisation step. We directly use lambda to compute weighted
  # log-likelihood using Tukey's bisquare function as a weight.

  if(!weight_method %in% c("tailed", "huber", "tukey", "step")){
    stop(paste0("DEV: The weight_method argument was not recognised: ", weight_method))
  }

  if(weight_method == "tailed" && is.null(central_weight)){
    stop(paste0("The central_weight argument cannot be NULL if the weigh_method is \"tailed\"."))
  }

  shift <- parameters[1]
  lambda <- parameters[2]

  # Check boundary conditions.
  if(shift < shift_range[1] || shift > shift_range[2] || lambda < lambda_range[1] || lambda > lambda_range[2]) return(NA_real_)

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

  # Make sure that x is sorted.
  if(is.unsorted(x)) x <- sort(x)

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
  residual <- (y - robust_estimates$mu) / robust_estimates$sigma - z

  if(weight_method == "tailed"){
    # Use Tukey window (also cosine-tapered window) for downweighting based on
    # expected quantile (not on residual!). The Tukey window has a flat top, and
    # then tapers off. The alpha parameter is related to central_weight: alpha =
    # 1 - central_weight.
    weights <- numeric(length(z)) + 1.0

    # Generate tails.
    tail_width <- ceiling(0.5 * length(z) * (1.0 - central_weight))
    tail_weights <- 0.5 * (1.0 - cos(2.0 * pi * seq_len(tail_width) / (length(z) * (1.0 - central_weight))))
    weights[seq_len(tail_width)] <- tail_weights
    weights[length(z) - seq_len(tail_width) + 1L] <- tail_weights

  } else if(weight_method == "huber"){
    # Huber weights using k=0.5 after Rousseuw and Raymaekers. This
    weights <- numeric(length(residual)) + 1.0
    outlier_residuals <- which(abs(residual) > k)

    if(length(outlier_residuals) > 0){
      weights[outlier_residuals] <- k / abs(residual[outlier_residuals])
    }

  } else if(weight_method == "tukey"){
    # Compute Tukey bisquare function to truncate weights of outlier residuals.
    # We use k=0.5 after Rousseeuw and Raymaekers.
    weights <- numeric(length(residual)) + 1.0
    valid_residuals <- which(abs(residual) <= k)

    if(length(valid_residuals) > 0){
      weights[valid_residuals] <- 1.0 - (1.0 - (residual[valid_residuals] / k)^2)^3
    }

    # We want to maximise log-likelihood, so we need to invert the weights.
    weights <- 1.0 - weights

  } else if(weight_method == "step"){
    weights <- as.numeric(abs(residual) < stats::qnorm(0.90))
  }

  # Assign full weight to central elements. We do this to ensure that poor fits
  # don't get down-weighted in the centre, and in that way produce
  # log-likelihood values that are too optimistic.
  if(!is.null(central_weight) && weight_method %in% c("huber", "tukey", "step")){
    if(central_weight < 0.0 || central_weight > 1.0){
      stop(paste0("DEV: central_weight should be NULL or between 0.0 and 1.0. Found: ", central_weight))
    }

    central_values <- which(abs(z) < stats::qnorm(0.50 + central_weight / 2))
    weights[central_values] <- 1.0
  }

  if(!all(is.finite(weights))) return(NA_real_)
  if(sum(weights) == 0.0) return(NA_real_)

  # Compute log-likelihood.
  llf <- suppressWarnings(do.call(
    loglik_FUN,
    args=list(
      "lambda"=lambda,
      "x"=x,
      "w"=weights)))

  return(llf)
}
