tukey_tapered_cosine_window <- function(x, k=0.80, ...){
  # Use Tukey window (also cosine-tapered window) for downweighting based on
  # expected quantile (not on residual!). The Tukey window has a flat top, with
  # width k (or 2 alpha), and then tapers off. The alpha parameter is related to
  # k: alpha = 1 - k.

  if(is.unsorted(x)) stop("DEV: x should be sorted.")

  # Iterate along x
  i <- seq_along(x)

  # Initialise weights.
  w <- numeric(length(i)) + 1.0

  # Left tapered cosine lobe.
  left_lobe <- which(i - 0.5 < length(i) * (1 - k) / 2)
  w[left_lobe] <- 0.5 * (1 + cos((2 * pi  * (left_lobe - (1 + length(i) * (1 - k)) /2) / (length(i) * (1 - k)))))

  # Right tapered cosine lobe.
  right_lobe <- which(i - 0.5 > length(i) * (1 + k) /2)
  w[right_lobe] <- 0.5 * (1 + cos((2 * pi * (right_lobe - (1 + length(i) * (1 + k)) /2) / (length(i) * (1 - k)))))

  return(w)
}



transformed_step_weighting <- function(x, lambda, type, tau=2.58, ...){
  # Step weighting used by Rademaekers and Rousseeuw (2021). This basically sets
  # weights for the outerlying transformed values to 0. Default value is to
  # qnorm(0.995).

  transform_FUN <- switch(
    type,
    "box_cox"=..box_cox_transform,
    "yeo_johnson"=..yeo_johnson_transform
  )

  # Find transformed feature values.
  y <- do.call(
    transform_FUN,
    args=list(
      "lambda"=lambda,
      "x"=x))

  # Approximate Huber's M-estimates for locality and scale.
  robust_estimates <- huber_estimate(y, tol=1E-3)

  # Check problematic values.
  if(!is.finite(robust_estimates$sigma)) return(NA_real_)
  if(robust_estimates$sigma == 0.0) return(NA_real_)

  # Compute weights.
  w <- as.numeric(abs(y - robust_estimates$mu) / robust_estimates$sigma <= tau)

  return(w)
}



residual_step_weighting <- function(x, lambda, type, tau=0.50, ...){
  # Step weighting based on residuals.

  # Compute residuals.
  r <- ..compute_residuals(
    x=x,
    lambda=lambda,
    type=type)

  # Set weights.
  w <- as.numeric(abs(r) < tau)

  return(w)
}



residual_tukey_biweight <- function(x, lambda, type, tau=0.50, ...){
  # Set weights based on Tukey's Biweigths.

  # Compute residuals.
  r <- ..compute_residuals(
    x=x,
    lambda=lambda,
    type=type)

  # Initialise 0-weights.
  w <- as.numeric(length(x))

  # Set non-zero weights for residual errors smaller than rho.
  ii <- which(abs(r) <= tau)
  if(length(ii) > 0) w[ii] <- (1 - (r[ii] / tau)^2 )^2

  return(w)
}



residual_huber_weight <- function(x, lambda, type, tau=0.20, ...){
  # Set weights based on Hubers weight.

  # Compute residuals.
  r <- ..compute_residuals(
    x=x,
    lambda=lambda,
    type=type)

  # Initialise 1-weights.
  w <- as.numeric(length(x))

  # Set weights. Values > 0 rho are down-weighted.
  ii <- which(abs(r) <= tau)
  if(length(ii) > 0)  w[ii] <- tau / r[ii]

  return(w)
}



set_full_central_weight <- function(x, w, k=0.80){
  # Assign a full weight to central elements.

  if(is.unsorted(x)) stop("DEV: x should be sorted.")

  # Indices to x.
  ii <- seq_along(x)

  # Find central elements.
  central_elements <- which(length(ii) * (1 - k) / 2 <= ii - 0.5 & ii - 0.5 <= length(ii) * (1 + k) /2)

  # Set central weights to 1
  w[central_elements] <- 1.0

  return(w)
}



..compute_residuals <- function(x, lambda, type){

  # Set transformation function
  transform_FUN <- switch(
    type,
    "box_cox"=..box_cox_transform,
    "yeo_johnson"=..yeo_johnson_transform
  )

  # Find transformed feature values.
  y <- do.call(
    transform_FUN,
    args=list(
      "lambda"=lambda,
      "x"=x))

  # Compute the expected z-score.
  z_expected <- compute_expected_z(x=x)

  # Approximate Huber's M-estimates for locality and scale.
  robust_estimates <- huber_estimate(y, tol=1E-3)

  # Check problematic values.
  if(!is.finite(robust_estimates$sigma)) return(NA_real_)
  if(robust_estimates$sigma == 0.0) return(NA_real_)

  # Compute the observed z-score.
  z_observed <- (y - robust_estimates$mu) / robust_estimates$sigma

  # Compute residuals.
  residual <- z_observed - z_expected

  return(residual)
}
