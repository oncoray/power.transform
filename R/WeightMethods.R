original_step <- function(x, k1=0.80, ...){
  return(.step_window(
    x = .compute_weight_original(x = x),
    k1 = k1))
}

original_triangle <- function(x, k1=0.70, k2=0.90, ...){
  return(.triangular_window(
    x = .compute_weight_original(x = x),
    k1 = k1,
    k2 = k2))
}

original_cosine <- function(x, k1=0.70, k2=0.90, ...){
  return(.tapered_cosine_window(
    x = .compute_weight_original(x = x),
    k1 = k1,
    k2 = k2))
}

transformed_step <- function(x, lambda, type, k1=stats::qnorm(0.90), ...){
  return(.step_window(
    x = .compute_weight_transformed(
      x = x,
      lambda = lambda,
      type = type),
    k1 = k1))
}

transformed_triangle <- function(x, lambda, type, k1=stats::qnorm(0.85), k2=stats::qnorm(0.95), ...){
  return(.triangular_window(
    x = .compute_weight_transformed(
      x = x,
      lambda = lambda,
      type = type),
    k1 = k1,
    k2 = k2))
}

transformed_cosine <- function(x, lambda, type, k1=stats::qnorm(0.85), k2=stats::qnorm(0.95), ...){
  return(.tapered_cosine_window(
    x = .compute_weight_transformed(
      x = x,
      lambda = lambda,
      type = type),
    k1 = k1,
    k2 = k2))
}

residual_step <- function(x, lambda, type, k1=0.40, ...){
  return(.step_window(
    x = .compute_weight_residual(
      x = x,
      lambda = lambda,
      type = type),
    k1 = k1))
}

residual_triangle <- function(x, lambda, type, k1=0.30, k2=0.50, ...){
  return(.triangular_window(
    x = .compute_weight_residual(
      x = x,
      lambda = lambda,
      type = type),
    k1 = k1,
    k2 = k2))
}

residual_cosine <- function(x, lambda, type, k1=0.30, k2=0.50, ...){
  return(.tapered_cosine_window(
    x = .compute_weight_residual(
      x = x,
      lambda = lambda,
      type = type),
    k1 = k1,
    k2 = k2))
}



.compute_weight_original <- function(x){
  # Check if x is sorted.
  if(is.unsorted(x)) stop(paste0("DEV: x is expected to be sorted in ascending order."))

  # Compute expected quantile.
  q <- (seq_along(x) - 1/3) / (length(x) + 1/3)

  # Centralise and map to [-1, 1] range.
  q <- 2.0 * (q - 0.5)

  return(q)
}



.compute_weight_transformed <- function(x, lambda, type){

  # Set transformer.
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

  return((y - robust_estimates$mu) / robust_estimates$sigma)
}



.compute_weight_residual <- function(x, lambda, type){
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



.step_window <- function(x, k1, ...){
  # Set weights using a step window.

  # k1 should be 0 or greater.
  if(k1 < 0) return(NA_real_)

  # Initialise weights.
  w <- numeric(length(x)) + 1.0

  # Set step window.
  w[abs(x) > k1] <- 0.0

  return(w)
}



.triangular_window <- function(x, k1, k2, ...){
  # Set weights using a triangular window.

  # k1 should be 0 or greater, and k2 cannot be smaller than k1.
  if(k2 < k1) return(NA_real_)
  if(k1 < 0) return(NA_real_)

  # Initialise weights.
  w <- numeric(length(x)) + 1.0

  # Set weights of elements between k1 and k2.
  if(k1 != k2){
    lobe_elements <- which(abs(x) >= k1 & abs(x) <= k2)
    w[lobe_elements] <- 1.0 - (abs(x[lobe_elements]) - k1) / (k2 - k1)
  }

  # Set weights of elements greater than k2.
  w[abs(x) > k2] <- 0.0

  return(w)
}



.tapered_cosine_window <- function(x, k1, k2, ...){
  # Set weights using a tapered cosine window.

  # k1 should be 0 or greater, and k2 cannot be smaller than k1.
  if(k2 < k1) return(NA_real_)
  if(k1 < 0) return(NA_real_)

  # Initialise weights.
  w <- numeric(length(x)) + 1.0

  # Set weights of elements between k1 and k2.
  if(k1 != k2){
    lobe_elements <- which(abs(x) >= k1 & abs(x) <= k2)
    w[lobe_elements] <- 0.5 + 0.5 * cos((abs(x[lobe_elements]) - k1) / (k2 - k1) * pi)
  }

  # Set weights of elements greater than k2.
  w[abs(x) > k2] <- 0.0

  return(w)
}
