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





# tukey_tapered_cosine_window <- function(x, k=0.00, ...){
#   # Use Tukey window (also cosine-tapered window) for downweighting based on
#   # expected quantile (not on residual!). The Tukey window has a flat top, with
#   # width k (or 2 alpha), and then tapers off. The alpha parameter is related to
#   # k: alpha = 1 - k.
#
#   if(is.unsorted(x)) stop("DEV: x should be sorted.")
#
#   # Iterate along x
#   i <- seq_along(x)
#
#   # Initialise weights.
#   w <- numeric(length(i)) + 1.0
#
#   # Left tapered cosine lobe.
#   left_lobe <- which(i - 0.5 < length(i) * (1 - k) / 2)
#   w[left_lobe] <- 0.5 * (1 + cos((2 * pi  * (left_lobe - (1 + length(i) * (1 - k)) /2) / (length(i) * (1 - k)))))
#
#   # Right tapered cosine lobe.
#   right_lobe <- which(i - 0.5 > length(i) * (1 + k) /2)
#   w[right_lobe] <- 0.5 * (1 + cos((2 * pi * (right_lobe - (1 + length(i) * (1 + k)) /2) / (length(i) * (1 - k)))))
#
#   return(w)
# }
#
#
# step_window <- function(x, k=0.80, ...){
#   # Use a step window for downweighting based on expected quantile. The window
#   # has a flat top with width k, that drops to 0 outside of the central k
#   # values.
#
#   if(is.unsorted(x)) stop("DEV: x should be sorted.")
#
#   # Iterate along x
#   i <- seq_along(x)
#
#   # Initialise weights.
#   w <- numeric(length(i)) + 1.0
#
#   # Left lobe.
#   left_lobe <- which(i - 0.5 < length(i) * (1 - k) / 2)
#   w[left_lobe] <- 0.0
#
#   # Right lobe.
#   right_lobe <- which(i - 0.5 > length(i) * (1 + k) /2)
#   w[right_lobe] <- 0.0
#
#   return(w)
# }
#
#
#
# tapered_step_window <- function(x, k_step=0.90, k=0.20, ...){
#
#   if(is.unsorted(x)) stop("DEV: x should be sorted.")
#
#   # Iterate along x
#   i <- seq_along(x)
#
#   # Initialise weights.
#   w <- numeric(length(i)) + 1.0
#
#   # Left lobe for step function.
#   left_lobe <- which(i - 0.5 < length(i) * (1 - k_step) / 2)
#   w[left_lobe] <- 0.0
#
#   # Right lobe for step function.
#   right_lobe <- which(i - 0.5 > length(i) * (1 + k_step) /2)
#   w[right_lobe] <- 0.0
#
#   # Left tapered lobe.
#   left_lobe <- which(i - 0.5 < length(i) * (1 - k) / 2 &
#                        i - 0.5 >= length(i) * (1 - k_step) / 2)
#
#   w[left_lobe] <- (left_lobe - min(left_lobe)) / (max(left_lobe) - min(left_lobe))
#
#   # Right tapered lobe.
#   right_lobe <- which(i - 0.5 > length(i) * (1 + k) /2 &
#                         i - 0.5 <= length(i) * (1 + k_step) / 2)
#   w[right_lobe] <- 1.0 - (right_lobe - min(right_lobe)) / (max(right_lobe) - min(right_lobe))
#
#   return(w)
# }
#
#
#
# transformed_step_weighting <- function(x, lambda, type, tau=2.58, ...){
#   # Step weighting used by Rademaekers and Rousseeuw (2021). This basically sets
#   # weights for the outerlying transformed values to 0. Default value is to
#   # qnorm(0.995).
#
#   transform_FUN <- switch(
#     type,
#     "box_cox"=..box_cox_transform,
#     "yeo_johnson"=..yeo_johnson_transform
#   )
#
#   # Find transformed feature values.
#   y <- do.call(
#     transform_FUN,
#     args=list(
#       "lambda"=lambda,
#       "x"=x))
#
#   # Approximate Huber's M-estimates for locality and scale.
#   robust_estimates <- huber_estimate(y, tol=1E-3)
#
#   # Check problematic values.
#   if(!is.finite(robust_estimates$sigma)) return(NA_real_)
#   if(robust_estimates$sigma == 0.0) return(NA_real_)
#
#   # Compute weights.
#   w <- as.numeric(abs(y - robust_estimates$mu) / robust_estimates$sigma <= tau)
#
#   return(w)
# }
#
#
#
# residual_step_weighting <- function(x, lambda, type, tau=0.50, ...){
#   # Step weighting based on residuals.
#
#   # Compute residuals.
#   r <- ..compute_residuals(
#     x=x,
#     lambda=lambda,
#     type=type)
#
#   # Set weights.
#   w <- as.numeric(abs(r) < tau)
#
#   return(w)
# }
#
#
#
# transformed_tukey_biweight <- function(x, lambda, type, tau=2.58, ...){
#   # Set weights based on Tukey's Biweights.
#
#   transform_FUN <- switch(
#     type,
#     "box_cox"=..box_cox_transform,
#     "yeo_johnson"=..yeo_johnson_transform
#   )
#
#   # Find transformed feature values.
#   y <- do.call(
#     transform_FUN,
#     args=list(
#       "lambda"=lambda,
#       "x"=x))
#
#   # Approximate Huber's M-estimates for locality and scale.
#   robust_estimates <- huber_estimate(y, tol=1E-3)
#
#   # Check problematic values.
#   if(!is.finite(robust_estimates$sigma)) return(NA_real_)
#   if(robust_estimates$sigma == 0.0) return(NA_real_)
#
#   # Initialise 0-weights.
#   w <- numeric(length(x))
#
#   # Compute standardised values.
#   z <- (y - robust_estimates$mu) / robust_estimates$sigma
#
#   # Set non-zero weights for z-scores smaller than tau.
#   ii <- which(abs(z) <= tau)
#   if(length(ii) > 0) w[ii] <- (1 - (z[ii] / tau)^2 )^2
#
#   # Compute weights.
#   return(w)
# }
#
#
#
# residual_tukey_biweight <- function(x, lambda, type, tau=0.50, ...){
#   # Set weights based on Tukey's Biweigths.
#
#   # Compute residuals.
#   r <- ..compute_residuals(
#     x=x,
#     lambda=lambda,
#     type=type)
#
#   # Initialise 0-weights.
#   w <- numeric(length(x))
#
#   # Set non-zero weights for residual errors smaller than tau.
#   ii <- which(abs(r) <= tau)
#   if(length(ii) > 0) w[ii] <- (1 - (r[ii] / tau)^2 )^2
#
#   return(w)
# }
#
#
#
#
# transformed_huber_weight <- function(x, lambda, type, tau=1.96, ...){
#   # Set weights based on Tukey's Biweights.
#
#   transform_FUN <- switch(
#     type,
#     "box_cox"=..box_cox_transform,
#     "yeo_johnson"=..yeo_johnson_transform
#   )
#
#   # Find transformed feature values.
#   y <- do.call(
#     transform_FUN,
#     args=list(
#       "lambda"=lambda,
#       "x"=x))
#
#   # Approximate Huber's M-estimates for locality and scale.
#   robust_estimates <- huber_estimate(y, tol=1E-3)
#
#   # Check problematic values.
#   if(!is.finite(robust_estimates$sigma)) return(NA_real_)
#   if(robust_estimates$sigma == 0.0) return(NA_real_)
#
#   # Compute standardised values.
#   z <- (y - robust_estimates$mu) / robust_estimates$sigma
#
#   # Initialise 1-weights.
#   w <- numeric(length(x)) + 1.0
#
#   # Set weights. Standardised values larger than tau are down-weighted.
#   ii <- which(abs(z) > tau)
#   if(length(ii) > 0)  w[ii] <- tau / abs(z[ii])
#
#   # Compute weights.
#   return(w)
# }
#
#
#
# residual_huber_weight <- function(x, lambda, type, tau=0.20, ...){
#   # Set weights based on Hubers weight.
#
#   # Compute residuals.
#   r <- ..compute_residuals(
#     x=x,
#     lambda=lambda,
#     type=type)
#
#   # Initialise 1-weights.
#   w <- numeric(length(x)) + 1.0
#
#   # Set weights. Values > tau are down-weighted.
#   ii <- which(abs(r) > tau)
#   if(length(ii) > 0)  w[ii] <- tau / abs(r[ii])
#
#   return(w)
# }
#
#
#
# set_full_central_weight <- function(x, w, k=0.80){
#   # Assign a full weight to central elements.
#
#   if(is.unsorted(x)) stop("DEV: x should be sorted.")
#
#   # Indices to x.
#   ii <- seq_along(x)
#
#   # Find central elements.
#   central_elements <- which(length(ii) * (1 - k) / 2 <= ii - 0.5 & ii - 0.5 <= length(ii) * (1 + k) /2)
#
#   # Set central weights to 1
#   w[central_elements] <- 1.0
#
#   return(w)
# }
#
#
#
# ..compute_residuals <- function(x, lambda, type){
#
#
# }
