box_cox_shift_range <- function(x){
  # Default range would be any shift between all-positive, limit on 0 (shift by lowest
  # value) and default (0).

  # Find the value that brings the entire distribution to 0. We need to
  # increment slightly to avoid x containing 0s.
  max_value <- min(x, na.rm=TRUE)
  max_value_offset <- 1.0
  min_value_offset <- 1.0

  # Find the typical, non-zero, distance between values. NA values should be
  # removed.
  dx <- unique(diff(sort(x, na.last=NA)))
  dx <- dx[dx > 0.0]

  # It shouldn't happen that all values are the same - but better check it.
  if(length(dx) > 0) max_value_offset <- stats::median(dx)

  # The increment should not grow too much.
  if(max_value_offset > 0.5) max_value_offset <- 0.5

  # Update min_value_offset.
  min_value_offset <- stats::median(x - max_value, na.rm=TRUE)
  if(min_value_offset < 1.0) min_value_offset <- 1.0

  # Set minimum shift value.
  min_value <- max_value - min_value_offset

  # Set maximum shift value.
  max_value <- max_value - max_value_offset

  # Set shift range. Occasionally, the values not be sorted.
  shift_range <- sort(c(min_value, max_value))

  return(shift_range)
}



box_cox_parameter_grid <- function(x, lambda_range){

  # Set up x-range.
  x_range <- box_cox_shift_range(x)

  # Set up grid positions.
  points_x <- x_range[1] + (seq_len(5)-1.0) * (x_range[2] - x_range[1]) / 4
  points_lambda <- lambda_range[1] + (seq_len(11L) - 1.0) * (lambda_range[2] - lambda_range[1]) / 10.0

  # Create parameter pairs that form the grid nodes.
  parameters <- mapply(
    function(x, lambda) (c(x, lambda)),
    x=rep(points_x, each=length(points_lambda)),
    lambda=rep(points_lambda, times=length(points_x)),
    SIMPLIFY=FALSE,
    USE.NAMES=FALSE)

  return(list(
    "x"=points_x,
    "x_range"=c(min(points_x), max(points_x)),
    "lambda"=points_lambda,
    "lambda_range"=c(min(points_lambda), max(points_lambda)),
    "parameter"=parameters))
}



..box_cox_transform <- function(lambda, x, invert=FALSE){
  # After Box, G. E., & Cox, D. R. (1964). An analysis of transformations.
  # Journal of the Royal Statistical Society. Series B (Methodological),
  # 211-252.

  if(invert){
    # Inverse transformations: From transformed value to original value.
    if(lambda==0){
      y <- exp(x)

    } else {
      y <- (x * lambda + 1)^(1/lambda)
    }

  } else {
    # From original value to transformed value.
    if(lambda==0){
      y <- log(x)

    } else {
      y <- (x^lambda - 1) / lambda
    }
  }

  return(y)
}



..box_cox_dev <- function(lambda, x){
  # First order derivative of the Yeo-Johnson transformation with respect to x.
  return(x^(lambda - 1))
}



..box_cox_loglik <- function(lambda, x, w=NULL){

  # Set w
  if(is.null(w)) w <- numeric(length(x)) + 1.0

  # Transform x under the provided lambda.
  y <- ..box_cox_transform(lambda=lambda, x=x)

  if(any(!is.finite(y))) return(NA_real_)

  # Compute the sum of the weights.
  sum_w <- sum(w)
  if(sum_w == 0) return(NA_real_)

  # Compute the weighted estimates of the mean mu and variance sigma squared for
  # y.
  mu_hat <- sum(w * y) / sum_w
  sigma_hat_squared <- sum(w * (y - mu_hat)^2) / sum_w

  # Log-likelihood cannot be determined if the sigma estimate equals 0.0
  if(sigma_hat_squared == 0) return(NA_real_)

  # Compute the log likelihood under the assumption that the transformed
  # variable y follows the normal distribution.
  llf <- (lambda - 1.0) * sum(w * log(x)) - sum_w / 2.0 * log(sigma_hat_squared)

  return(llf)
}



..box_cox_transform_rectified <- function(lambda, x){
  # The rectified transform replaces part of the transformed values by a first
  # order (linear) approximation. Linear approximation of a function f(x) at
  # point a is defined as y = f(a) + (x - a) * f'(a), with f'(a) being the
  # derivative of f(x=a). Here function f is the Box-Cox transformation, and
  # point a is the first or third quartile, depending on lambda.

  # Find first and third quartiles.
  cut_off <- stats::quantile(x, probs=c(0.25, 0.75), names=FALSE)

  y <- numeric(length(x))

  # Perform rectified transformation
  if(lambda == 1.0){
    # For lambda equal to 1, the mapping is linear, and no elements are
    # out-of-range and require rectification.
    out_of_range <- logical(length(x))

  } else if(lambda > 1.0){
    # Select the cut-off value, i.e. the first quartile.
    cut_off <- cut_off[1]

    # Elements that have value below Cl (1st quartile) are rectified.
    out_of_range <- x < cut_off

  } else {
    # Lambda < 1.0.

    # Select the cut-off value, i.e. the third quartile.
    cut_off <- cut_off[2]

    # Elements that have value above Cu (3rd quartile) are rectified.
    out_of_range <- x > cut_off
  }

  if(any(out_of_range)){
    # Linear approximation to out-of-range elements.
    y[out_of_range] <- ..box_cox_transform(lambda=lambda, x=cut_off) + (x[out_of_range]-cut_off) * ..box_cox_dev(lambda=lambda, x=cut_off)
  }

  # Map elements that do not require rectification using the normal Box-Cox
  # transformation.
  if(any(!out_of_range)){
    y[!out_of_range] <- ..box_cox_transform(lambda=lambda, x=x[!out_of_range])
  }

  return(y)
}
