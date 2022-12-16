
..box_cox_transform <- function(lambda, x, invert=FALSE){
  # After Box, G. E., & Cox, D. R. (1964). An analysis of transformations.
  # Journal of the Royal Statistical Society. Series B (Methodological),
  # 211-252.

  if(invert){
    # Inverse transformations: From transformed value to original value
    if(lambda==0){
      y <- exp(x)

    } else {
      y <- (x * lambda + 1)^(1/lambda)
    }

  } else {
    # From original value to transformed value

    # Find any non-positive entries and replace them (this may happen in new
    # applications).
    neg_index <- x <= 0 & is.finite(x)
    if(any(neg_index)) x[neg_index] <- min(x[x>0 & is.finite(x)])

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
