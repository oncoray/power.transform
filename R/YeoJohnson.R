yeo_johnson_shift_range <- function(x){
  # Default range would be any shift between all-positive (shift by lowest
  # value) and all-negative (shift by highest value).

  # Find the (negative or zero) minimum value. We need to increment slightly
  # to avoid x containing 0s.
  min_value <- min(x, na.rm=TRUE)
  max_value <- max(x, na.rm=TRUE)

  return(c(min_value, max_value))
}



yeo_johnson_parameter_grid <- function(x, lambda_range){

  # Set up grid positions.
  points_x <- unique(stats::quantile(x, c(0.05, 0.2, 0.35, 0.5, 0.65, 0.8, 0.95), names=FALSE))
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



..yeo_johnson_transform <- function(lambda, x, invert=FALSE){
  # After Yeo, I. K., & Johnson, R. A. (2000). A new family of power
  # transformations to improve normality or symmetry. Biometrika, 87(4),
  # 954-959.

  # Copy output
  y <- x

  # Determine positive and negative elements of the input vector
  pos_index <- x >= 0 & is.finite(x)
  neg_index <- x < 0 & is.finite(x)

  if(invert) {
    # Inverse transformations: From transformed value to original value
    if(any(pos_index)){
      if(lambda != 0){
        y[pos_index] <- ((x[pos_index] * lambda + 1)^(1/lambda) - 1)

      } else {
        y[pos_index] <- exp(x[pos_index]) - 1
      }
    }

    if(any(neg_index)){
      if(lambda != 2) {
        y[neg_index] <- 1 - (x[neg_index] * (lambda-2) + 1)^(1/(2-lambda))

      } else {
        y[neg_index] <- 1 - exp(-x[neg_index])
      }
    }

  } else {

    # From original value to transformed value
    if(any(pos_index)){
      if(lambda == 0.0){
        y[pos_index] <- log1p(x[pos_index])

      } else {
        y[pos_index] <- ((x[pos_index] + 1)^lambda - 1) / lambda
      }
    }

    if(any(neg_index)){
      if(lambda == 2.0){
        y[neg_index] <- -log1p(-x[neg_index])

      } else {
        y[neg_index] <- -((-x[neg_index] + 1)^(2-lambda) - 1) / (2-lambda)
      }
    }
  }

  return(y)
}



..yeo_johnson_dev <- function(lambda, x){
  # First order derivative of the Yeo-Johnson transformation with respect to x.
  return((1 + abs(x))^(sign(x) * (lambda - 1)))
}



..yeo_johnson_loglik <- function(lambda, x, w=NULL){

  # Set w
  if(is.null(w)) w <- numeric(length(x)) + 1.0

  # Transform x under the provided lambda.
  y <- ..yeo_johnson_transform(lambda=lambda, x=x)

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
  llf <- (lambda - 1.0) * sum(w * sign(x) * log1p(abs(x))) - sum_w/2.0 * log(sigma_hat_squared)

  return(llf)
}



..yeo_johnson_transform_rectified <- function(lambda, x){
  # The rectified transform replaces part of the transformed values by a first
  # order (linear) approximation. Linear approximation of a function f(x) at
  # point a is defined as y = f(a) + (x - a) * f'(a), with f'(a) being the
  # derivative of f(x=a). Here function f is the Yeo-Johnson transformation, and
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
    y[out_of_range] <- ..yeo_johnson_transform(lambda=lambda, x=cut_off) + (x[out_of_range]-cut_off) * ..yeo_johnson_dev(lambda=lambda, x=cut_off)
  }

  # Map elements that do not require rectification using the normal Box-Cox
  # transformation.
  if(any(!out_of_range)){
    y[!out_of_range] <- ..yeo_johnson_transform(lambda=lambda, x=x[!out_of_range])
  }

  return(y)
}
