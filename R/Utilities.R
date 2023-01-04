huber_estimate <- function(x, k=1.28, tol=1E-4){
  # k=1.28 is based on Wilcox RR. Introduction to Robust Estimation and
  # Hypothesis Testing. Academic Press; 2011.

  # Filter missing values
  x <- x[is.finite(x)]

  if(length(x) == 0) return(list("mu"=NA_real_, "sigma"=NA_real_))

  # Initial estimates for the estimate mu and scale sigma
  mu_0 <- stats::median(x)
  sigma_0 <- stats::mad(x)

  # Check that there is an initial estimate for scale.
  if(sigma_0 == 0.0) return(list("mu"=mu_0, "sigma"=0.0))

  # Determine theta, i.e. the probability corresponding to k, according to a
  # Gaussian distribution.
  theta <- 2.0 * stats::pnorm(k) - 1.0
  beta <- theta + k^2 * (1 - theta) - 2.0 * k * stats::dnorm(k)

  for(ii in seq_len(50)){
    # Apply Huber's rho, that basically winsorizes values |x| >= k * sigma_0
    xx <- pmin(pmax(mu_0 - k * sigma_0, x), mu_0 + k * sigma_0)

    # Compute updated mean and scale estimate.
    mu_1 <- sum(xx) / length(xx)
    sigma_1 <- sqrt(sum((xx - mu_1)^2) / (beta * (length(xx) - 1)))

    # Check convergence.
    if(abs(mu_0 - mu_1) < tol * mu_0 && abs(sigma_0 - sigma_1) < tol * sigma_0) break

    # Values for next iteration
    mu_0 <- mu_1
    sigma_0 <- sigma_1
  }

  return(list("mu"=mu_0, "sigma"=sigma_0))
}



compute_expected_z <- function(x){

  # Check if x is sorted, and sort otherwise.
  if(is.unsorted(x)) stop(paste0("DEV: x is expected to be sorted in ascending order."))

  # Compute expected quantile.
  q <- (seq_along(x) - 1/3) / (length(x) + 1/3)

  # Set up a data.table.
  data <- data.table::data.table(
    "x"=x,
    "q"=q)

  # Average quantile for when x has multiple values. Though this necessitates
  # using the data.table package, this is by far the fastest implementation.
  data[, "q_group":=mean(q), by="x"]

  # Compute z-scores.
  z <- stats::qnorm(p=data$q_group)

  return(z)
}



apply_boundary <- function(
    llf,
    shift,
    shift_range,
    lambda,
    lambda_range){
  # Apply boundary conditions to steer optimisation algorithms away from the
  # boundaries.

  if(shift < shift_range[1] || shift > shift_range[2] || lambda < lambda_range[1] || lambda > lambda_range[2]){
    # Build a gradient that increases the further shift and/or lambda parameters
    # deviate from the boundary.

    return(NA_real_)
  }

  return(llf)
}



select_neighbourhood <- function(x, x_range){
  neighbourhood_range <- c(0.0, 0.0)

  # Find coordinates.
  ii <- which(x_range == x)

  # Lower bound.
  if(ii == 1){
    neighbourhood_range[1] <- x

  } else {
    neighbourhood_range[1] <- x_range[ii - 1L]
  }

  # Upper bound.
  if(ii == length(x_range)){
    neighbourhood_range[2] <- x

  } else {
    neighbourhood_range[2] <- x_range[ii + 1L]
  }

  return(neighbourhood_range)
}



compute_expected_z <- function(x){

  # Check if x is sorted, and sort otherwise.
  if(is.unsorted(x)) stop(paste0("DEV: x is expected to be sorted in ascending order."))

  # Compute expected quantile.
  q <- (seq_along(x) - 1/3) / (length(x) + 1/3)

  # Set up a data.table.
  data <- data.table::data.table(
    "x"=x,
    "q"=q)

  # Average quantile for when x has multiple values. Though this necessitates
  # using the data.table package, this is by far the fastest implementation.
  data[, "q_group":=mean(q), by="x"]

  # Compute z-scores.
  z <- stats::qnorm(p=data$q_group)

  return(z)
}
