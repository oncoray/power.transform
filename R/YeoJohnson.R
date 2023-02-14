#' @include TransformationObjects.R
NULL

# Yeo-Johnson transformation class ---------------------------------------------

#' Yeo-Johnson transformation object
#'
#' This class is used for Yeo-Johnson transformations.
#'
#' @slot method Main transformation method, i.e. `"yeo_johnson"`.
#' @slot robust Indicates whether a robust version of the Yeo-Johnson
#'   transformation is used to set transformation parameters. The value depends
#'   on the `robust` argument of the `find_transformation_parameters` function.
#' @slot lambda Numeric lambda parameter for the Yeo-Johnson transformation.
#' @slot shift Numeric shift parameter for the Yeo-Johnson transformation.If
#'   `shift=TRUE` in the `find_transformation_parameters` function, `lambda` and
#'   `shift` parameters are optimised simultaneously. Otherwise, the `shift`
#'   parameter has a value of `0.0`.
#' @slot complete Indicates whether transformation parameters were set.
#'
#' @seealso [find_transformation_parameters]
#' @export
#' @rdname transformation_yeo_johnson

setClass(
  "transformationYeoJohnson",
  contains="transformationPowerTransform",
  slots=list(
    "method" = "character",
    "robust" = "logical",
    "lambda" = "numeric",
    "shift" = "numeric"),
  prototype=list(
    "method" = "yeo_johnson",
    "robust" = FALSE,
    "lambda" = NA_real_,
    "shift" = 0.0))

#' @rdname transformation_yeo_johnson
#' @export
setClass(
  "transformationYeoJohnsonShift",
  contains="transformationYeoJohnson")



# .set_transformation_parameters (Yeo-Johnson) ---------------------------------
setMethod(
  ".set_transformation_parameters",
  signature("transformationYeoJohnson"),
  function(object, x, lambda, ...){

    if(is.null(lambda)){
      # Pick very wide lambda-range, if lambda is NULL.
      lambda_range <- c(-100.0, 100.0)

    } else if(length(lambda) == 2){
      lambda_range <- lambda
    }

    if(length(lambda) == 1){
      # Set lambda parameter. Since we don't use any robust algorithm, set
      # robust to FALSE.
      object@lambda <- lambda
      object@robust <- FALSE
      object@complete <- TRUE

    } else if(object@robust){
      # Robust method based on Raymaekers J, Rousseeuw PJ. Transforming
      # variables to central normality. Mach Learn. 2021.
      # doi:10.1007/s10994-021-05960-5
      lambda <- .transformation_robust_optimisation(
        x = x,
        type = "yeo_johnson",
        lambda_range = lambda_range)

      # Set lambda parameter.
      object@lambda <- lambda
      object@complete <- TRUE

    } else {
      # Standard method based on optimising log-likelihood of the normal
      # distribution.
      optimal_lambda <- suppressWarnings(
        stats::optimise(
          ..yeo_johnson_loglik,
          interval=lambda_range,
          x=x,
          maximum=TRUE))

      lambda <- ifelse(
        is.finite(optimal_lambda$objective),
        optimal_lambda$maximum,
        1.0)

      # Set lambda parameter.
      object@lambda <- lambda
      object@complete <- TRUE
    }

    return(object)
  }
)


# .set_transformation_parameters (Yeo-Johnson (shift)) -------------------------
setMethod(
  ".set_transformation_parameters",
  signature("transformationYeoJohnsonShift"),
  function(object, x, lambda, optimiser="subplex", backup_use_default=TRUE, ...){

    if(is.null(lambda)){
      # Pick very wide lambda-range, if lambda is NULL.
      lambda_range <- c(-100.0, 100.0)

    } else if(length(lambda) == 2){
      lambda_range <- lambda

    } else {
      return(methods::callNextMethod())
    }

    # Get optimisation parameters to initialise and configure the optimiser.
    optimisation_parameters <- ..optimisation_parameters(
      object=object,
      x = x,
      lambda = lambda_range)

    if(object@robust){
      # Sort values
      x <- sort(x)

      # Optimisation function.
      opt_fun <- function(...) -.transformation_robust_shifted_optimisation(...)

    } else {
      # Optimisation function.
      # opt_fun <- .transform_shifted_optimisation
      opt_fun <- function(...) -.transform_shifted_optimisation(...)
    }

    if(!is_package_installed("nloptr") & optimiser %in% c("direct-l", "subplex", "nelder-mead")){
      warning(paste0(
        "The nloptr package is required to optimise power transformation parameters using the ",
        optimiser, " algoritm. stats::optim is used as a fallback option."))

      optimiser <- "optim-nelder-mead"
    }

    if(optimiser == "direct-l"){
      # DIRECT-L algorithm
      #
      # D. R. Jones, C. D. Perttunen, and B. E. Stuckmann, “Lipschitzian
      # optimization without the lipschitz constant,” J. Optimization Theory and
      # Applications, vol. 79, p. 157 (1993).
      #
      # J. M. Gablonsky and C. T. Kelley, “A locally-biased form of the DIRECT
      # algorithm," J. Global Optimization, vol. 21 (1), p. 27-37 (2001).

      results <- nloptr::directL(
        fn=opt_fun,
        lower=optimisation_parameters$lower,
        upper=optimisation_parameters$upper,
        control=list("xtol_rel"=1e-3, ftol_rel=1e-4),
        x=x,
        ...,
        type="yeo_johnson")

    } else if(optimiser == "subplex"){
      # SUBPLEX algorithm
      #
      # T. Rowan, “Functional Stability Analysis of Numerical
      # Algorithms”, Ph.D. thesis, Department of Computer Sciences, University
      # of Texas at Austin, 1990.

      results <- tryCatch(
        nloptr::sbplx(
          x0=optimisation_parameters$initial,
          fn=opt_fun,
          lower=optimisation_parameters$lower,
          upper=optimisation_parameters$upper,
          control=list("xtol_rel"=1e-3, ftol_rel=1e-4),
          x=x,
          ...,
          type="yeo_johnson"),
        error = identity)

    } else if(optimiser == "nelder-mead"){
      # Nelder-Mead simplex algorithm
      #
      # J. A. Nelder and R. Mead, “A simplex method for function minimization,”
      # The Computer Journal 7, p. 308-313 (1965).
      #
      # M. J. Box, “A new method of constrained optimization and a comparison
      # with other methods,” Computer J. 8 (1), 42-52 (1965).

      results <- nloptr::neldermead(
        x0=optimisation_parameters$initial,
        fn=opt_fun,
        lower=optimisation_parameters$lower,
        upper=optimisation_parameters$upper,
        control=list("xtol_rel"=1e-3, ftol_rel=1e-4),
        x=x,
        ...,
        type="yeo_johnson")

    } else if(optimiser == "optim-nelder-mead"){
      # Fall-back optimiser in case nloptr is not available. The Nelder-Mead
      # algorithm in stats::optim does not yield results as consistent as the
      # nloptr optimisers.
      results <- stats::optim(
        par=optimisation_parameters$initial,
        fn=opt_fun,
        gr=NULL,
        x=x,
        ...,
        type="yeo_johnson",
        control=list(
          "abstol"=1E-5,
          "reltol"=1E-5))

    } else {
      stop(paste0("Optimiser not recognised: ", optimiser))
    }

    if(inherits(results, "error")){
      if(results$message == "objective in x0 returns NA"){
        shift <- NA_real_
        lambda <- NA_real_

      } else {
        stop(results)
      }

    } else {
      # Extract optimal values.
      shift <- results$par[1]
      lambda <- results$par[2]
    }

    if(backup_use_default && (!is.finite(results$value) || !is.finite(shift) || !is.finite(lambda))){
      shift <- 0.0
      lambda <- 1.0
    }

    # Set parameters.
    object@shift <- shift
    object@lambda <- lambda
    object@complete <- TRUE

    return(object)
  }
)



# .transform (Yeo-Johnson) -----------------------------------------------------
setMethod(
  ".transform",
  signature(object="transformationYeoJohnson"),
  function(object, x, ...){

    # Set up vector to copy new values to.
    y <- numeric(length(x))

    # Find non-finite instances.
    na_entries <- which(!is.finite(x))
    if(length(na_entries) > 0){
      warning("One or more NA or infinite values were found.")

      y[na_entries] <- NA_real_
    }

    # Find remaining valid instances.
    valid_entries <- setdiff(seq_along(x), na_entries)

    # Perform transformation.
    if(length(valid_entries) > 0) y[valid_entries] <- ..transform(
      object = object,
      x = x[valid_entries])

    return(y)
  })



# ..transform (Yeo-Johnson) ----------------------------------------------------
setMethod(
  "..transform",
  signature(object = "transformationYeoJohnson"),
  function(object, x, ...){

    # Subtract shift.
    y <- x - object@shift

    # Determine positive and negative elements of the input vector
    pos_index <- y >= 0
    neg_index <- y < 0

    if(any(pos_index)){
      if(object@lambda == 0.0){
        y[pos_index] <- log1p(x[pos_index])

      } else {
        y[pos_index] <- ((x[pos_index] + 1)^object@lambda - 1) / object@lambda
      }
    }

    if(any(neg_index)){
      if(object@lambda == 2.0){
        y[neg_index] <- -log1p(-x[neg_index])

      } else {
        y[neg_index] <- -((-x[neg_index] + 1)^(2 - object@lambda) - 1) / (2 - object@lambda)
      }
    }

    return(y)
  }
)



# .revert_transform (Yeo-Johnson) ----------------------------------------------
setMethod(
  ".revert_transform",
  signature(object="transformationYeoJohnson"),
  function(object, x, ...){

    # Copy output
    y <- x

    # Determine positive and negative elements of the input vector
    pos_index <- x >= 0 & is.finite(x)
    neg_index <- x < 0 & is.finite(x)

    if(any(pos_index)){
      if(object@lambda != 0){
        y[pos_index] <- ((x[pos_index] * object@lambda + 1)^(1 / object@lambda) - 1)

      } else {
        y[pos_index] <- exp(x[pos_index]) - 1
      }
    }

    if(any(neg_index)){
      if(object@lambda != 2) {
        y[neg_index] <- 1 - (x[neg_index] * (object@lambda - 2) + 1)^(1 / (2 - object@lambda))

      } else {
        y[neg_index] <- 1 - exp(-x[neg_index])
      }
    }

    return(x)
  })



# ..optimisation_parameters (Yeo-Johnson) --------------------------------------
setMethod(
  "..optimisation_parameters",
  signature(object = "transformationYeoJohnson"),
  function(object, lambda, ...){

    return(
      list(
        "initial"=mean(lambda),
        "lower"=min(lambda),
        "upper"=max(lambda)))
  }
)


# ..optimisation_parameters (Yeo-Johnson (shift)) ------------------------------
setMethod(
  "..optimisation_parameters",
  signature(object = "transformationYeoJohnsonShift"),
  function(object, x, lambda, ...){

    # Find the (negative or zero) minimum value. We need to increment slightly
    # to avoid x containing 0s.
    min_value <- min(x, na.rm=TRUE)
    max_value <- max(x, na.rm=TRUE)

    return(
      list(
        "initial"=c(min_value, mean(lambda)),
        "lower"=c(min_value, min(lambda)),
        "upper"=c(max_value, max(lambda))))
  }
)



# ..log_likelihood (Yeo-Johnson) -----------------------------------------------
setMethod(
  "..log_likelihood",
  signature(object = "transformationYeoJohnson"),
  function(
    object,
    x,
    w,
    sigma_hat_squared,
    ...){

    # Compute the log likelihood under the assumption that the transformed
    # variable y follows the normal distribution.
    return(object@lambda - 1.0) * sum(w * sign(x) * log1p(abs(x))) - sum(w)/2.0 * log(sigma_hat_squared)
  }
)



# ..first_derivative (Yeo-Johnson) ---------------------------------------------
setMethod(
  "..first_derivative",
  signature(object = "transformationYeoJohnson"),
  function(object, x){
    return((1 + abs(x))^(sign(x) * (object@lambda - 1)))
  }
)



# ..get_available_estimators (Yeo-Johnson) -------------------------------------
setMethod(
  "..get_available_estimators",
  signature(object = "transformationYeoJohnson"),
  function(object, ...){

    available_estimators <- ..estimators_all()

    # Only allow Raymaekers and Rousseeuw's method for robust optimisation.
    if(!object@robust) available_estimators <- setdiff(available_estimators, ..estimators_raymaekers_robust())

    return(available_estimators)
  }
)



# ..get_available_estimators (Yeo-Johnson (shift)) -----------------------------
setMethod(
  "..get_available_estimators",
  signature(object = "transformationYeoJohnsonShift"),
  function(object, ...){

    available_estimators <- ..estimators_all()

    # Raymaekers and Rousseeuw's method for robust optimisation is not suited
    # for simultaneous estimation of shift and lambda parameters.
    available_estimators <- setdiff(available_estimators, ..estimators_raymaekers_robust())

    return(available_estimators)
  }
)




# show (Yeo-Johnson) -----------------------------------------------------------
setMethod(
  "show",
  signature("transformationYeoJohnson"),
  function(object){

    str <- paste0(
      "A ", ifelse(object@robust, "robust ", ""),
      ifelse(object@shift != 0.0, "shifted ", ""),
      "Yeo-Johnson transformation object"
    )

    if(object@complete){
      cat(paste0(str, ".\n"))
      cat("  lambda: ", object@lambda, "\n")

      if(object@shift != 0.0) cat(paste0("  shift: ", object@shift, "\n"))

    } else {
      cat(paste0(str, " with unset transformation parameters.\n"))
    }
  }
)



# show (Yeo-Johnson (shift)) ---------------------------------------------------
setMethod(
  "show",
  signature("transformationYeoJohnsonShift"),
  function(object){

    str <- paste0(
      "A ", ifelse(object@robust, "robust ", ""),
      ifelse(object@shift != 0.0, "shifted ", ""),
      "Yeo-Johnson transformation object"
    )

    if(object@complete){
      cat(paste0(str, ".\n"))
      cat("  lambda: ", object@lambda, "\n")
      cat("  shift: ", object@shift, "\n")

    } else {
      cat(paste0(str, " with unset transformation parameters.\n"))
    }
  }
)




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

  # Log-likelihood cannot be estimated if sigma is NaN.
  if(!is.finite(sigma_hat_squared)) return(NA_real_)

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
