#' @include TransformationObjects.R
NULL

# transformationYeoJohnson definition ------------------------------------------

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



# transformationYeoJohnsonShift definition -------------------------------------

#' @rdname transformation_yeo_johnson
#' @export
setClass(
  "transformationYeoJohnsonShift",
  contains="transformationYeoJohnson")



# .set_transformation_parameters (Yeo-Johnson) ---------------------------------
setMethod(
  ".set_transformation_parameters",
  signature("transformationYeoJohnson"),
  function(
    object,
    x,
    lambda,
    estimation_method = "mle",
    weighting_function = NULL,
    weighting_function_parameters = NULL,
    optimiser="subplex",
    backup_use_default=TRUE,
    ...){

    # Set lambda range. If lambda is NULL, set a very wide range.
    if(is.null(lambda)) lambda <- ..get_default_lambda_range()

    # Set lambda, in case a fixed lambda is provided.
    object <- ..set_lambda(
      object = object,
      lambda = lambda)

    # Get optimisation parameters to initialise and configure the optimiser.
    optimisation_parameters <- ..optimisation_parameters(
      object = object,
      x = x,
      lambda = lambda)

    # Skip optimisation if there is nothing to optimise.
    if(is.null(optimisation_parameters)){
      object@complete <- TRUE

      return(object)
    }

    # Initialise the estimator.
    estimator <- .set_estimator(
      transformer = object,
      estimation_method = estimation_method,
      weighting_function = weighting_function,
      weighting_function_parameters = weighting_function_parameters)

    # Optimise transformation parameters.
    optimised_parameters <- .optimise_transformation_parameters(
      object = estimator,
      transformer = object,
      x = x,
      optimiser = optimiser,
      optimisation_parameters = optimisation_parameters)

    # Update shift and lambda values with the optimised parameters.
    if(!is.null(optimised_parameters$shift)){
      if(is.finite(optimised_parameters$shift)){
        object@shift <- optimised_parameters$shift

      } else if(!backup_use_default){
        object@shift <- optimised_parameters$shift
      }
    }

    if(!is.null(optimised_parameters$lambda)){
      if(is.finite(optimised_parameters$lambda)){
        object@lambda <- optimised_parameters$lambda

      } else if(!backup_use_default){
        object@lambda <- optimised_parameters$lambda
      }
    }

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



# ..get_default_shift_range (Yeo-Johnson) --------------------------------------
setMethod(
  "..get_default_shift_range",
  signature(object = "transformationYeoJohnson"),
  function(object, x, ...){
    # Set shift range.
    shift_range <- c(min(x), max(x))

    return(shift_range)
  }
)



# ..get_default_lambda_range (Yeo-Johnson) -------------------------------------
setMethod(
  "..get_default_lambda_range",
  signature(object = "transformationYeoJohnson"),
  function(object, ...){
    return(c(-100.0, 100.0))
  }
)



# ..set_lambda (Yeo-Johnson) ---------------------------------------------------
setMethod(
  "..set_lambda",
  signature(object = "transformationYeoJohnson"),
  function(object, lambda, ...){

    # Only update lambda if it is not a range.
    if(length(lambda) != 1) return(object)

    # Update lambda.
    object@lambda <- lambda

    return(object)
  }
)



# ..optimisation_parameters (Yeo-Johnson) --------------------------------------
setMethod(
  "..optimisation_parameters",
  signature(object = "transformationYeoJohnson"),
  function(object, lambda, ...){

    if(length(lambda) == 1) return(NULL)

    return(
      list(
        "initial"=mean(lambda),
        "lower"=min(lambda),
        "upper"=max(lambda),
        "parameter_type" = "lambda"))
  }
)



# ..optimisation_parameters (Yeo-Johnson (shift)) ------------------------------
setMethod(
  "..optimisation_parameters",
  signature(object = "transformationYeoJohnsonShift"),
  function(object, x, lambda, ...){

    # Set up x-range.
    x_range <- ..get_default_shift_range(
      object = object,
      x = x)

    if(length(lambda) == 1){
      return(
        list(
          "initial" = x_range[1],
          "lower" = x_range[1],
          "upper" = x_range[2],
          "parameter_type" = c("shift")))

    } else {
      return(
        list(
          "initial" = c(x_range[1], mean(lambda)),
          "lower" = c(x_range[1], min(lambda)),
          "upper" = c(x_range[2], max(lambda)),
          "parameter_type" = c("shift", "lambda")))
    }
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



#
# yeo_johnson_shift_range <- function(x){
#   # Default range would be any shift between all-positive (shift by lowest
#   # value) and all-negative (shift by highest value).
#
#   # Find the (negative or zero) minimum value. We need to increment slightly
#   # to avoid x containing 0s.
#   min_value <- min(x, na.rm=TRUE)
#   max_value <- max(x, na.rm=TRUE)
#
#   return(c(min_value, max_value))
# }
#
#
#
# yeo_johnson_parameter_grid <- function(x, lambda_range){
#
#   # Set up grid positions.
#   points_x <- unique(stats::quantile(x, c(0.05, 0.2, 0.35, 0.5, 0.65, 0.8, 0.95), names=FALSE))
#   points_lambda <- lambda_range[1] + (seq_len(11L) - 1.0) * (lambda_range[2] - lambda_range[1]) / 10.0
#
#   # Create parameter pairs that form the grid nodes.
#   parameters <- mapply(
#     function(x, lambda) (c(x, lambda)),
#     x=rep(points_x, each=length(points_lambda)),
#     lambda=rep(points_lambda, times=length(points_x)),
#     SIMPLIFY=FALSE,
#     USE.NAMES=FALSE)
#
#   return(list(
#     "x"=points_x,
#     "x_range"=c(min(points_x), max(points_x)),
#     "lambda"=points_lambda,
#     "lambda_range"=c(min(points_lambda), max(points_lambda)),
#     "parameter"=parameters))
# }
#
#
#
# ..yeo_johnson_transform <- function(lambda, x, invert=FALSE){
#   # After Yeo, I. K., & Johnson, R. A. (2000). A new family of power
#   # transformations to improve normality or symmetry. Biometrika, 87(4),
#   # 954-959.
#
#   # Copy output
#   y <- x
#
#   # Determine positive and negative elements of the input vector
#   pos_index <- x >= 0 & is.finite(x)
#   neg_index <- x < 0 & is.finite(x)
#
#   if(invert) {
#     # Inverse transformations: From transformed value to original value
#     if(any(pos_index)){
#       if(lambda != 0){
#         y[pos_index] <- ((x[pos_index] * lambda + 1)^(1/lambda) - 1)
#
#       } else {
#         y[pos_index] <- exp(x[pos_index]) - 1
#       }
#     }
#
#     if(any(neg_index)){
#       if(lambda != 2) {
#         y[neg_index] <- 1 - (x[neg_index] * (lambda-2) + 1)^(1/(2-lambda))
#
#       } else {
#         y[neg_index] <- 1 - exp(-x[neg_index])
#       }
#     }
#
#   } else {
#
#     # From original value to transformed value
#     if(any(pos_index)){
#       if(lambda == 0.0){
#         y[pos_index] <- log1p(x[pos_index])
#
#       } else {
#         y[pos_index] <- ((x[pos_index] + 1)^lambda - 1) / lambda
#       }
#     }
#
#     if(any(neg_index)){
#       if(lambda == 2.0){
#         y[neg_index] <- -log1p(-x[neg_index])
#
#       } else {
#         y[neg_index] <- -((-x[neg_index] + 1)^(2-lambda) - 1) / (2-lambda)
#       }
#     }
#   }
#
#   return(y)
# }
#
#
#
# ..yeo_johnson_dev <- function(lambda, x){
#   # First order derivative of the Yeo-Johnson transformation with respect to x.
#   return((1 + abs(x))^(sign(x) * (lambda - 1)))
# }
#
#
#
# ..yeo_johnson_loglik <- function(lambda, x, w=NULL){
#
#   # Set w
#   if(is.null(w)) w <- numeric(length(x)) + 1.0
#
#   # Transform x under the provided lambda.
#   y <- ..yeo_johnson_transform(lambda=lambda, x=x)
#
#   # Compute the sum of the weights.
#   sum_w <- sum(w)
#   if(sum_w == 0) return(NA_real_)
#
#   # Compute the weighted estimates of the mean mu and variance sigma squared for
#   # y.
#   mu_hat <- sum(w * y) / sum_w
#   sigma_hat_squared <- sum(w * (y - mu_hat)^2) / sum_w
#
#   # Log-likelihood cannot be estimated if sigma is NaN.
#   if(!is.finite(sigma_hat_squared)) return(NA_real_)
#
#   # Log-likelihood cannot be determined if the sigma estimate equals 0.0
#   if(sigma_hat_squared == 0) return(NA_real_)
#
#   # Compute the log likelihood under the assumption that the transformed
#   # variable y follows the normal distribution.
#   llf <- (lambda - 1.0) * sum(w * sign(x) * log1p(abs(x))) - sum_w/2.0 * log(sigma_hat_squared)
#
#   return(llf)
# }
#
#
#
# ..yeo_johnson_transform_rectified <- function(lambda, x){
#   # The rectified transform replaces part of the transformed values by a first
#   # order (linear) approximation. Linear approximation of a function f(x) at
#   # point a is defined as y = f(a) + (x - a) * f'(a), with f'(a) being the
#   # derivative of f(x=a). Here function f is the Yeo-Johnson transformation, and
#   # point a is the first or third quartile, depending on lambda.
#
#   # Find first and third quartiles.
#   cut_off <- stats::quantile(x, probs=c(0.25, 0.75), names=FALSE)
#
#   y <- numeric(length(x))
#
#   # Perform rectified transformation
#   if(lambda == 1.0){
#     # For lambda equal to 1, the mapping is linear, and no elements are
#     # out-of-range and require rectification.
#     out_of_range <- logical(length(x))
#
#   } else if(lambda > 1.0){
#     # Select the cut-off value, i.e. the first quartile.
#     cut_off <- cut_off[1]
#
#     # Elements that have value below Cl (1st quartile) are rectified.
#     out_of_range <- x < cut_off
#
#   } else {
#     # Lambda < 1.0.
#
#     # Select the cut-off value, i.e. the third quartile.
#     cut_off <- cut_off[2]
#
#     # Elements that have value above Cu (3rd quartile) are rectified.
#     out_of_range <- x > cut_off
#   }
#
#   if(any(out_of_range)){
#     # Linear approximation to out-of-range elements.
#     y[out_of_range] <- ..yeo_johnson_transform(lambda=lambda, x=cut_off) + (x[out_of_range]-cut_off) * ..yeo_johnson_dev(lambda=lambda, x=cut_off)
#   }
#
#   # Map elements that do not require rectification using the normal Box-Cox
#   # transformation.
#   if(any(!out_of_range)){
#     y[!out_of_range] <- ..yeo_johnson_transform(lambda=lambda, x=x[!out_of_range])
#   }
#
#   return(y)
# }
