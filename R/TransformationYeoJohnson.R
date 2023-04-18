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
  contains = "transformationPowerTransform",
  slots = list(
    "method" = "character",
    "robust" = "logical",
    "lambda" = "numeric",
    "shift" = "numeric"),
  prototype = list(
    "method" = "yeo_johnson",
    "robust" = FALSE,
    "lambda" = NA_real_,
    "shift" = 0.0))



# transformationYeoJohnsonShift definition -------------------------------------

#' @rdname transformation_yeo_johnson
#' @export
setClass(
  "transformationYeoJohnsonShift",
  contains = "transformationYeoJohnson")



# .set_transformation_parameters (Yeo-Johnson) ---------------------------------
setMethod(
  ".set_transformation_parameters",
  signature(object = "transformationYeoJohnson"),
  function(
    object,
    x,
    lambda,
    estimation_method = "mle",
    weighting_function = NULL,
    weighting_function_parameters = NULL,
    optimiser = NULL,
    backup_use_default = TRUE,
    ...) {

    # Set lambda range. If lambda is NULL, set a very wide range.
    if (is.null(lambda)) lambda <- ..get_default_lambda_range(object = object)

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
    if (is.null(optimisation_parameters)) {
      # Add package version.
      object <- .set_version(object = object)

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
      optimisation_parameters = optimisation_parameters,
      ...)

    if (!is.finite(optimised_parameters$lambda) & estimation_method != "mle") {
      # Fallback option in case a non-MLE method fails.
      estimator <- .set_estimator(
        transformer = object,
        estimation_method = "mle",
        weighting_function = weighting_function,
        weighting_function_parameters = weighting_function_parameters)

      # Optimise transformation parameters.
      optimised_parameters <- .optimise_transformation_parameters(
        object = estimator,
        transformer = object,
        x = x,
        optimiser = optimiser,
        optimisation_parameters = optimisation_parameters,
        ...)
    }

    # Update shift and lambda values with the optimised parameters.
    if (!is.null(optimised_parameters$shift)) {
      if (is.finite(optimised_parameters$shift)) {
        object@shift <- optimised_parameters$shift

      } else if (!backup_use_default) {
        object@shift <- optimised_parameters$shift

      } else {
        object@shift <- 0.0
      }
    }

    if (!is.null(optimised_parameters$lambda)) {
      if (is.finite(optimised_parameters$lambda)) {
        object@lambda <- optimised_parameters$lambda

      } else if (!backup_use_default) {
        object@lambda <- optimised_parameters$lambda

      } else {
        object@lambda <- 1.0
      }
    }

    # Add package version.
    object <- .set_version(object = object)

    object@complete <- TRUE

    return(object)
  }
)



# .transform (Yeo-Johnson) -----------------------------------------------------
setMethod(
  ".transform",
  signature(object = "transformationYeoJohnson"),
  function(object, x, ...) {

    # Set up vector to copy new values to.
    y <- numeric(length(x))

    # Find non-finite instances.
    na_entries <- which(!is.finite(x))
    if (length(na_entries) > 0) {
      rlang::warn(
        message = paste0(
          "Yeo-Johnson transforms are only defined for finite values. ",
          length(na_entries), " NA or infinite values were found."),
        class = "power_transform_transform_invalid_values"
      )

      y[na_entries] <- NA_real_
    }

    # Find remaining valid instances.
    valid_entries <- setdiff(seq_along(x), na_entries)

    # Perform transformation.
    if (length(valid_entries) > 0) y[valid_entries] <- ..transform(
      object = object,
      x = x[valid_entries])

    return(y)
  }
)



# ..transform (Yeo-Johnson) ----------------------------------------------------
setMethod(
  "..transform",
  signature(object = "transformationYeoJohnson"),
  function(object, x, ...) {

    # Subtract shift.
    x <- x - object@shift

    # Determine positive and negative elements of the input vector
    pos_index <- x >= 0
    neg_index <- x < 0

    # Initialise y.
    y <- rep(NA_real_, length(x))

    if (any(pos_index)) {
      if (object@lambda == 0.0) {
        y[pos_index] <- log1p(x[pos_index])

      } else {
        y[pos_index] <- ((x[pos_index] + 1)^object@lambda - 1) / object@lambda
      }
    }

    if (any(neg_index)) {
      if (object@lambda == 2.0) {
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
  signature(object = "transformationYeoJohnson"),
  function(object, x, ...) {

    # Copy output
    y <- x

    # Determine positive and negative elements of the input vector
    pos_index <- x >= 0 & is.finite(x)
    neg_index <- x < 0 & is.finite(x)

    if (any(pos_index)) {
      if (object@lambda != 0) {
        y[pos_index] <- ((x[pos_index] * object@lambda + 1)^(1 / object@lambda) - 1)

      } else {
        y[pos_index] <- exp(x[pos_index]) - 1
      }
    }

    if (any(neg_index)) {
      if (object@lambda != 2) {
        y[neg_index] <- 1 - (x[neg_index] * (object@lambda - 2) + 1)^(1 / (2 - object@lambda))

      } else {
        y[neg_index] <- 1 - exp(-x[neg_index])
      }
    }

    # Apply shift.
    y <- y + object@shift

    return(y)
  })



# ..requires_shift_optimisation (Yeo-Johnson (shift)) --------------------------
setMethod(
  "..requires_shift_optimisation",
  signature(object = "transformationYeoJohnsonShift"),
  function(object, ...) {
    return(TRUE)
  }
)



# ..get_default_shift_range (Yeo-Johnson) --------------------------------------
setMethod(
  "..get_default_shift_range",
  signature(object = "transformationYeoJohnson"),
  function(object, x, ...) {
    # Set shift range.
    shift_range <- c(min(x), max(x))

    return(shift_range)
  }
)



# ..get_default_lambda_range (Yeo-Johnson) -------------------------------------
setMethod(
  "..get_default_lambda_range",
  signature(object = "transformationYeoJohnson"),
  function(object, ...) {
    return(c(-100.0, 100.0))
  }
)



# ..set_lambda (Yeo-Johnson) ---------------------------------------------------
setMethod(
  "..set_lambda",
  signature(object = "transformationYeoJohnson"),
  function(object, lambda, ...) {

    # Only update lambda if it is not a range.
    if (length(lambda) != 1) return(object)

    # Update lambda.
    object@lambda <- lambda

    return(object)
  }
)



# ..optimisation_parameters (Yeo-Johnson) --------------------------------------
setMethod(
  "..optimisation_parameters",
  signature(object = "transformationYeoJohnson"),
  function(object, lambda, ...) {

    if (length(lambda) == 1) return(NULL)

    return(
      list(
        "initial" = mean(lambda),
        "lower" = min(lambda),
        "upper" = max(lambda),
        "parameter_type" = "lambda"))
  }
)



# ..optimisation_parameters (Yeo-Johnson (shift)) ------------------------------
setMethod(
  "..optimisation_parameters",
  signature(object = "transformationYeoJohnsonShift"),
  function(object, x, lambda, ...) {

    # Set up x-range.
    x_range <- ..get_default_shift_range(
      object = object,
      x = x)

    # Compute skewness. This prevent local optimisers from becoming stuck in a
    # local optimum at the edge of the search grid.
    mu <- sum(x) / length(x)
    sigma_squared <- sum((x - mu)^2) / length(x)

    # Compute (weighted) skewness and kurtosis.
    skewness <- (1.0 / sigma_squared^(3 / 2)) * (sum((x - mu)^3)) / length(x)
    if (is.na(skewness)) skewness <- 0.0

    if (skewness < 0.0) {
      x_initial <- stats::quantile(x, 0.95)

    } else {
      x_initial <- stats::quantile(x, 0.05)
    }

    if (length(lambda) == 1) {
      return(
        list(
          "initial" = x_initial,
          "lower" = x_range[1],
          "upper" = x_range[2],
          "parameter_type" = c("shift")))

    } else {
      return(
        list(
          "initial" = c(x_initial, mean(lambda)),
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
    ...) {

    # Ensure that data is shifted
    x <- x - object@shift

    # Compute the log likelihood under the assumption that the transformed
    # variable y follows the normal distribution.
    return((object@lambda - 1.0) * sum(w * sign(x) * log1p(abs(x))) - sum(w) / 2.0 * log(sigma_hat_squared))
  }
)



# ..first_derivative (Yeo-Johnson) ---------------------------------------------
setMethod(
  "..first_derivative",
  signature(object = "transformationYeoJohnson"),
  function(object, x) {
    return((1 + abs(x))^(sign(x) * (object@lambda - 1)))
  }
)



# ..get_available_estimators (Yeo-Johnson) -------------------------------------
setMethod(
  "..get_available_estimators",
  signature(object = "transformationYeoJohnson"),
  function(object, ...) {

    available_estimators <- ..estimators_all()

    # Only allow Raymaekers and Rousseeuw's method for robust optimisation.
    if (!object@robust) {
      available_estimators <- setdiff(available_estimators, ..estimators_raymaekers_robust())
    }

    return(available_estimators)
  }
)



# ..get_available_estimators (Yeo-Johnson (shift)) -----------------------------
setMethod(
  "..get_available_estimators",
  signature(object = "transformationYeoJohnsonShift"),
  function(object, ...) {

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
  signature(object = "transformationYeoJohnson"),
  function(object) {

    str <- paste0(
      "A ", ifelse(object@robust, "robust ", ""),
      ifelse(object@shift != 0.0, "shifted ", ""),
      "Yeo-Johnson transformation object"
    )

    if (object@complete) {
      cat(paste0(str, ".\n"))
      cat("  lambda: ", object@lambda, "\n")

      if (object@shift != 0.0) cat(paste0("  shift: ", object@shift, "\n"))

    } else {
      cat(paste0(str, " with unset transformation parameters.\n"))
    }
  }
)



# show (Yeo-Johnson (shift)) ---------------------------------------------------
setMethod(
  "show",
  signature(object = "transformationYeoJohnsonShift"),
  function(object) {

    str <- paste0(
      "A ", ifelse(object@robust, "robust ", ""),
      ifelse(object@shift != 0.0, "shifted ", ""),
      "Yeo-Johnson transformation object"
    )

    if (object@complete) {
      cat(paste0(str, ".\n"))
      cat("  lambda: ", object@lambda, "\n")
      cat("  shift: ", object@shift, "\n")

    } else {
      cat(paste0(str, " with unset transformation parameters.\n"))
    }
  }
)
