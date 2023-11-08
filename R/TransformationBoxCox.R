#' @include TransformationObjects.R
NULL

# transformationBoxCox definition ----------------------------------------------

#' Box-Cox transformation object
#'
#' This class is used for Box-Cox transformations.
#'
#' @slot method Main transformation method, i.e. `"box_cox"`.
#' @slot robust Indicates whether a robust version of the Box-Cox transformation
#'   is used to set transformation parameters. The value depends on the `robust`
#'   argument of the `find_transformation_parameters` function.
#' @slot lambda Numeric lambda parameter for the Box-Cox transformation.
#' @slot shift Numeric shift parameter for the Box-Cox transformation. The value
#'   depends on the data used for setting transformation parameters. If all data
#'   are strictly positive, `shift` has a value of `0.0`. When negative or zero
#'   values are present, data are shifted to be strictly positive. If
#'   `invariant=TRUE` in the `find_transformation_parameters` function,
#'   `lambda`, `shift` and `scale` parameters are optimised simultaneously.
#' @slot scale Numeric scale parameter for the Box-Cox transformation. If
#'   `invariant=TRUE` in the `find_transformation_parameters` function,
#'   `lambda`, `shift` and `scale` parameters are optimised simultaneously.
#'   Otherwise, the `scale` parameter has a value of `1.0`.
#' @slot complete Indicates whether transformation parameters were set.
#'
#' @seealso [find_transformation_parameters]
#' @export
#' @rdname transformation_box_cox

setClass(
  "transformationBoxCox",
  contains = "transformationPowerTransform",
  slots = list(
    "method" = "character",
    "robust" = "logical",
    "lambda" = "numeric",
    "shift" = "numeric",
    "scale" = "numeric"
  ),
  prototype = list(
    "method" = "box_cox",
    "robust" = FALSE,
    "lambda" = 1.0,
    "shift" = 0.0,
    "scale" = 1.0
  )
)



# transformationBoxCoxInvariant definition -------------------------------------

#' @rdname transformation_box_cox
#' @export
setClass(
  "transformationBoxCoxInvariant",
  contains = "transformationBoxCox")



# .set_transformation_parameters (Box-Cox) -------------------------------------
setMethod(
  ".set_transformation_parameters",
  signature(object = "transformationBoxCox"),
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

    # Pre-scale and shift x. For conventional Box-Cox transformations, this does
    # nothing. For scale and shift-invariant Box-Cox transformations, the data
    # are standardised, and corresponding scale and shift parameters are
    # exported.
    init_transform_data <- ..standardise_data(
      object = object,
      x = x
    )
    x <- init_transform_data$x
    init_shift <- init_transform_data$shift
    init_scale <- init_transform_data$scale

    # Set minimum shift, in case any negative values are present.
    object <- ..set_minimum_shift(
      object = object,
      x = x
    )

    # Set lambda, in case a fixed lambda is provided.
    object <- ..set_lambda(
      object = object,
      lambda = lambda
    )

    # Get optimisation parameters to initialise and configure the optimiser.
    optimisation_parameters <- ..optimisation_parameters(
      object = object,
      x = x,
      lambda = lambda
    )

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

    if (!is.finite(optimised_parameters$lambda) && estimation_method != "mle") {
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

    # Update shift, scale and lambda values with the optimised parameters.
    if (!is.null(optimised_parameters$shift)) {
      if (is.finite(optimised_parameters$shift)) {
        object@shift <- optimised_parameters$shift * init_scale + init_shift

      } else if (!backup_use_default) {
        object@shift <- optimised_parameters$shift * init_scale + init_shift

      } else {
        object@shift <- 0.0
      }
    }

    if (!is.null(optimised_parameters$scale)) {
      if (is.finite(optimised_parameters$scale)) {
        object@scale <- optimised_parameters$scale * init_scale

      } else if (!backup_use_default) {
        object@scale <- optimised_parameters$scale * init_scale

      } else {
        object@scale <- 1.0
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



# .transform (Box-Cox) ---------------------------------------------------------
setMethod(
  ".transform",
  signature(object = "transformationBoxCox"),
  function(object, x, oob_action = "na", ...) {

    # Set up vector to copy new values to.
    y <- numeric(length(x))

    # Check oob_action
    .check_oob_action(oob_action)

    # Find non-positive and non-finite instances.
    na_entries <- which(x <= object@shift | !is.finite(x))

    if (length(na_entries) > 0 && oob_action == "na") {
      rlang::warn(
        message = paste0(
          "Box-cox power transforms are only defined for strictly positive values. ",
          length(na_entries), " zero, negative, NA or infinite values were found."),
        class = "power_transform_transform_invalid_values"
      )

      y[na_entries] <- NA_real_

    } else if (length(na_entries) > 0 && oob_action == "valid") {
      lower_bound <- ..get_value_bounds(object = object)[1]
      invalid_values <- is.finite(x) & x < lower_bound

      if (any(invalid_values)) {
        rlang::warn(
          message = paste0(
            "Box-cox power transforms are only defined for strictly positive values. ",
            length(na_entries), " zero or negative values were found."),
          class = "power_transform_transform_invalid_values"
        )

        y[invalid_values] <- lower_bound
      }
    }

    # Find remaining valid instances.
    valid_entries <- setdiff(seq_along(x), na_entries)

    # Perform transformation.
    if (length(valid_entries) > 0) {
      y[valid_entries] <- ..transform(
        object = object,
        x = x[valid_entries])
    }

    return(y)
  }
)



# ..transform (Box-Cox) --------------------------------------------------------
setMethod(
  "..transform",
  signature(object = "transformationBoxCox"),
  function(object, x, ...) {

    # Subtract shift and divide by scale.
    x <- (x - object@shift) / object@scale

    if (object@lambda == 0) {
      y <- log(x)

    } else {
      y <- (x^object@lambda - 1) / object@lambda
    }

    return(y)
  }
)


# .revert_transformation (Box-Cox) ---------------------------------------------
setMethod(
  ".revert_transform",
  signature(object = "transformationBoxCox"),
  function(object, x, ...) {

    # Revert transformation.
    if (object@lambda == 0) {
      y <- exp(x)

    } else {
      y <- (x * object@lambda + 1)^(1 / object@lambda)
    }

    # Apply scale and shift.
    y <- y * object@scale + object@shift

    return(y)
  })



# ..requires_shift_scale_optimisation (Box-Cox (invariant)) --------------------
setMethod(
  "..requires_shift_scale_optimisation",
  signature(object = "transformationBoxCoxInvariant"),
  function(object, ...) {
    return(TRUE)
  }
)



# ..get_default_shift_range (Box-Cox) ------------------------------------------
setMethod(
  "..get_default_shift_range",
  signature(object = "transformationBoxCox"),
  function(object, x, ...) {
    # Find the value that brings the entire distribution to 0. We need to
    # increment slightly to avoid x containing 0s.
    max_value <- min(x, na.rm = TRUE)
    max_value_offset <- 1.0
    min_value_offset <- 1.0

    # Find the typical, non-zero, distance between values. NA values should be
    # removed.
    dx <- unique(diff(sort(x, na.last = NA)))
    dx <- dx[dx > 0.0]

    # It shouldn't happen that all values are the same - but better check it.
    if (length(dx) > 0) max_value_offset <- stats::median(dx)

    # The increment should not grow too much.
    if (max_value_offset > 0.5) max_value_offset <- 0.5

    # Update min_value_offset.
    min_value_offset <- stats::median(x - max_value, na.rm = TRUE)
    if (min_value_offset < 1.0) min_value_offset <- 1.0

    # Set minimum shift value.
    min_value <- max_value - min_value_offset

    # Set maximum shift value.
    max_value <- max_value - max_value_offset

    # Set shift range. Occasionally, the values not be sorted.
    shift_range <- sort(c(min_value, max_value))

    return(shift_range)
  }
)



# ..get_default_scale_range (Box-Cox) ------------------------------------------
setMethod(
  "..get_default_scale_range",
  signature(object = "transformationBoxCox"),
  function(object, x, ...) {

    scale <- stats::IQR(x)
    if (scale == 0.0) scale <- stats::sd(x)
    if (!is.finite(scale)) scale <- 1.0
    if (scale == 0.0) scale <- 1.0

    # Limit scaling to half the scale up to twice the scale.
    scale_range <- c(scale / 2.0, scale * 2.0)

    return(scale_range)
  }
)



# ..get_default_lambda_range (Box-Cox) -----------------------------------------
setMethod(
  "..get_default_lambda_range",
  signature(object = "transformationBoxCox"),
  function(object, ...) {
    return(c(-100.0, 100.0))
  }
)



# ..standardise_data (Box-Cox (invariant)) -------------------------------------
setMethod(
  "..standardise_data",
  signature(object = "transformationBoxCoxInvariant"),
  function(object, x, ...) {

    shift <- stats::median(x)
    if (!is.finite(shift)) shift <- 0.0

    scale <- stats::IQR(x)
    if (scale == 0.0) scale <- stats::sd(x)
    if (!is.finite(scale)) scale <- 1.0
    if (scale == 0.0) scale <- 1.0

    # Shift and scale x.
    x <- (x - shift) / scale

    return(list(
      "x" = x,
      "shift" = shift,
      "scale" = scale
    ))
  }
)



# ..set_minimum_shift (Box-Cox) ------------------------------------------------
setMethod(
  "..set_minimum_shift",
  signature(object = "transformationBoxCox"),
  function(object, x, ...) {

    # Shift is necessary to avoid transformation with negative or zero values.
    if (all(x > 0.0)) return(object)

    # Warn to notify that zero or negative values are present. Avoid warning if
    # shifts are optimised.
    if (!is(object, "transformationBoxCoxInvariant")) {
      rlang::warn(
        message = paste0(
          "Box-cox power transforms are only defined for strictly positive values. ",
          "One or more zero or negative values were found in \"x\". ",
          "The values are shifted to induce strictly positive values."),
        class = "power_transform_transform_invalid_values"
      )
    }

    # Find shift range.
    x_range <- ..get_default_shift_range(
      object = object,
      x = x)

    # Set shift to its maximum value.
    object@shift <- max(x_range)

    return(object)
  }
)



# ..get_value_bounds (Box-Cox) -------------------------------------------------
setMethod(
  "..get_value_bounds",
  signature(object = "transformationBoxCox"),
  function(object, ...) {
    # By default, the minimum valid value is indicated by shift attribute. All
    # values that are equal to or smaller than the shift attribute, are invalid.
    minimum_shift <- object@shift

    # Increase the minimum shift by a small amount to create the valid lower
    # bound.
    return(c(minimum_shift + 1E-8, Inf))
  }
)



# ..set_lambda (Box-Cox) ------------------------------------------------------
setMethod(
  "..set_lambda",
  signature(object = "transformationBoxCox"),
  function(object, lambda, ...) {

    # Only update lambda if it is not a range.
    if (length(lambda) != 1) return(object)

    # Update lambda.
    object@lambda <- lambda

    return(object)
  }
)



# ..optimisation_parameters (Box-Cox) ------------------------------------------
setMethod(
  "..optimisation_parameters",
  signature(object = "transformationBoxCox"),
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



# ..optimisation_parameters (Box-Cox (invariant)) ------------------------------
setMethod(
  "..optimisation_parameters",
  signature(object = "transformationBoxCoxInvariant"),
  function(object, x, lambda, ...) {

    # Set up x-range.
    x_range <- ..get_default_shift_range(
      object = object,
      x = x
    )

    # Set up scale.
    s_range <- ..get_default_scale_range(
      object = object,
      x = x
    )

    if (length(lambda) == 1) {
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



# ..log_likelihood (Box-Cox) ---------------------------------------------------
setMethod(
  "..log_likelihood",
  signature(object = "transformationBoxCox"),
  function(
    object,
    x,
    w,
    sigma_hat_squared,
    ...) {

    # Compute the log likelihood under the assumption that the transformed
    # variable y follows the normal distribution. Note that shift is not
    # explicitly taken into account here, because this is handled by the
    # shifting x prior to computing the log-likelihood.
    return((object@lambda - 1.0) * sum(w * log(x - object@shift)) - sum(w) / 2.0 * log(sigma_hat_squared))
  }
)



# ..first_derivative (Box-Cox) -------------------------------------------------
setMethod(
  "..first_derivative",
  signature(object = "transformationBoxCox"),
  function(object, x) {
    return(x^(object@lambda - 1.0))
  }
)



# ..get_available_estimators (Box-Cox) -----------------------------------------
setMethod(
  "..get_available_estimators",
  signature(object = "transformationBoxCox"),
  function(object, ...) {

    available_estimators <- ..estimators_all()

    # Only allow Raymaekers and Rousseeuw's method for robust optimisation.
    if (!object@robust) {
      available_estimators <- setdiff(available_estimators, ..estimators_raymaekers_robust())
    }

    return(available_estimators)
  }
)


# ..get_available_estimators (Box-Cox (invariant)) -----------------------------
setMethod(
  "..get_available_estimators",
  signature(object = "transformationBoxCoxInvariant"),
  function(object, ...) {

    available_estimators <- ..estimators_all()

    # Raymaekers and Rousseeuw's method for robust optimisation is not suited
    # for simultaneous estimation of shift, scale and lambda parameters.
    available_estimators <- setdiff(available_estimators, ..estimators_raymaekers_robust())

    return(available_estimators)
  }
)



# show (Box-Cox) ---------------------------------------------------------------
setMethod(
  "show",
  signature("transformationBoxCox"),
  function(object) {

    str <- paste0(
      "A ", ifelse(object@robust, "robust ", ""),
      ifelse(object@shift != 0.0, "shifted ", ""),
      ifelse(object@scale != 1.0, "scaled ", ""),
      "Box-Cox transformation object"
    )

    if (object@complete) {
      cat(paste0(str, ".\n"))
      cat("  lambda: ", object@lambda, "\n")

      if (object@shift != 0.0) cat(paste0("  shift: ", object@shift, "\n"))
      if (object@scale != 1.0) cat(paste0("  scale: ", object@scale, "\n"))

    } else {
      cat(paste0(str, " with unset transformation parameters.\n"))
    }
  }
)
