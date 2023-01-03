#' Generic transformation object
#'
#' This is the superclass for transformation objects.
#'
#' @slot method Main transformation method.
#' @slot complete Indicates whether transformation parameters were set.
#'
#' @export

setClass(
  "transformationPowerTransform",
  slots=list(
    "method" = "character",
    "complete" = "logical"),
  prototype=list(
    "method" = "none",
    "complete" = FALSE))



#' No transformation object
#'
#' This class is for transformers that do not alter the data.
#'
#' @slot method Main transformation method, i.e. `"none"`.
#' @slot complete Indicates whether transformation parameters were set.
#'
#' @export

setClass(
  "transformationNone",
  contains="transformationPowerTransform")



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
#'   `shift=TRUE` in the `find_transformation_parameters` function, `lambda` and
#'   `shift` parameters are optimised simultaneously.
#' @slot complete Indicates whether transformation parameters were set.
#'
#' @seealso [find_transformation_parameters]
#' @export
#' @rdname transformation_box_cox

setClass(
  "transformationBoxCox",
  contains="transformationPowerTransform",
  slots=list(
    "method" = "character",
    "robust" = "logical",
    "lambda" = "numeric",
    "shift" = "numeric"),
  prototype=list(
    "method" = "box_cox",
    "robust" = FALSE,
    "lambda" = NA_real_,
    "shift" = 0.0))

#' @rdname transformation_box_cox
setClass(
  "transformationBoxCoxShift",
  contains="transformationBoxCox")



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
setClass(
  "transformationYeoJohnsonShift",
  contains="transformationYeoJohnson")



#### .set_transformation_parameters (generic) ----------------------------------
setGeneric(
  ".set_transformation_parameters",
  function(object, ...) standardGeneric(".set_transformation_parameters"))



#### .set_transformation_parameters (general) ----------------------------------
setMethod(
  ".set_transformation_parameters",
  signature("transformationPowerTransform"),
  function(object, x, ...){
    # Default method.

    # Mark as complete.
    object@complete <- TRUE

    return(object)
  }
)



#### .set_transformation_parameters (Box-Cox) ----------------------------------
setMethod(
  ".set_transformation_parameters",
  signature("transformationBoxCox"),
  function(object, x, ...){

    if(any(x <= 0.0)){
      warning(paste0(
        "Box-cox power transforms are only defined for strictly positive values. ",
        "One or more zero or negative values are present in \"x\". ",
        "The values are shifted to induce strictly positive values."))

      # Find the (negative or zero) minimum value. We need to increment slightly
      # to avoid x containing 0s.
      min_value <- min(x, na.rm=TRUE)
      min_value_dx <- 1.0

      # Find the typical, non-zero, distance between values. NA values should be
      # removed.
      dx <- unique(diff(sort(x, na.last=NA)))
      dx <- dx[dx > 0.0]

      # It shouldn't happen that all values are the same - but better check it.
      if(length(dx) > 0){
        min_value_dx <- stats::median(dx)
      }

      # The increment should not grow too much.
      if(min_value_dx > 1.0) min_value_dx <- 1.0

      # Update shift value.
      object@shift <- min_value - min_value_dx

      # Shift data.
      x <- x - object@shift
    }

    # Optimise lambda for Box-Cox transformations.
    if(object@robust){
      # Robust method based on Raymaekers J, Rousseeuw PJ. Transforming
      # variables to central normality. Mach Learn. 2021.
      # doi:10.1007/s10994-021-05960-5
      lambda <- .transformation_robust_optimisation(
        x=x,
        type="box_cox")

    } else {
      # Standard method based on optimising log-likelihood of the normal
      # distribution.
      optimal_lambda <- suppressWarnings(
        stats::optimise(
          ..box_cox_loglik,
          interval=c(-4.0, 4.0),
          x=x,
          maximum=TRUE))

      lambda <- ifelse(
        is.finite(optimal_lambda$objective),
        optimal_lambda$maximum,
        1.0)
    }

    # Set lambda parameter.
    object@lambda <- lambda
    object@complete <- TRUE

    return(object)
  }
)



#### .set_transformation_parameters (Box-Cox (shift)) --------------------------
setMethod(
  ".set_transformation_parameters",
  signature("transformationBoxCoxShift"),
  function(object, x, ...){

    # Set up initial search grid for shift and optimisation parameters to narrow
    # down the search area.
    search_grid <- box_cox_parameter_grid(x)

    if(object@robust){
      # Sort values
      x <- sort(x)

      # Compute z-values according to the inverse cumulative density function.
      z_expected <- stats::qnorm(p=(seq_along(x) - 1/3) / (length(x) + 1/3))

      # Optimisation function.
      opt_fun <- .transformation_robust_shifted_optimisation

    } else {
      # No z required.
      z_expected <- NULL

      # Optimisation function.
      opt_fun <- .transform_shifted_optimisation
    }

    # Then compute the log-likelihood at each grid node.
    llf <- sapply(
      search_grid$parameter,
      opt_fun,
      x=x,
      z=z_expected,
      type="box_cox",
      shift_range=search_grid$x_range,
      lambda_range=search_grid$lambda_range)

    # Select the node with the highest llf.
    ii <- which.max(llf)

    # Select the direct neighbourhood of the llf (neighbouring lambda, x).
    initial_parameter <- search_grid$parameter[[ii]]

    shift_range <- select_neighbourhood(initial_parameter[1], search_grid$x)
    lambda_range <- select_neighbourhood(initial_parameter[2], search_grid$lambda)

    # Local approximation.
    results <- stats::optim(
      par=initial_parameter,
      fn=opt_fun,
      gr=NULL,
      x=x,
      z=z_expected,
      type="box_cox",
      shift_range=shift_range,
      lambda_range=lambda_range,
      control=list(
        "fnscale"=-1.0,
        "abstol"=1E-5,
        "reltol"=1E-5))

    # Extract optimal values.
    shift <- results$par[1]
    lambda <- results$par[2]

    if(!is.finite(results$value) || !is.finite(shift) || !is.finite(lambda)){
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



#### .set_transformation_parameters (Yeo-Johnson) ------------------------------
setMethod(
  ".set_transformation_parameters",
  signature("transformationYeoJohnson"),
  function(object, x, ...){

    # Optimise lambda for Yeo-Johnson transformations.
    if(object@robust){
      # Robust method based on Raymaekers J, Rousseeuw PJ. Transforming
      # variables to central normality. Mach Learn. 2021.
      # doi:10.1007/s10994-021-05960-5
      lambda <- .transformation_robust_optimisation(
        x=x,
        type="yeo_johnson")

    } else {
      # Standard method based on optimising log-likelihood of the normal
      # distribution.
      optimal_lambda <- suppressWarnings(
        stats::optimise(
          ..yeo_johnson_loglik,
          interval=c(-4.0, 4.0),
          x=x,
          maximum=TRUE))

      lambda <- ifelse(
        is.finite(optimal_lambda$objective),
        optimal_lambda$maximum,
        1.0)
    }

    # Set lambda parameter.
    object@lambda <- lambda
    object@complete <- TRUE

    return(object)
  }
)


#### .set_transformation_parameters (Yeo-Johnson (shift)) ----------------------
setMethod(
  ".set_transformation_parameters",
  signature("transformationYeoJohnsonShift"),
  function(object, x, ...){

    # Set up initial search grid for shift and optimisation parameters to narrow
    # down the search area.
    search_grid <- yeo_johnson_parameter_grid(x)

    if(object@robust){
      # Sort values
      x <- sort(x)

      # Compute z-values according to the inverse cumulative density function.
      z_expected <- stats::qnorm(p=(seq_along(x) - 1/3) / (length(x) + 1/3))

      # Optimisation function.
      opt_fun <- .transformation_robust_shifted_optimisation

    } else {
      # No z required.
      z_expected <- NULL

      # Optimisation function.
      opt_fun <- .transform_shifted_optimisation
    }

    # Then compute the log-likelihood at each grid node.
    llf <- sapply(
      search_grid$parameter,
      opt_fun,
      x=x,
      z=z_expected,
      type="yeo_johnson",
      shift_range=search_grid$x_range,
      lambda_range=search_grid$lambda_range)

    # Select the node with the highest llf.
    ii <- which.max(llf)

    # Select the direct neighbourhood of the llf (neighbouring lambda, x).
    initial_parameter <- search_grid$parameter[[ii]]

    shift_range <- select_neighbourhood(initial_parameter[1], search_grid$x)
    lambda_range <- select_neighbourhood(initial_parameter[2], search_grid$lambda)

    # Local approximation.
    results <- stats::optim(
      par=initial_parameter,
      fn=opt_fun,
      gr=NULL,
      x=x,
      z=z_expected,
      type="yeo_johnson",
      shift_range=shift_range,
      lambda_range=lambda_range,
      control=list(
        "fnscale"=-1.0,
        "abstol"=1E-5,
        "reltol"=1E-5))

    # Extract optimal values.
    shift <- results$par[1]
    lambda <- results$par[2]

    if(!is.finite(results$value) || !is.finite(shift) || !is.finite(lambda)){
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



#### .apply_transformation_parameters (generic) --------------------------------
setGeneric(
  ".apply_transformation_parameters",
  function(object, ...) standardGeneric(".apply_transformation_parameters"))



#### .apply_transformation_parameters (general) --------------------------------
setMethod(
  ".apply_transformation_parameters",
  signature("transformationPowerTransform"),
  function(object, x, ...){
    # Default method.

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
    if(length(valid_entries) > 0) y[valid_entries] <- x[valid_entries]

    return(y)
  }
)



#### .apply_transformation_parameters (Box-Cox) --------------------------------
setMethod(
  ".apply_transformation_parameters",
  signature(object="transformationBoxCox"),
  function(object, x, ...){

    # Apply shift.
    x <- x - object@shift

    # Set up vector to copy new values to.
    y <- numeric(length(x))

    # Find non-positive and non-finite instances.
    na_entries <- which(x <= 0.0 | !is.finite(x))
    if(length(na_entries) > 0){
      warning(paste0(
        "Box-cox power transforms are only defined for strictly positive values. ",
        "One or more zero, negative, NA or infinite values were found."))

      y[na_entries] <- NA_real_
    }

    # Find remaining valid instances.
    valid_entries <- setdiff(seq_along(x), na_entries)

    # Perform transformation.
    if(length(valid_entries) > 0){
      y[valid_entries] <- ..box_cox_transform(
        lambda=object@lambda,
        x=x[valid_entries],
        invert=FALSE)
    }

    return(y)
  })



#### .apply_transformation_parameters (Yeo-Johnson) ----------------------------
setMethod(
  ".apply_transformation_parameters",
  signature(object="transformationYeoJohnson"),
  function(object, x, ...){

    # Apply shift.
    x <- x - object@shift

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
    if(length(valid_entries) > 0){
      y[valid_entries] <- ..yeo_johnson_transform(
        lambda=object@lambda,
        x=x[valid_entries],
        invert=FALSE)
    }

    return(y)
  })



#### .invert_transformation (generic) ------------------------------------------
setGeneric(
  ".invert_transformation",
  function(object, ...) standardGeneric(".invert_transformation"))



#### .invert_transformation (general) ------------------------------------------
setMethod(
  ".invert_transformation",
  signature("transformationPowerTransform"),
  function(object, x, ...){
    # Default method.

    return(x)
  }
)



#### .invert_transformation (Box-Cox) ------------------------------------------
setMethod(
  ".invert_transformation",
  signature(object="transformationBoxCox"),
  function(object, x, ...){

    # Perform transformation.
    x <- ..box_cox_transform(
      lambda=object@lambda,
      x=x,
      invert=TRUE)

    # Apply shift.
    x <- x + object@shift

    return(x)
  })



#### .invert_transformation (Yeo-Johnson) --------------------------------------
setMethod(
  ".invert_transformation",
  signature(object="transformationYeoJohnson"),
  function(object, x, ...){

    # Perform transformation.
    x <- ..yeo_johnson_transform(
      lambda=object@lambda,
      x=x,
      invert=TRUE)

    # Apply shift.
    x <- x + object@shift

    return(x)
  })



#### show (general) ------------------------------------------------------------
setMethod(
  "show",
  signature("transformationPowerTransform"),
  function(object){
    cat(paste0("A generic power transformation object."))
  }
)



#### show (Box-Cox) ------------------------------------------------------------
setMethod(
  "show",
  signature("transformationBoxCox"),
  function(object){

    str <- paste0(
      "A ", ifelse(object@robust, "robust ", ""),
      ifelse(object@shift != 0.0, "shifted ", ""),
      "Box-Cox transformation object"
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



#### show (Yeo-Johnson) --------------------------------------------------------
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



#### show (Yeo-Johnson (shift)) ------------------------------------------------
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
