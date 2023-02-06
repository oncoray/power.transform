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
#' @export
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
#' @export
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
  function(object, x, lambda, ...){

    if(is.null(lambda)){
      # Pick very wide lambda-range, if lambda is NULL.
      lambda_range <- c(-100.0, 100.0)

    } else if(length(lambda) == 2){
      lambda_range <- lambda
    }

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
        x=x,
        type="box_cox",
        lambda_range=lambda_range)

      # Set lambda parameter.
      object@lambda <- lambda
      object@complete <- TRUE

    } else {
      # Standard method based on optimising log-likelihood of the normal
      # distribution.
      optimal_lambda <- suppressWarnings(
        stats::optimise(
          ..box_cox_loglik,
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



#### .set_transformation_parameters (Box-Cox (shift)) --------------------------
setMethod(
  ".set_transformation_parameters",
  signature("transformationBoxCoxShift"),
  function(object, x, lambda, optimiser="subplex", ...){

    if(is.null(lambda)){
      # Pick very wide lambda-range, if lambda is NULL.
      lambda_range <- c(-100.0, 100.0)

    } else if(length(lambda) == 2){
      lambda_range <- lambda

    } else {
      return(methods::callNextMethod())
    }

    # Set up initial search grid for shift and optimisation parameters to narrow
    # down the search area.
    search_grid <- box_cox_parameter_grid(
      x=x,
      lambda_range=lambda_range)

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

    # Check fall-back option.
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
        lower=c(search_grid$x_range[1], search_grid$lambda_range[1]),
        upper=c(search_grid$x_range[2], search_grid$lambda_range[2]),
        control=list("xtol_rel"=1e-3, ftol_rel=1e-4),
        x=x,
        type="box_cox",
        ...)

    } else if(optimiser == "subplex"){
      # SUBPLEX algorithm
      #
      # T. Rowan, “Functional Stability Analysis of Numerical
      # Algorithms”, Ph.D. thesis, Department of Computer Sciences, University
      # of Texas at Austin, 1990.

      results <- nloptr::sbplx(
        x0=c(search_grid$x_range[1], mean(search_grid$lambda_range)),
        fn=opt_fun,
        lower=c(search_grid$x_range[1], search_grid$lambda_range[1]),
        upper=c(search_grid$x_range[2], search_grid$lambda_range[2]),
        control=list("xtol_rel"=1e-3, ftol_rel=1e-4),
        x=x,
        type="box_cox",
        ...)

    } else if(optimiser == "nelder-mead"){
      # Nelder-Mead simplex algorithm
      #
      # J. A. Nelder and R. Mead, “A simplex method for function minimization,”
      # The Computer Journal 7, p. 308-313 (1965).
      #
      # M. J. Box, “A new method of constrained optimization and a comparison
      # with other methods,” Computer J. 8 (1), 42-52 (1965).

      results <- nloptr::neldermead(
        x0=c(search_grid$x_range[1], mean(search_grid$lambda_range)),
        fn=opt_fun,
        lower=c(search_grid$x_range[1], search_grid$lambda_range[1]),
        upper=c(search_grid$x_range[2], search_grid$lambda_range[2]),
        control=list("xtol_rel"=1e-3, ftol_rel=1e-4),
        x=x,
        type="box_cox",
        ...)

    } else if(optimiser == "optim-nelder-mead"){
      # Fall-back optimiser in case nloptr is not available. The Nelder-Mead
      # algorithm in stats::optim does not yield results as consistent as the
      # nloptr optimisers.
      results <- stats::optim(
        par=c(search_grid$x_range[1], mean(search_grid$lambda_range)),
        fn=opt_fun,
        gr=NULL,
        x=x,
        type="box_cox",
        ...,
        control=list(
          "abstol"=1E-5,
          "reltol"=1E-5))

    } else {
      stop(paste0("Optimiser not recognised: ", optimiser))
    }

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


#### .set_transformation_parameters (Yeo-Johnson (shift)) ----------------------
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

    # Set up initial search grid for shift and optimisation parameters to narrow
    # down the search area.
    search_grid <- yeo_johnson_parameter_grid(
      x = x,
      lambda_range = lambda_range)

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
        lower=c(search_grid$x_range[1], search_grid$lambda_range[1]),
        upper=c(search_grid$x_range[2], search_grid$lambda_range[2]),
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
          x0=c(search_grid$x_range[1], mean(search_grid$lambda_range)),
          fn=opt_fun,
          lower=c(search_grid$x_range[1], search_grid$lambda_range[1]),
          upper=c(search_grid$x_range[2], search_grid$lambda_range[2]),
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
        x0=c(search_grid$x_range[1], mean(search_grid$lambda_range)),
        fn=opt_fun,
        lower=c(search_grid$x_range[1], search_grid$lambda_range[1]),
        upper=c(search_grid$x_range[2], search_grid$lambda_range[2]),
        control=list("xtol_rel"=1e-3, ftol_rel=1e-4),
        x=x,
        ...,
        type="yeo_johnson")

    } else if(optimiser == "optim-nelder-mead"){
      # Fall-back optimiser in case nloptr is not available. The Nelder-Mead
      # algorithm in stats::optim does not yield results as consistent as the
      # nloptr optimisers.
      results <- stats::optim(
        par=c(search_grid$x_range[1], mean(search_grid$lambda_range)),
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
