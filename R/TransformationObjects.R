setClass(
  "transformationPowerTransform",
  slots=list(
    "method" = "character",
    "complete" = "logical"),
  prototype=list(
    "method" = "none",
    "complete" = FALSE))

setClass(
  "transformationNone",
  contains="transformationPowerTransform")

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

setClass(
  "transformationBoxCoxShift",
  contains="transformationBoxCox")

setClass(
  "transformationYeoJohnson",
  contains="transformationPowerTransform",
  slots=list(
    "method" = "character",
    "robust" = "logical",
    "lambda" = "numeric"),
  prototype=list(
    "method" = "yeo_johnson",
    "robust" = FALSE,
    "lambda" = NA_real_))

setClass(
  "transformationYeoJohnsonShift",
  contains="transformationYeoJohnson",
  slots=list(
    "shift" = "numeric"),
  prototype=list(
    "shift" = NA_real_))



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
          interval=c(-10, 10),
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

    # Get range for shift parameters.
    shift_range <- box_cox_shift_range(x)

    # Optimise shift and lambda for Box-Cox transformations. Choose mean
    # shift and no transformation as initial values.
    if(object@robust){

      # Robust transformation.
      results <- stats::optim(
        par=c(mean(shift_range), 1.0),
        fn=.transformation_robust_shifted_optimisation,
        gr=NULL,
        x=x,
        type="box_cox",
        shift_range=shift_range,
        lambda_range=c(-6, 4),
        control=list(
          "fnscale"=-1.0,
          "abstol"=1E-5,
          "reltol"=1E-5))

    } else {

      # Conventional transformation.
      results <- stats::optim(
        par=c(mean(shift_range), 1.0),
        fn=.transform_shifted_optimisation,
        gr=NULL,
        x=x,
        type="box_cox",
        shift_range=shift_range,
        lambda_range=c(-6, 4),
        control=list(
          "fnscale"=-1.0,
          "abstol"=1E-5,
          "reltol"=1E-5))
    }

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
          interval=c(-10, 10),
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

    # Get range for shift parameters.
    shift_range <- yeo_johnson_shift_range(x)

    # Optimise shift and lambda for Yeo-Johnson transformations. Choose mean
    # shift and no transformation as initial values.
    if(object@robust){

      # Robust transformation.
      results <- stats::optim(
        par=c(mean(shift_range), 1.0),
        fn=.transformation_robust_shifted_optimisation,
        gr=NULL,
        x=x,
        type="yeo_johnson",
        shift_range=shift_range,
        lambda_range=c(-6, 4),
        control=list(
          "fnscale"=-1.0,
          "abstol"=1E-5,
          "reltol"=1E-5))

    } else {

      # Conventional transformation.
      results <- stats::optim(
        par=c(mean(shift_range), 1.0),
        fn=.transform_shifted_optimisation,
        gr=NULL,
        x=x,
        type="yeo_johnson",
        shift_range=shift_range,
        lambda_range=c(-6, 4),
        control=list(
          "fnscale"=-1.0,
          "abstol"=1E-5,
          "reltol"=1E-5))
    }

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

    return(x)
  }
)



#### .apply_transformation_parameters (Box-Cox) --------------------------------
setMethod(
  ".apply_transformation_parameters",
  signature(object="transformationBoxCox"),
  function(object, x, ...){

    # Apply shift.
    x <- x - object@shift

    # Perform transformation.
    x <- ..box_cox_transform(
      lambda=object@lambda,
      x=x,
      invert=FALSE)

    return(x)
  })



#### .apply_transformation_parameters (Yeo-Johnson) ----------------------------
setMethod(
  ".apply_transformation_parameters",
  signature(object="transformationYeoJohnson"),
  function(object, x, ...){

    # Perform transformation.
    x <- ..yeo_johnson_transform(
      lambda=object@lambda,
      x=x,
      invert=FALSE)

    return(x)
  })



#### .apply_transformation_parameters (Yeo-Johnson (shift)) --------------------
setMethod(
  ".apply_transformation_parameters",
  signature(object="transformationYeoJohnsonShift"),
  function(object, x, ...){

    # Apply shift.
    x <- x - object@shift

    # Perform transformation.
    x <- ..yeo_johnson_transform(
      lambda=object@lambda,
      x=x,
      invert=FALSE)

    return(x)
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

    return(x)
  })



#### .invert_transformation (Yeo-Johnson (shift)) ------------------------------
setMethod(
  ".invert_transformation",
  signature(object="transformationYeoJohnsonShift"),
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
      "Yeo-Johnson transformation object"
    )

    if(object@complete){
      cat(paste0(str, ".\n"))
      cat("  lambda: ", object@lambda, "\n")

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
