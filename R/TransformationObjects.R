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
    "method" = NA_character_,
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
    "method" = NA_character_,
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



#### .set_transformation_parameters (transformationPowerTransform) -------------
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



#### .set_transformation_parameters (transformationBoxCox) ---------------------
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
      min_value <- min(x)
      min_value_dx <- 1.0

      # Find the typical, non-zero, distance between values.
      dx <- unique(diff(sort(x)))
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
        x=data,
        type="box_cox")

    } else {
      # Standard method based on optimising log-likelihood of the normal
      # distribution.
      optimal_lambda <- suppressWarnings(
        stats::optimise(
          ..box_cox_loglik,
          interval=c(-10, 10),
          x=data,
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



#### .set_transformation_parameters (transformationYeoJohnson) -----------------
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
        x=data,
        type="yeo_johnson")

    } else {
      # Standard method based on optimising log-likelihood of the normal
      # distribution.
      optimal_lambda <- suppressWarnings(
        stats::optimise(
          ..yeo_johnson_loglik,
          interval=c(-10, 10),
          x=data,
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


#### .set_transformation_parameters (transformationYeoJohnson) -----------------
setMethod(
  ".set_transformation_parameters",
  signature("transformationYeoJohnsonShift"),
  function(object, x, ...){

    # Optimise lambda for Yeo-Johnson transformations.
    if(object@robust){
      # Robust method based on Raymaekers J, Rousseeuw PJ. Transforming
      # variables to central normality. Mach Learn. 2021.
      # doi:10.1007/s10994-021-05960-5
      lambda <- .transformation_robust_optimisation(
        x=data,
        type="yeo_johnson")

    } else {
      # Standard method based on optimising log-likelihood of the normal
      # distribution.
      optimal_lambda <- suppressWarnings(
        stats::optimise(
          ..yeo_johnson_loglik,
          interval=c(-10, 10),
          x=data,
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




#### .apply_transformation_parameters (generic) --------------------------------
setGeneric(
  ".apply_transformation_parameters",
  function(object, ...) standardGeneric(".apply_transformation_parameters"))



#### .apply_transformation_parameters (transformationPowerTransform) -----------
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
  signature(object="transformationYeoJohnson"),
  function(object, x, ...){

    # Apply shift.
    x <- x - object@shift

    # Perform transformation.
    x <- ..box_cox_transform(
      lambda=object@lambda,
      x=x)

    # Remove shift.
    x <- x + object@shift

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
      x=x)

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
      x=x)

    # Remove shift.
    x <- x + object@shift

    return(x)
  })
