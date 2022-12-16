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

    # Optimise lambda for Box-Cox transformations.
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












.transformation_robust_optimisation <- function(x, type){
  # This follows the algorithm from Raymaekers J, Rousseeuw PJ. Transforming
  # variables to central normality. Mach Learn. 2021.
  # doi:10.1007/s10994-021-05960-5

  # Sort x.
  x <- x[order(x)]

  # Compute z-values according to the inverse cumulative density function.
  z_expected <- stats::qnorm(p=(seq_along(x) - 1/3) / (length(x) + 1/3))

  # Step 1: Compute initial estimate for lambda.
  # Standard method based on optimising log-likelihood of the normal
  # distribution.

  optimal_lambda <- suppressWarnings(
    stats::optimise(
      ..transformation_rectified_optimisation,
      interval=c(-10, 10),
      x=x,
      z=z_expected,
      type=type,
      maximum=FALSE))

  optimal_lambda <- ifelse(
    is.finite(optimal_lambda$objective),
    optimal_lambda$minimum,
    1.0)

  # Step 2: Compute lambda from reweighted maximum likelihood.
  optimal_lambda <- ..transformation_reweighted_optimisation(
    lambda_0=optimal_lambda,
    x=x,
    type=type,
    ii=1L)

  # Step 3: Compute lambda from reweighted maximum likelihood again.
  optimal_lambda <- ..transformation_reweighted_optimisation(
    lambda_0=optimal_lambda,
    x=x,
    type=type,
    ii=2L)

  if(!is.finite(optimal_lambda)) return(1.0)

  return(optimal_lambda)
}




..transformation_rectified_optimisation <- function(lambda, x, z, type){

  rectifier_FUN <- switch(
    type,
    "box_cox"=..box_cox_transform_rectified,
    "yeo_johnson"=..yeo_johnson_transform_rectified
  )

  # Compute values after power transformation with the appropriate lambda value.
  y <- do.call(
    rectifier_FUN,
    args=list(
      "lambda"=lambda,
      "x"=x))

  # Compute M-estimates for locality and scale
  robust_estimates <- huber_estimate(y)

  # Check problematic values.
  if(!is.finite(robust_estimates$sigma)) return(NA_real_)
  if(robust_estimates$sigma == 0.0) return(NA_real_)

  # Compute residuals.
  residual <- (y - robust_estimates$mu) / robust_estimates$sigma - z

  # Compute Tukey bisquare function to truncate weights of outlier residuals.
  truncated_weights <- numeric(length(residual)) + 1.0
  valid_residuals <- which(abs(residual) <= 0.5)

  if(length(valid_residuals) > 0){
    truncated_weights[valid_residuals] <- 1.0 - (1.0 - (residual[valid_residuals] / 0.5)^2)^3
  }

  return(sum(truncated_weights))
}



..transformation_reweighted_optimisation <- function(lambda_0, x, type, ii){

  if(!is.finite(lambda_0)) return(NA_real_)

  # Set transformation function.
  if(ii == 1){
    # For the initial step, use the rectified transformations
    transform_FUN <- switch(
      type,
      "box_cox"=..box_cox_transform_rectified,
      "yeo_johnson"=..yeo_johnson_transform_rectified
    )

  } else {
    transform_FUN <- switch(
      type,
      "box_cox"=..box_cox_transform,
      "yeo_johnson"=..yeo_johnson_transform
    )
  }

  # Set log-likelihood function.
  loglik_FUN <- switch(
    type,
    "box_cox"=..box_cox_loglik,
    "yeo_johnson"=..yeo_johnson_loglik
  )

  # Perform transformation for lambda_0
  y <- do.call(
    transform_FUN,
    args=list(
      "lambda"=lambda_0,
      "x"=x
    )
  )

  # Compute M-estimates for locality and scale
  robust_estimates <- huber_estimate(y)

  # Check problematic values.
  if(!is.finite(robust_estimates$sigma)) return(NA_real_)
  if(robust_estimates$sigma == 0.0) return(NA_real_)

  # Compute weights.
  weights <- as.numeric(abs(y - robust_estimates$mu) / robust_estimates$sigma <= stats::qnorm(0.99))
  if(sum(weights) == 0.0) return(NA_real_)

  # Compute optimal lambda.
  optimal_lambda <- suppressWarnings(
    stats::optimise(
      loglik_FUN,
      interval=c(-10, 10),
      x=x,
      w=weights,
      maximum=TRUE))

  lambda <- ifelse(
    is.finite(optimal_lambda$objective),
    optimal_lambda$maximum,
    1.0)

  return(lambda)
}
