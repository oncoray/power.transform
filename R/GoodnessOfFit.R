#' @include TransformationObjects.R



#### .assess_transformation (generic) ------------------------------------------
setGeneric(
  ".assess_transformation",
  function(object, ...) standardGeneric(".assess_transformation"))


setMethod(
  ".assess_transformation",
  signature("transformationPowerTransform"),
  function(
    object,
    x,
    threshold,
    ...){

    # Sort x, if required.
    if(is.unsorted(x)) x <- sort(x)

    # Perform transformation.
    y <- .apply_transformation_parameters(
      object=object,
      x=x)

    # Compute expected z-score.
    z <- compute_expected_z(x=x)

    # Compute M-estimates for locality and scale
    robust_estimates <- huber_estimate(y)


  }
)
