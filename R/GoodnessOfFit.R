#' @include TransformationObjects.R

#' Compute residuals of transformation to normality
#'
#' @inheritParams assess_transformation
#'
#' @return A `data.table` containing the expected (according to a normal
#'   distribution) and observed z-scores, and their difference as residuals.
#' @export
#'
#' @examples
#' x <- exp(stats::rnorm(1000))
#' transformer <- find_transformation_parameters(
#'   x = x,
#'   method = "box_cox")
#'
#' residual_data <- get_residuals(
#'   x = x,
#'   transformer = transformer)
get_residuals <- function(
    x,
    transformer,
    ...){

  # Perform checks on x.
  .check_data(x)

  # Check the transformer.
  .check_transformer(transformer)

  return(.get_fit_data(
    object = transformer,
    x = x))
}



#### .get_fit_data (generic) ---------------------------------------------------
setGeneric(
  ".get_fit_data",
  function(object, ...) standardGeneric(".get_fit_data"))



#### .get_fit_data (general) ---------------------------------------------------
setMethod(
  ".get_fit_data",
  signature("transformationPowerTransform"),
  function(
    object,
    x,
    ...){

    # Sort x, if required.
    if(is.unsorted(x)) x <- sort(x)

    # Perform transformation.
    y <- .apply_transformation_parameters(
      object=object,
      x=x)

    # Compute the expected z-score.
    z_expected <- compute_expected_z(x=x)

    # Compute M-estimates for locality and scale
    robust_estimates <- huber_estimate(y)

    # Compute the observed z-score.
    z_observed <- (y - robust_estimates$mu) / robust_estimates$sigma

    # Compute residuals.
    residual <- z_observed - z_expected

    return(data.table::data.table(
      "z_expected" = z_expected,
      "z_observed" = z_observed,
      "residual" = residual))
  }
)


