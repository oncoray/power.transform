#' @include TransformationObjects.R
NULL


#' Assess normality of transformed data
#'
#' Not all data allows for a reasonable transformation to normality using power
#' transformation. For example, uniformly distributed data or multi-modal data
#' cannot be transformed to normality. This function computes a score expressing
#' deviation from normal distributions.
#'
#' @param transformer A transformer object created using
#'   `find_transformation_parameters`.
#' @param threshold A numeric value corresponding to the upper limit of the
#'   difference between the z-score of one instance of the transformed vector
#'   `x` and its expected z-score according to the normal distribution.
#' @param centre_width A numeric value between 0.0 and 1.0 that describes the
#'   width of the centre of the data. Only instances within this centre are
#'   compared against the threshold. Can be `NULL` to use all instances.
#' @param ... Unused arguments.
#'
#' @inheritParams power_transform
#'
#' @return Fraction of instances whose residual exceeds the threshold.
#' @export
#'
#' @examples
#' x <- exp(stats::rnorm(1000))
#' transformer <- find_transformation_parameters(
#'   x = x,
#'   method = "box_cox")
#'
#' assess_transformation(
#'   x = x,
#'   transformer = transformer)
assess_transformation <- function(
    x,
    transformer,
    threshold = 0.15,
    centre_width = 0.80,
    ...){

  # Prevent CRAN NOTE due to non-standard use of variables by data.table.
  z_expected <- NULL

  # Compute fit data.
  residual_data <- get_residuals(
    x = x,
    transformer = transformer)

  # Assign full weight to central elements. We do this to ensure that poor fits
  # don't get down-weighted in the centre, and in that way produce
  # log-likelihood values that are too optimistic.
  if(!is.null(centre_width)){
    if(centre_width < 0.0 || centre_width > 1.0){
      stop(paste0("centre_width should be NULL or between 0.0 and 1.0. Found: ", centre_width))
    }

    residual_data <- residual_data[abs(z_expected) <= stats::qnorm(0.50 + centre_width / 2), ]
    if(nrow(residual_data) == 0){
      stop(paste0("centre_width may be too small: no residuals were selected."))
    }
  }

  return(sum(abs(residual_data$residual) >= threshold) / nrow(residual_data))
}



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
    robust_estimates <- huber_estimate(y, tol=1E-3)

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


