
#' Set transformation parameters
#'
#' `find_transformation_parameters` is used to find optimal parameters for
#' univariate transformation to normality.
#'
#' @param x A vector with numeric values.
#' @param method One of the following methods for power transformation:
#'
#'   * `box_cox`: Transformation using the Box-Cox transformation (Box and Cox,
#'   1964). The Box-Cox transformation requires that all data are strictly
#'   positive. Features that contain zero or negative values cannot be
#'   transformed using this transformation. In their work, Box and Cox define a
#'   shifted variant. We use this variant to shift values to a strictly positive
#'   range, when negative values are present. The Box-Cox transformation relies
#'   on a single parameter lambda, which is estimated through maximisation of
#'   the log-likelihood function corresponding to a normal distribution.
#'
#'   * `yeo_johnson`:Transformation using the Yeo-Johnson
#'   transformation (Yeo and Johnson, 2000). Unlike the Box-Cox transformation,
#'   the Yeo-Johnson transformation allows for negative and positive values.
#'   Like the Box-Cox transformation, this transformation relies on a single
#'   parameter lambda, which is estimated through maximisation of the
#'   log-likelihood function corresponding to a normal distribution.
#'
#'   * `none`: A fall-back method that will not transform values.
#'
#' @param robust Flag for using a robust version of Box-Cox or Yeo-Johnson
#'   transformation, as defined by Raymaekers and Rousseeuw (2021). This version
#'   is less sensitive in the presence outliers.
#' @param shift Flag for using a version of Box-Cox or Yeo-Johnson
#'   transformation that simultaneously optimises location in addition to the
#'   lambda parameter.
#' @param lambda Single lambda value, or range of lambda values that should be
#'   considered. Default: c(4.0, 6.0). Can be `NULL` to force optimisation
#'   without a constraint in lambda values.
#' @param ... Unused parameters.
#'
#' @return A transformer object that can be used to transform values.
#' @export
#'
#' @references 1. Yeo, I. & Johnson, R. A. A new family of power transformations
#'   to improve normality or symmetry. Biometrika 87, 954–959 (2000).
#'
#'   1. Box, G. E. P. & Cox, D. R. An analysis of transformations. J. R. Stat.
#'   Soc. Series B Stat. Methodol. 26, 211–252 (1964).
#'
#'   1. Raymaekers, J., Rousseeuw,  P. J. Transforming variables to central
#'   normality. Mach Learn. (2021).
#'
#' @examples
#' x <- exp(stats::rnorm(1000))
#' transformer <- find_transformation_parameters(
#'   x = x,
#'   method = "box_cox")
find_transformation_parameters <- function(
    x,
    method = "yeo_johnson",
    robust = TRUE,
    shift = TRUE,
    lambda = c(-4.0, 6.0),
    ...){

  # Check transformation methods.
  if(!method %in% c("box_cox", "yeo_johnson", "none")){
    stop(paste0(
      "The method argument should be one of \"box_cox\", \"yeo_johnson\" or \"none\". ",
      "Found: ", method))
  }

  # Perform checks on x.
  .check_data(x)

  # Remove NA or inf values.
  x <- x[is.finite(x)]

  # Check number of unique values.
  n_unique_values <- length(unique(x))
  if(n_unique_values <= 3 && method != "none"){
    warning("x contains three or fewer unique values, and power transformation is not performed.")
    method <- "none"
  }

  if(n_unique_values > 3 && n_unique_values <= 10){
    warning("x contains ten or fewer unique values. Power transformation may be difficult.")
  }

  # Create transformation objects.
  if(method == "none"){
    object <- methods::new("transformationNone")

  } else if(method == "box_cox"){
    object <- methods::new(
      "transformationBoxCox",
      robust=robust)

    if(shift){
      object <- methods::new(
        "transformationBoxCoxShift",
        object)
    }

    # Check lambda.
    .check_lambda(lambda)

  } else if(method == "yeo_johnson"){
    object <- methods::new(
      "transformationYeoJohnson",
      robust=robust)

    if(shift){
      object <- methods::new(
        "transformationYeoJohnsonShift",
        object)
    }

    # Check lambda.
    .check_lambda(lambda)

  } else {
    stop(
      paste0("Encountered an unknown transformation method: ",
             method)
    )
  }

  # Sort x.
  x <- sort(x)

  # Set transformation parameters.
  object <- .set_transformation_parameters(
    object = object,
    x = x,
    lambda = lambda,
    ...)

  return(object)
}



#' Transform values
#'
#' `power_transform` transforms numeric values to normality.
#'
#' @param x A vector with numeric values that should be transformed to
#'   normality.
#' @param transformer A transformer object created using
#'   `find_transformation_parameters`. If `NULL`, a transformer is generated
#'   internally.
#' @param ... Parameters passed on to `find_transformation_parameters`.
#'
#' @inheritDotParams find_transformation_parameters
#'
#' @return A vector of transformed values of `x`.
#' @export
#'
#' @seealso [find_transformation_parameters]
#' @examples
#' x <- exp(stats::rnorm(1000))
#' y <- power_transform(
#'   x = x,
#'   method = "box_cox")
power_transform <- function(
    x,
    transformer = NULL,
    ...){

  # Create a transformer.
  if(is.null(transformer)){
    transformer <- do.call(
      find_transformation_parameters,
      c(list("x"=x), list(...)))
  }

  # Check the transformer.
  .check_transformer(transformer)

  # Check that x is numeric.
  if(!is.numeric(x)){
    stop("x does not contain numeric values.")
  }

  # Transform data using the transformer.
  y <- .transform(
    object=transformer,
    x=x)

  return(y)
}



#' Revert transformation
#'
#' `revert_power_transform` reverts the transformation of numeric values to
#' normality.
#'
#' @param y A vector with numeric values that was previously transformed to
#'   normality.
#' @param transformer A transformer object created using
#'   `find_transformation_parameters` that was used to transform the values to
#'   normality previously. Cannot be `NULL`.
#'
#' @return A vector of values.
#' @export
#'
#' @examples
#' x0 <- exp(stats::rnorm(1000))
#'
#' transformer <- find_transformation_parameters(
#'   x = x0,
#'   method = "box_cox")
#'
#' y <- power_transform(
#'   x = x0,
#'   transformer = transformer)
#'
#' x1 <- revert_power_transform(
#'   y = y,
#'   transformer = transformer)
revert_power_transform <- function(
    y,
    transformer){

  if(missing(transformer)){
    stop(paste0(
      "A transformer object is required to revert the transformation."
    ))
  }

  # Check the transformer.
  .check_transformer(transformer)

  # Revert transformation.
  x <- .revert_transform(
    object=transformer,
    x=y)

  return(x)
}
