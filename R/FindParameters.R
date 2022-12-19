#' @include TransformationObjects.R
NULL

#' Set transformation parameters
#'
#' @param x A numerical vector.
#' @param transformation_method One of the following methods for power
#'   transformation:
#'
#'   * `yeo_johnson`:
#'
#'   * `box_cox`:
#'
#'   * `none`: A fall-back method that will not transform values.
#'
#' @param robust
#' @param shift
#'
#' @return
#' @export
#'
#' @examples
find_transformation_parameters <- function(
    x,
    transformation_method="yeo_johnson",
    robust=TRUE,
    shift=TRUE){

  # Check transformation methods.
  if(!transformation_method %in% c("box_cox", "yeo_johnson", "none")){
    stop(paste0(
      "The transformation_method argument should be one of \"box_cox\", \"yeo_johnson\" or \"none\". ",
      "Found: ", transformation_method))
  }

  # Perform checks on x.
  if(is.factor(x) && transformation_method != "none"){
    warning("Power transformations are not applicable to categorical data.")
    transformation_method <- "none"
  }

  if(length(x) == 0){
    stop("x does not contain any values.")
  }

  if(!is.numeric(x)){
    stop("x does not contain numeric values.")
  }

  # Remove NA or inf values.
  x <- x[is.finite(x)]

  # Check length again.
  if(length(x) == 0){
    stop("x only contained NA or inf values.")
  }

  # Check number of unique values.
  n_unique_values <- length(unique(x))
  if(n_unique_values <= 3 && transformation_method != "none"){
    warning("x contains three or fewer unique values, and power transformation is not performed.")
    transformation_method <- "none"
  }

  if(n_unique_values > 3 && n_unique_values <= 10){
    warning("x contains ten or fewer unique values. Power transformation may be difficult.")
  }

  # Create transformation objects.
  if(transformation_method == "none"){
    object <- methods::new("transformationNone")

  } else if(transformation_method == "box_cox"){
    object <- methods::new(
      "transformationBoxCox",
      robust=robust)

    if(shift){
      object <- methods::new(
        "transformationBoxCoxShift",
        object)
    }

  } else if(transformation_method == "yeo_johnson"){
    object <- methods::new(
      "transformationYeoJohnson",
      robust=robust)

    if(shift){
      object <- methods::new(
        "transformationYeoJohnsonShift",
        object)
    }

  } else {
    stop(
      paste0("Encountered an unknown transformation method: ",
             transformation_method)
    )
  }

  # Set transformation parameters.
  object <- .set_transformation_parameters(object)

  return(object)
}



power_transform <- function(
    x,
    transformer=NULL,
    ...){

  # Create a transformer.
  if(is.null(transformer)){
    transformer <- do.call(
      find_transformation_parameters,
      c(list("x"=x),
        list(...)))
  }

  # Check that the transformer is a transformer.
  if(!is(transformer, "transformationPowerTransform")){
    stop(paste0(
      "The transformer object does not have the expected class. ",
      "Expected: transformationPowerTransform (or subclass). ",
      "Found: ", class(transformer)[1]))
  }

  # Check that transformer is complete.
  if(!transformer@complete) stop(paste0("Parameters for the transformer object were not fully set."))

  y <- .apply_transformation_parameters(
    object=transformer,
    x=x)

  return(y)
}



revert_power_transform <- function(
    y,
    transformer){

  if(missing(transformer)){
    stop(paste0(
      "A transformer object is required to revert the transformation."
    ))
  }

  # Check that the transformer is a transformer.
  if(!is(transformer, "transformationPowerTransform")){
    stop(paste0(
      "The transformer object does not have the expected class. ",
      "Expected: transformationPowerTransform (or subclass). ",
      "Found: ", class(transformer)[1]))
  }

  # Check that transformer is complete.
  if(!transformer@complete) stop(paste0("Parameters for the transformer object were not fully set."))

  # Revert transformation.
  x <- .invert_transformation(
    object=transformer,
    x=y)

  return(x)
}
