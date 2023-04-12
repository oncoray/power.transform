#' @title Create transformation object skeleton
#'
#' @description Creates skeleton objects. This generates objects without fitting
#'   parameters. This is primarily intended for creating transformers
#'   externally, where fitting parameters are known.
#'
#' @param method Transformation method. Can be `"none"`, `"box_cox"` or
#'   `"yeo_johnson"`.
#' @param lambda Value of the transformation parameter lambda. Can also be
#'   changed using the `set_lambda` method.
#' @param shift Value of the shift parameter. Can also be changed using the
#'   `set_shift` method.
#'
#' @return A transformer object
#' @export
create_transformer_skeleton <- function(
    method,
    lambda = 1.0,
    shift = 0.0) {

  if (method == "none") {
    object <- methods::new("transformationNone")

  } else if (method == "box_cox") {
    object <- methods::new("transformationBoxCox")
    object <- set_lambda(object = object, lambda = lambda)
    object <- set_shift(object = object, shift = shift)

  } else if (method == "yeo_johnson") {
    object <- methods::new("transformationYeoJohnson")
    object <- set_lambda(object = object, lambda = lambda)
    object <- set_shift(object = object, shift = shift)

  } else {
    stop(paste0("The method was not recognised: ", method))
  }

  object <- .set_version(object)
  object@complete <- TRUE

  return(object)
}
