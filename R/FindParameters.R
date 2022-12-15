#' @include TransformationObjects.R

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

    }

  } else if(transformation_method == "yeo_johnson"){

  } else {
    stop(
      paste0("Encountered an unknown transformation method: ",
             transformation_method)
    )
  }
}
