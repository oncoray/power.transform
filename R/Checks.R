.check_data <- function(x){
  # Perform checks on x.
  if(length(x) == 0){
    stop("x does not contain any values.")
  }

  if(is.factor(x)){
    stop("x is categorical, and power transformations are not applicable.")
  }

  if(!is.numeric(x)){
    stop("x does not contain numeric values.")
  }

  if(all(!is.finite(x))){
    stop("x only contains NA or inf values.")
  }

  return(invisible(TRUE))
}



.check_transformer <- function(x){

  if(!is(x, "transformationPowerTransform")){
    stop(paste0(
      "The transformer object does not have the expected class. ",
      "Expected: transformationPowerTransform (or subclass). ",
      "Found: ", class(x)[1]))
  }

  if(!x@complete){
    stop(paste0(
      "The transformer object did not have all fitting parameters set."))
  }

  return(invisible(TRUE))
}
