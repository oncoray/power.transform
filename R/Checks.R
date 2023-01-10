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



.check_lambda <- function(x){

  # NULL is a valid value.
  if(is.null(x)) return(invisible(TRUE))

  # Otherwise, must be length 2, numeric, finite and sorted.
  if(!length(x) %in% c(1L, 2L)){
    stop(paste0("lambda should consist of 1 or 2 numeric values. ", length(x), " values were found."))
  }

  if(!is.numeric(x)){
    stop("lambda should consist of 1 or 2 numeric values. The values are not numeric.")
  }

  if(any(!is.finite(x))){
    stop("lambda should consist of 1 or 2 numeric values. One or both values are not finite.")
  }

  if(length(x) == 2){
    if(diff(x) == 0.0){
      stop("If lambda consists of 2 numeric values, these can not be identical.")
    }

    if(is.unsorted(x)){
      stop("If lambda consists of 2 numeric values, the values should be ordered by increasing value.")
    }
  }

  return(invisible(TRUE))
}
