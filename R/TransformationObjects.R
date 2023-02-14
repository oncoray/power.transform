# transformationPowerTransform  definition -------------------------------------

#' Generic transformation object
#'
#' This is the superclass for transformation objects.
#'
#' @slot method Main transformation method.
#' @slot complete Indicates whether transformation parameters were set.
#'
#' @export

setClass(
  "transformationPowerTransform",
  slots=list(
    "method" = "character",
    "complete" = "logical"),
  prototype=list(
    "method" = "none",
    "complete" = FALSE))

# transformationNone definition ------------------------------------------------

#' No transformation object
#'
#' This class is for transformers that do not alter the data.
#'
#' @slot method Main transformation method, i.e. `"none"`.
#' @slot complete Indicates whether transformation parameters were set.
#'
#' @export

setClass(
  "transformationNone",
  contains="transformationPowerTransform")



# .set_transformation_parameters (generic) -------------------------------------
setGeneric(
  ".set_transformation_parameters",
  function(object, ...) standardGeneric(".set_transformation_parameters"))



# .set_transformation_parameters (general) -------------------------------------
setMethod(
  ".set_transformation_parameters",
  signature(object = "transformationPowerTransform"),
  function(object, x, ...){
    # Default method.

    # Mark as complete.
    object@complete <- TRUE

    return(object)
  }
)



# .transform (generic) ---------------------------------------------------------
setGeneric(
  ".transform",
  function(object, ...) standardGeneric(".transform"))



# .transform (general) ---------------------------------------------------------
setMethod(
  ".transform",
  signature(object = "transformationPowerTransform"),
  function(object, x, ...){
    # Default method.

    # Set up vector to copy new values to.
    y <- numeric(length(x))

    # Find non-finite instances.
    na_entries <- which(!is.finite(x))
    if(length(na_entries) > 0){
      warning("One or more NA or infinite values were found.")

      y[na_entries] <- NA_real_
    }

    # Find remaining valid instances.
    valid_entries <- setdiff(seq_along(x), na_entries)

    # Perform transformation.
    if(length(valid_entries) > 0) y[valid_entries] <- ..transform(
      object = object,
      x = x[valid_entries])

    return(y)
  }
)



# ..transform (generic) --------------------------------------------------------
setGeneric(
  "..transform",
  function(object, ...) standardGeneric("..transform"))



# ..transform (general) --------------------------------------------------------
setMethod(
  "..transform",
  signature(object = "transformationPowerTransform"),
  function(object, x, ...){
    return(x)
  }
)



# .revert_transform (generic) --------------------------------------------------
setGeneric(
  ".revert_transform",
  function(object, ...) standardGeneric(".revert_transform"))



# .revert_transform (general) --------------------------------------------------
setMethod(
  ".revert_transform",
  signature(object = ".revert_transform"),
  function(object, x, ...){
    # Default method.

    return(x)
  }
)



# ..set_minimum_shift (generic) ------------------------------------------------
setGeneric(
  "..set_minimum_shift",
  function(object, ...) standardGeneric("..set_minimum_shift"))



# ..set_minimum_shift (general) -----------------------------------------------
setMethod(
  "..set_minimum_shift",
  signature(object = "transformationPowerTransform"),
  function(object, ...){
    # Default action is to not change the shift attribute.
    return(object)
  }
)



# ..set_lambda (generic) -------------------------------------------------------
setGeneric(
  "..set_lambda",
  function(object, ...) standardGeneric("..set_lambda"))



# ..set_lambda (general) ------------------------------------------------------
setMethod(
  "..set_lambda",
  signature(object = "transformationPowerTransform"),
  function(object, ...){
    # Default action is to not change the lambda attribute.
    return(object)
  }
)



# ..optimisation_parameters (generic) ------------------------------------------
setGeneric(
  "..optimisation_parameters",
  function(object, ...) standardGeneric("..optimisation_parameters"))



# ..optimisation_parameters (general) ------------------------------------------
setMethod(
  "..optimisation_parameters",
  signature(object = "transformationPowerTransform"),
  function(object, ...){
    # Default method

    return(NULL)
  }
)



# ..get_default_lambda_range (generic) -----------------------------------------
setGeneric(
  "..get_default_lambda_range",
  function(object, ...) standardGeneric("..get_default_lambda_range"))



# ..get_default_shift_range (generic) ------------------------------------------
setGeneric(
  "..get_default_shift_range",
  function(object, ...) standardGeneric("..get_default_shift_range"))




# ..log_likelihood (generic) ---------------------------------------------------
setGeneric(
  "..log_likelihood",
  function(object, ...) standardGeneric("..log_likelihood"))



# ..log_likelihood (general) ---------------------------------------------------
setMethod(
  "..log_likelihood",
  signature(object = "transformationPowerTransform"),
  function(object, ...){
    return(NA_real_)
  }
)



# ..first_derivative (generic) -------------------------------------------------
setGeneric(
  "..first_derivative",
  function(object, ...) standardGeneric("..first_derivative"))



# show (general) ---------------------------------------------------------------
setMethod(
  "show",
  signature("transformationPowerTransform"),
  function(object){
    cat(paste0("A generic power transformation object."))
  }
)



# ..get_available_estimators (generic) -----------------------------------------
setGeneric(
  "..get_available_estimators",
  function(object, ...) standardGeneric("..get_available_estimators"))
