# transformationPowerTransform  definition -------------------------------------

#' Generic transformation object
#'
#' This is the superclass for transformation objects.
#'
#' @slot method Main transformation method.
#' @slot complete Indicates whether transformation parameters were set.
#' @slot version Version of the power.transform package that was used to create
#'   the transformation objecst.
#'
#' @export

setClass(
  "transformationPowerTransform",
  slots = list(
    "method" = "character",
    "complete" = "logical",
    "version" = "ANY"),
  prototype = list(
    "method" = "none",
    "complete" = FALSE,
    "version" = NULL))


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
  contains = "transformationPowerTransform")




# .set_transformation_parameters (generic) -------------------------------------
setGeneric(
  ".set_transformation_parameters",
  function(object, ...) standardGeneric(".set_transformation_parameters"))



# .set_transformation_parameters (general) -------------------------------------
setMethod(
  ".set_transformation_parameters",
  signature(object = "transformationPowerTransform"),
  function(object, x, ...) {
    # Default method.

    # Add package version.
    object <- .set_version(object = object)

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
  function(object, x, ...) {
    # Default method.

    # Set up vector to copy new values to.
    y <- numeric(length(x))

    # Find non-finite instances.
    na_entries <- which(!is.finite(x))
    if (length(na_entries) > 0) {
      rlang::warn(
        message = paste0(
          "Power transformations are only defined for finite values. ",
          length(na_entries), " NA or infinite values were found."),
        class = "power_transform_transform_invalid_values"
      )

      y[na_entries] <- NA_real_
    }

    # Find remaining valid instances.
    valid_entries <- setdiff(seq_along(x), na_entries)

    # Perform transformation.
    if (length(valid_entries) > 0) y[valid_entries] <- ..transform(
      object = object,
      x = x[valid_entries])

    return(y)
  }
)



# .transform (none) ------------------------------------------------------------
setMethod(
  ".transform",
  signature(object = "transformationNone"),
  function(object, x, ...) {
    # In case no transformation takes place, we can just return x.
    return(x)
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
  function(object, x, ...) {
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
  signature(object = "transformationPowerTransform"),
  function(object, x, ...) {
    # Default method.

    return(x)
  }
)



# ..requires_shift_scale_optimisation (generic) --------------------------------
setGeneric(
  "..requires_shift_scale_optimisation",
  function(object, ...) standardGeneric("..requires_shift_scale_optimisation"))




# ..requires_shift_scale_optimisation (general) --------------------------------
setMethod(
  "..requires_shift_scale_optimisation",
  signature(object = "transformationPowerTransform"),
  function(object, ...) {
    return(FALSE)
  }
)



# ..standardise_data (generic) -------------------------------------------------
setGeneric(
  "..standardise_data",
  function(object, ...) standardGeneric("..standardise_data")
)



# ..standardise_data (general) -------------------------------------------------
setMethod(
  "..standardise_data",
  signature(object = "transformationPowerTransform"),
  function(object, x, ...) {
    return(list(
      "x" = x,
      "shift" = 0.0,
      "scale" = 1.0
    ))
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
  function(object, ...) {
    # Default action is to not change the shift attribute.
    return(object)
  }
)



# ..get_value_bounds (generic) -------------------------------------------------
setGeneric(
  "..get_value_bounds",
  function(object, ...) standardGeneric("..get_value_bounds"))



# ..get_value_bounds (general) -------------------------------------------------
setMethod(
  "..get_value_bounds",
  signature(object = "transformationPowerTransform"),
  function(object, ...) {
    # By default, there is no minimum shift.
    return(c(-Inf, Inf))
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
  function(object, ...) {
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
  function(object, ...) {
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



# ..get_default_scale_range (generic) ------------------------------------------
setGeneric(
  "..get_default_scale_range",
  function(object, ...) standardGeneric("..get_default_scale_range")
)



# ..log_likelihood (generic) ---------------------------------------------------
setGeneric(
  "..log_likelihood",
  function(object, ...) standardGeneric("..log_likelihood"))



# ..log_likelihood (general) ---------------------------------------------------
setMethod(
  "..log_likelihood",
  signature(object = "transformationPowerTransform"),
  function(object, ...) {
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
  signature(object = "transformationPowerTransform"),
  function(object) {
    cat(paste0("A generic power transformation object."))
  }
)



# ..get_available_estimators (generic) -----------------------------------------
setGeneric(
  "..get_available_estimators",
  function(object, ...) standardGeneric("..get_available_estimators"))



# .set_version (generic) -------------------------------------------------------
setGeneric(
  ".set_version",
  function(object, ...) standardGeneric(".set_version")
)



# .set_version (general) -------------------------------------------------------
setMethod(
  ".set_version",
  signature(object = "transformationPowerTransform"),
  function(object, ...) {
    # Set package version.
    object@version <- utils::packageVersion("power.transform")

    return(object)
  }
)
