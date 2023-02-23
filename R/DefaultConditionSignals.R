.warn_poor_weighting_method <- function(x){

  # Throw a warning.
  rlang::warn(
    message = paste0(
      "The ", x, " weighting function leads to worse overall parameter estimations than non-robust transformations ",
      "and should be avoided. Use one of the weighting functions based on empirical probability instead."),
    .frequency = "once",
    .frequency_id = "poor_weight_method"
  )

  return(invisible(TRUE))
}
