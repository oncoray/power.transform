parameter_list <- list()
ii <- 1
for(method in c("box_cox", "yeo_johnson", "none")){
  for(robust in c(FALSE, TRUE)){
    for(shift in c(FALSE, TRUE)){
      parameter_list[[ii]] <- list(
        "method"=method,
        "robust"=robust,
        "shift"=shift
      )

      ii <- ii + 1
    }
  }
}

# Set seed.
set.seed(19L)

# Draw 10000 samples according to a normal distribution.
x_normal <- stats::rnorm(10000)

# Make exponentially distributed.
x_exponential <- exp(x_normal)

# Draw bimodal with large separation
x_bimodal_large <- c(
  stats::rnorm(5000, mean=0.0),
  stats::rnorm(5000, mean=6.0))

# Iterate over all parameter sets.
for(ii in seq_along(parameter_list)){
  # Get parameter set.
  parameter_set <- parameter_list[[ii]]

  #### Normally distributed data -----------------------------------------------

  # Create the transformer.
  transformer <- suppressWarnings(do.call(
    power.transform::find_transformation_parameters,
    args=c(list("x"=x_normal),
           parameter_set)))

  # Compute the fraction of instances which exceed the threshold.
  poor_fit_value <- power.transform::assess_transformation(
    x = x_exponential,
    transformer = transformer)

  #### Exponentially distributed data ------------------------------------------

  # Create the transformer.
  transformer <- suppressWarnings(do.call(
    power.transform::find_transformation_parameters,
    args=c(list("x"=x_exponential),
           parameter_set)))

  # Compute the fraction of instances which exceed the threshold.
  poor_fit_value <- power.transform::assess_transformation(
    x = x_exponential,
    transformer = transformer)

  #### Bi-modal distribution (large separation) --------------------------------

  # Create the transformer.
  transformer <- suppressWarnings(do.call(
    power.transform::find_transformation_parameters,
    args=c(list("x"=x_bimodal_large),
           parameter_set)))

  # Compute the fraction of instances which exceed the threshold.
  poor_fit_value <- power.transform::assess_transformation(
    x = x_bimodal_large,
    transformer = transformer)

}
