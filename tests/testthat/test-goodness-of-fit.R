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

# Draw bimodal with intermediate separation
x_bimodal_intermediate <- c(
  stats::rnorm(5000, mean=0.0),
  stats::rnorm(5000, mean=3.0))

# Draw bimodal with small separation
x_bimodal_small <- c(
  stats::rnorm(5000, mean=0.0),
  stats::rnorm(5000, mean=1.0))

# Draw uniform.
x_uniform <- stats::runif(10000, min=1E-5)

# Draw Gamma (left skewed)
x_gamma_left <- stats::rgamma(10000, shape=2.0, scale=2.0)

# Draw Gamma (monotonic decreasing)
x_gamma_monotonic <- stats::rgamma(10000, shape=1.0, scale=2.0)

# Iterate over all parameter sets.
for(ii in seq_along(parameter_list)){

  #### Normally distributed data -----------------------------------------------

  # Create the transformer.
  transformer <- suppressWarnings(do.call(
    power.transform::find_transformation_parameters,
    args=c(list("x"=x_normal),
           parameter_list[[ii]])))

  # Compute the fraction of instances which exceed the threshold.
  power.transform::assess_transformation(
    x = x_exponential,
    transformer = transformer)

  #### Exponentially distributed data ------------------------------------------

  # Create the transformer.
  transformer <- suppressWarnings(do.call(
    power.transform::find_transformation_parameters,
    args=c(list("x"=x_exponential),
           parameter_list[[ii]])))

  # Compute the fraction of instances which exceed the threshold.
  power.transform::assess_transformation(
    x = x_exponential,
    transformer = transformer)

  #### Left-skewed gamma distribution ------------------------------------------

  # Create the transformer.
  transformer <- suppressWarnings(do.call(
    power.transform::find_transformation_parameters,
    args=c(list("x"=x_gamma_left),
           parameter_list[[ii]])))

  # Compute the fraction of instances which exceed the threshold.
  power.transform::assess_transformation(
    x = x_gamma_left,
    transformer = transformer)

  #### Monotonically decreasing gamma distribution -----------------------------

  # Create the transformer.
  transformer <- suppressWarnings(do.call(
    power.transform::find_transformation_parameters,
    args=c(list("x"=x_gamma_monotonic),
           parameter_list[[ii]])))

  # Compute the fraction of instances which exceed the threshold.
  power.transform::assess_transformation(
    x = x_gamma_monotonic,
    transformer = transformer)
}

  #### Bi-modal distribution (large separation) --------------------------------

  # Create the transformer.
  transformer <- suppressWarnings(do.call(
    power.transform::find_transformation_parameters,
    args=c(list("x"=x_bimodal_large),
           parameter_list[[ii]])))

  # Compute the fraction of instances which exceed the threshold.
  power.transform::assess_transformation(
    x = x_bimodal_large,
    transformer = transformer)

  #### Bi-modal distribution (intermediate separation) -------------------------

  # Create the transformer.
  transformer <- suppressWarnings(do.call(
    power.transform::find_transformation_parameters,
    args=c(list("x"=x_bimodal_intermediate),
           parameter_list[[ii]])))

  # Compute the fraction of instances which exceed the threshold.
  power.transform::assess_transformation(
    x = x_bimodal_intermediate,
    transformer = transformer)

  #### Bi-modal distribution (small separation) --------------------------------

  # Create the transformer.
  transformer <- suppressWarnings(do.call(
    power.transform::find_transformation_parameters,
    args=c(list("x"=x_bimodal_small),
           parameter_list[[ii]])))

  # Compute the fraction of instances which exceed the threshold.
  power.transform::assess_transformation(
    x = x_bimodal_small,
    transformer = transformer)

  #### Uniform distribution ----------------------------------------------------

  # Create the transformer.
  transformer <- suppressWarnings(do.call(
    power.transform::find_transformation_parameters,
    args=c(list("x"=x_uniform),
           parameter_list[[ii]])))

  # Compute the fraction of instances which exceed the threshold.
  power.transform::assess_transformation(
    x = x_uniform,
    transformer = transformer)


