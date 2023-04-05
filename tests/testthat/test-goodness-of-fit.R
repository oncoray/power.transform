parameter_list <- list()
ii <- 1
for (method in c("box_cox", "yeo_johnson", "none")) {
  for (robust in c(FALSE, TRUE)) {
    for (shift in c(FALSE, TRUE)) {
      parameter_list[[ii]] <- list(
        "method" = method,
        "robust" = robust,
        "shift" = shift
      )

      ii <- ii + 1
    }
  }
}

# Iterate over different sample numbers.
for (n in c(100, 1000, 10000)) {
  # Set seed.
  set.seed(19L)

  # Draw 10000 samples according to a normal distribution.
  x_normal <- stats::rnorm(n)

  # Make exponentially distributed.
  x_exponential <- exp(x_normal)

  # Draw bimodal with large separation
  x_bimodal_large <- c(
    stats::rnorm(n / 2, mean = 0.0),
    stats::rnorm(n / 2, mean = 6.0)
  )

  # Draw bimodal with intermediate separation
  x_bimodal_intermediate <- c(
    stats::rnorm(n / 2, mean = 0.0),
    stats::rnorm(n / 2, mean = 3.0)
  )

  # Draw bimodal with small separation
  x_bimodal_small <- c(
    stats::rnorm(n / 2, mean = 0.0),
    stats::rnorm(n / 2, mean = 1.0)
  )

  # Draw uniform.
  x_uniform <- stats::runif(n, min = 1E-5)

  # Draw Gamma (left skewed)
  x_gamma_left <- stats::rgamma(n, shape = 2.0, scale = 2.0)

  # Draw Gamma (monotonic decreasing)
  x_gamma_monotonic <- stats::rgamma(n, shape = 1.0, scale = 2.0)

  # Iterate over all parameter sets.
  for (ii in seq_along(parameter_list)) {
    testthat::test_that(
      paste0(
        "Assessing transformation goodness generates the correct results. ",
        "(", ii,
        "; n: ", n,
        "; method: ", parameter_list[[ii]]$method,
        "; robust: ", parameter_list[[ii]]$robust,
        "; shift: ", parameter_list[[ii]]$shift, ")"
      ),
      {
        # Normally distributed data --------------------------------------------

        # Create the transformer.
        transformer <- suppressWarnings(do.call(
          power.transform::find_transformation_parameters,
          args = c(
            list("x" = x_normal),
            parameter_list[[ii]]
          )
        ))

        # Compute the p-value.
        p_value <- suppressMessages(power.transform::assess_transformation(
          x = x_normal,
          transformer = transformer
        ))

        if (n == 100) {
          testthat::expect_gt(p_value, 0.10)

        } else {
          testthat::expect_gt(p_value, 0.50)
        }

        # Exponentially distributed data ---------------------------------------

        # Create the transformer.
        transformer <- suppressWarnings(do.call(
          power.transform::find_transformation_parameters,
          args = c(
            list("x" = x_exponential),
            parameter_list[[ii]]
          )
        ))

        # Compute the fraction of instances which exceed the threshold.
        p_value <- suppressMessages(power.transform::assess_transformation(
          x = x_exponential,
          transformer = transformer
        ))

        if (parameter_list[[ii]]$method == "none") {
          testthat::expect_lt(p_value, 0.01)
        } else {
          testthat::expect_gt(p_value, 0.05)
        }


        # Left-skewed gamma distribution ---------------------------------------

        # Create the transformer.
        transformer <- suppressWarnings(do.call(
          power.transform::find_transformation_parameters,
          args = c(
            list("x" = x_gamma_left),
            parameter_list[[ii]]
          )
        ))

        # Compute the fraction of instances which exceed the threshold.
        p_value <- suppressMessages(power.transform::assess_transformation(
          x = x_gamma_left,
          transformer = transformer
        ))

        if (parameter_list[[ii]]$method == "none") {
          testthat::expect_lt(p_value, 0.01)
        } else {
          testthat::expect_gt(p_value, 0.05)
        }

        # Monotonically decreasing gamma distribution --------------------------

        # Create the transformer.
        transformer <- suppressWarnings(do.call(
          power.transform::find_transformation_parameters,
          args = c(
            list("x" = x_gamma_monotonic),
            parameter_list[[ii]]
          )
        ))

        # Compute the fraction of instances which exceed the threshold.
        p_value <- suppressMessages(power.transform::assess_transformation(
          x = x_gamma_monotonic,
          transformer = transformer
        ))

        if (parameter_list[[ii]]$method == "none") {
          testthat::expect_lt(p_value, 0.01)
        } else {
          testthat::expect_gt(p_value, 0.05)
        }

        # Bi-modal distribution (large separation) -----------------------------

        # Create the transformer.
        transformer <- suppressWarnings(do.call(
          power.transform::find_transformation_parameters,
          args = c(
            list("x" = x_bimodal_large),
            parameter_list[[ii]]
          )
        ))

        # Compute the fraction of instances which exceed the threshold.
        p_value <- suppressMessages(power.transform::assess_transformation(
          x = x_bimodal_large,
          transformer = transformer
        ))

        testthat::expect_lt(p_value, 0.001)

        # Bi-modal distribution (intermediate separation) ----------------------

        # Create the transformer.
        transformer <- suppressWarnings(do.call(
          power.transform::find_transformation_parameters,
          args = c(
            list("x" = x_bimodal_intermediate),
            parameter_list[[ii]]
          )
        ))

        # Compute the fraction of instances which exceed the threshold.
        p_value <- suppressMessages(power.transform::assess_transformation(
          x = x_bimodal_intermediate,
          transformer = transformer
        ))

        testthat::expect_lt(p_value, 0.10)

        # Bi-modal distribution (small separation) -----------------------------

        # Create the transformer.
        transformer <- suppressWarnings(do.call(
          power.transform::find_transformation_parameters,
          args = c(
            list("x" = x_bimodal_small),
            parameter_list[[ii]]
          )
        ))

        # Compute the fraction of instances which exceed the threshold.
        p_value <- suppressMessages(power.transform::assess_transformation(
          x = x_bimodal_small,
          transformer = transformer
        ))

        testthat::expect_gt(p_value, 0.01)

        # Uniform distribution -------------------------------------------------

        # Create the transformer.
        transformer <- suppressWarnings(do.call(
          power.transform::find_transformation_parameters,
          args = c(
            list("x" = x_uniform),
            parameter_list[[ii]]
          )
        ))

        # Compute the fraction of instances which exceed the threshold.
        p_value <- suppressMessages(power.transform::assess_transformation(
          x = x_uniform,
          transformer = transformer
        ))

        if (parameter_list[[ii]]$method == "none" && n == 100) {
          testthat::expect_lt(p_value, 0.05)
        } else {
          testthat::expect_gt(p_value, 0.01)
        }
      }
    )
  }
}
