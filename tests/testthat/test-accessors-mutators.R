
for (method in c("box_cox", "yeo_johnson", "none")) {
  testthat::test_that(
    "Accessing and changing lambda parameters functions correctly.", {

      # Generate data.
      x <- stats::rnorm(n = 1000L)

      # Create the transformer. Note that the data are not normally distributed.
      transformer <- power.transform::find_transformation_parameters(
        x = x,
        method = method,
        estimation_method = "mle"
      )

      if (method == "none") {
        testthat::expect_warning(
          transformer_new <- power.transform::set_lambda(
            object = transformer,
            lambda = 2.0),
          class = "power_transform_no_attribute"
        )

        testthat::expect_equal(transformer, transformer_new)

        testthat::expect_warning(
          lambda_value <- power.transform::get_lambda(
            object = transformer_new),
          class = "power_transform_no_attribute"
        )

        testthat::expect_equal(lambda_value, NA_real_)

      } else {
        transformer_new <- power.transform::set_lambda(
          object = transformer,
          lambda = 2.0)

        lambda_value <- power.transform::get_lambda(
          object = transformer_new)

        testthat::expect_equal(lambda_value, 2.0)
      }
    }
  )

  testthat::test_that(
    "Accessing and changing shift parameters functions correctly.", {

      # Generate data.
      x <- stats::rnorm(n = 1000L)

      # Create the transformer. Note that the data are not normally distributed.
      transformer <- power.transform::find_transformation_parameters(
        x = x,
        method = method,
        estimation_method = "mle"
      )

      if (method == "none") {
        testthat::expect_warning(
          transformer_new <- power.transform::set_shift(
            object = transformer,
            shift = 3.0),
          class = "power_transform_no_attribute"
        )

        testthat::expect_equal(transformer, transformer_new)

        testthat::expect_warning(
          shift_value <- power.transform::get_shift(
            object = transformer_new),
          class = "power_transform_no_attribute"
        )

        testthat::expect_equal(shift_value, NA_real_)

      } else {
        transformer_new <- power.transform::set_shift(
          object = transformer,
          shift = 3.0)

        shift_value <- power.transform::get_shift(
          object = transformer_new)

        testthat::expect_equal(shift_value, 3.0)
      }
    }
  )
}

