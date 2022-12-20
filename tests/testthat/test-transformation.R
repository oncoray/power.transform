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

# Draw 1000 normally
x <- stats::rnorm(1000)

# Iterate over all parameter sets.
for(parameter_set in parameter_list){

  #### All-positive values -----------------------------------------------------
  x_positive <- exp(x)

  testthat::test_that(
    paste0("Transforming all-positive values generates the correct results. "),
    {
      # Create the transformer.
      transformer <- do.call(
        power.transform::find_transformation_parameters,
        args=c(list("x"=x_positive),
               parameter_set))

      if(parameter_set$method != "none"){
        testthat::expect_equal(transformer@lambda, 0.0, tolerance=0.1)
      }

      # Transform values.
      x_transformed <- power.transform::power_transform(
        x=x_positive,
        transformer = transformer)

      # Revert transformation.
      x_reverted <- power.transform::revert_power_transform(
        y=x_transformed,
        transformer = transformer)

      # Test that reverting the transform produces the original values.
      testthat::expect_equal(
        x_positive,
        x_reverted)
    }
  )

  #### Some negative values ----------------------------------------------------
  x_part_negative <- exp(x) - 1

  testthat::test_that(
    paste0("Transforming partially negative values generates the correct results. "),
    {
      # Create the transformer.
      if(parameter_set$method == "box_cox" && !parameter_set$shift){
        testthat::expect_warning(
          transformer <- do.call(
            power.transform::find_transformation_parameters,
            args=c(list("x"=x_part_negative),
                   parameter_set)),
          "Box-cox power transforms are only defined for strictly positive values.")

      } else {
        transformer <- do.call(
          power.transform::find_transformation_parameters,
          args=c(list("x"=x_part_negative),
                 parameter_set))
      }

      # Check lambda values.
      if(parameter_set$method != "none"){
        testthat::expect_equal(transformer@lambda, 0.0, tolerance=0.2)
      }

      # Transform values.
      x_transformed <- power.transform::power_transform(
        x=x_part_negative,
        transformer = transformer)

      # Revert transformation.
      x_reverted <- power.transform::revert_power_transform(
        y=x_transformed,
        transformer = transformer)

      # Test that reverting the transform produces the original values.
      testthat::expect_equal(
        x_part_negative,
        x_reverted)
    }
  )


  #### All-negative values -----------------------------------------------------
  x_negative <- exp(x) - exp(max(x)) - 1E-8

  testthat::test_that(
    paste0("Transforming partially negative values generates the correct results. "),
    {
      # Create the transformer.
      if(parameter_set$method == "box_cox" && !parameter_set$shift){
       testthat::expect_warning(
          transformer <- do.call(
            power.transform::find_transformation_parameters,
            args=c(list("x"=x_negative),
                   parameter_set)),
          "Box-cox power transforms are only defined for strictly positive values.")

      } else {
        transformer <- do.call(
          power.transform::find_transformation_parameters,
          args=c(list("x"=x_negative),
                 parameter_set))
      }

      # Check lambda values.
      if(parameter_set$method != "none"){
        testthat::expect_equal(transformer@lambda, 0.0, tolerance=0.2)
      }

      # Transform values.
      x_transformed <- power.transform::power_transform(
        x=x_negative,
        transformer = transformer)

      # Revert transformation.
      x_reverted <- power.transform::revert_power_transform(
        y=x_transformed,
        transformer = transformer)

      # Test that reverting the transform produces the original values.
      testthat::expect_equal(
        x_negative,
        x_reverted)
    }
  )

  # Some NA values.
  x_some_na <- x_positive

  # Some infinite values.

  # All NA values.

  # All infinite values.

  # Non-numerical (character) values.

  # Categorical values.

  # Single value.

  # Few unique values (<= 3)

  # Few unique values (< 10)
}
