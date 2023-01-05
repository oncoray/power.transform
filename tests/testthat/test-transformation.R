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

# Draw 10000 normally
x <- stats::rnorm(10000)

# Iterate over all parameter sets.
for(ii in seq_along(parameter_list)){
  # Get parameter set.
  parameter_set <- parameter_list[[ii]]

  #### All-positive values -----------------------------------------------------
  x_positive <- exp(x)

  testthat::test_that(
    paste0(
      "Transforming all-positive values generates the correct results. ",
      "(", ii,
      "; method: ", parameter_set$method,
      "; robust: ", parameter_set$robust,
      "; shift: ", parameter_set$shift, ")"),
    {
      # Create the transformer.
      transformer <- do.call(
        power.transform::find_transformation_parameters,
        args=c(list("x"=x_positive),
               parameter_set))

      if(parameter_set$method == "box_cox"){
        if(parameter_set$shift){
          testthat::expect_equal(transformer@lambda, 0.0, tolerance=0.1)
          testthat::expect_equal(transformer@shift, 0.0, tolerance=0.1)

        } else {
          testthat::expect_equal(transformer@lambda, 0.0, tolerance=0.2)
        }


      } else if(parameter_set$method == "yeo_johnson"){
        if(parameter_set$shift){
          if(parameter_set$robust){
            testthat::expect_equal(transformer@lambda, -1.0, tolerance=0.2)

          } else {
            testthat::expect_equal(transformer@lambda, -0.6, tolerance=0.2)
          }

        } else {
          testthat::expect_equal(transformer@lambda, -0.8, tolerance=0.2)
        }
      }

      # Transform values.
      x_transformed <- power.transform::power_transform(
        x = x_positive,
        transformer = transformer)

      # Revert transformation.
      x_reverted <- power.transform::revert_power_transform(
        y = x_transformed,
        transformer = transformer)

      # Test that reverting the transform produces the original values.
      testthat::expect_equal(
        x_positive,
        x_reverted)

      #### Single instance -----------------------------------------------------
      x_transformed_single <- power.transform::power_transform(
        x = 1.5,
        transformer = transformer)

      # Revert transformation.
      x_reverted_single <- power.transform::revert_power_transform(
        y = x_transformed_single,
        transformer = transformer)

      # Test that reverting the transform for the single instance produces the
      # original value.
      testthat::expect_equal(
        1.5,
        x_reverted_single)

      #### Single instance NA --------------------------------------------------
      testthat::expect_warning(
        x_transformed_single <- power.transform::power_transform(
          x = NA_real_,
          transformer = transformer),
        "NA or infinite values were found"
      )

      # Revert transformation.
      x_reverted_single <- power.transform::revert_power_transform(
        y = x_transformed_single,
        transformer = transformer)

      # Test that reverting the transform for the single instance produces the
      # original value.
      testthat::expect_equal(
        NA_real_,
        x_reverted_single)

      #### Single instance non-numeric -----------------------------------------
      testthat::expect_error(
        power.transform::power_transform(
          x = "a",
          transformer = transformer),
        "x does not contain numeric values.")

      #### Single instance negative --------------------------------------------

      if(parameter_set$method == "box_cox"){
        # Box-Cox cannot handle data that fall outside its range.
        testthat::expect_warning(
          x_transformed_single <- power.transform::power_transform(
            x = -100.0,
            transformer = transformer),
          "Box-cox power transforms are only defined for strictly positive values.")

        # Test that the transformed value is now NA.
        testthat::expect_equal(
          NA_real_,
          x_transformed_single)

        # Revert transformation.
        x_reverted_single <- power.transform::revert_power_transform(
          y = x_transformed_single,
          transformer = transformer)

        # Test that reverting the transform for the single instance propagates
        # the NA value.
        testthat::expect_equal(
          NA_real_,
          x_reverted_single)

      } else {
        x_transformed_single <- power.transform::power_transform(
          x = -100.0,
          transformer = transformer)

        # Revert transformation.
        x_reverted_single <- power.transform::revert_power_transform(
          y = x_transformed_single,
          transformer = transformer)

        # Test that reverting the transform for the single instance produces the
        # original value.
        testthat::expect_equal(
          -100.0,
          x_reverted_single)
      }
    }
  )


  #### Some negative values ----------------------------------------------------
  x_part_negative <- exp(x) - 1

  testthat::test_that(
    paste0(
      "Transforming partially negative values generates the correct results. ",
      "(", ii,
      "; method: ", parameter_set$method,
      "; robust: ", parameter_set$robust,
      "; shift: ", parameter_set$shift, ")"),
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
      if(parameter_set$method == "box_cox"){

        if(parameter_set$shift){
          testthat::expect_equal(transformer@lambda, 0.0, tolerance=0.1)
          testthat::expect_equal(transformer@shift, -1.0, tolerance=0.1)

        } else {
          testthat::expect_equal(transformer@lambda, 0.1, tolerance=0.2)
        }

      } else if(parameter_set$method == "yeo_johnson"){
        if(parameter_set$shift){
          if(parameter_set$robust){
            testthat::expect_equal(transformer@lambda, -1.0, tolerance=0.2)

          } else {
            testthat::expect_equal(transformer@lambda, -0.6, tolerance=0.2)
          }

        } else {
          testthat::expect_equal(transformer@lambda, -0.2, tolerance=0.2)
        }
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
    paste0(
      "Transforming completely negative values generates the correct results. ",
      "(", ii,
      "; method: ", parameter_set$method,
      "; robust: ", parameter_set$robust,
      "; shift: ", parameter_set$shift, ")"),
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
      if(parameter_set$method == "box_cox"){

        if(parameter_set$shift){
          testthat::expect_equal(transformer@lambda, 0.0, tolerance=0.1)
          testthat::expect_equal(transformer@shift, -exp(max(x)) - 1E-8, tolerance=1.0)

        } else {
          testthat::expect_equal(transformer@lambda, 0.1, tolerance=0.2)
        }

      } else if(parameter_set$method == "yeo_johnson"){
        if(parameter_set$shift){
          if(parameter_set$robust){
            testthat::expect_equal(transformer@lambda, -1.0, tolerance=0.2)

          } else {
            testthat::expect_equal(transformer@lambda, -0.6, tolerance=0.2)
          }

        } else {
          testthat::expect_equal(transformer@lambda, -4.0, tolerance=0.2)
        }
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


  #### Some NA values ----------------------------------------------------------
  x_some_na <- exp(x)
  x_some_na[c(1,2)] <- NA_real_

  testthat::test_that(
    paste0(
      "Transforming all-positive values, with some NA values, generates the correct results. ",
      "(", ii,
      "; method: ", parameter_set$method,
      "; robust: ", parameter_set$robust,
      "; shift: ", parameter_set$shift, ")"),
    {
      # Create the transformer.
      transformer <- do.call(
        power.transform::find_transformation_parameters,
        args=c(list("x"=x_some_na),
               parameter_set))

      # Check lambda values.
      if(parameter_set$method == "box_cox"){

        if(parameter_set$shift){
          testthat::expect_equal(transformer@lambda, 0.0, tolerance=0.1)
          testthat::expect_equal(transformer@shift, 0.0, tolerance=0.1)

        } else {
          testthat::expect_equal(transformer@lambda, 0.1, tolerance=0.2)
        }

      } else if(parameter_set$method == "yeo_johnson"){
        if(parameter_set$shift){
          if(parameter_set$robust){
            testthat::expect_equal(transformer@lambda, -1.0, tolerance=0.2)

          } else {
            testthat::expect_equal(transformer@lambda, -0.6, tolerance=0.2)
          }

        } else {
          testthat::expect_equal(transformer@lambda, -0.8, tolerance=0.2)
        }
      }

      # Transform values.
      testthat::expect_warning(
        x_transformed <- power.transform::power_transform(
          x = x_some_na,
          transformer = transformer),
        "NA or infinite values were found"
      )

      # Revert transformation.
      x_reverted <- power.transform::revert_power_transform(
        y=x_transformed,
        transformer = transformer)

      # Test that reverting the transform produces the original values.
      testthat::expect_equal(
        x_some_na,
        x_reverted)
    }
  )


  #### Some infinite values ----------------------------------------------------
  x_some_inf <- exp(x)
  x_some_inf[c(1,2)] <- Inf

  testthat::test_that(
    paste0(
      "Transforming all-positive values, with some Inf values, generates the correct results. ",
      "(", ii,
      "; method: ", parameter_set$method,
      "; robust: ", parameter_set$robust,
      "; shift: ", parameter_set$shift, ")"),
    {
      # Create the transformer.
      transformer <- do.call(
        power.transform::find_transformation_parameters,
        args=c(list("x"=x_some_inf),
               parameter_set))

      # Check lambda values.
      if(parameter_set$method == "box_cox"){

        if(parameter_set$shift){
          testthat::expect_equal(transformer@lambda, 0.0, tolerance=0.1)
          testthat::expect_equal(transformer@shift, 0.0, tolerance=0.1)

        } else {
          testthat::expect_equal(transformer@lambda, 0.1, tolerance=0.2)
        }

      } else if(parameter_set$method == "yeo_johnson"){
        if(parameter_set$shift){
          if(parameter_set$robust){
            testthat::expect_equal(transformer@lambda, -1.0, tolerance=0.2)

          } else {
            testthat::expect_equal(transformer@lambda, -0.6, tolerance=0.2)
          }

        } else {
          testthat::expect_equal(transformer@lambda, -0.8, tolerance=0.2)
        }
      }

      # Transform values.
      testthat::expect_warning(
        x_transformed <- power.transform::power_transform(
          x = x_some_inf,
          transformer = transformer),
        "NA or infinite values were found"
      )

      # Revert transformation.
      x_reverted <- power.transform::revert_power_transform(
        y=x_transformed,
        transformer = transformer)

      # The first two values should now be NA instead of Inf, because the
      # transformation routines replace this by NA.
      x_some_inf[c(1,2)] <- NA_real_

      # Test that reverting the transform produces the expected values.
      testthat::expect_equal(
        x_some_inf,
        x_reverted,
        tolerance=1E-8)
    }
  )


  #### All NA values -----------------------------------------------------------
  x_all_na <- rep_len(NA_real_, 1000L)

  testthat::test_that(
    paste0(
      "Transforming all-na values generates the correct results. ",
      "(", ii,
      "; method: ", parameter_set$method,
      "; robust: ", parameter_set$robust,
      "; shift: ", parameter_set$shift, ")"),
    {
      # Creating the transformer should throw an error.
      testthat::expect_error(
        do.call(
          power.transform::find_transformation_parameters,
          args=c(list("x"=x_all_na),
                 parameter_set)),
        "x only contains NA or inf values.")
    }
  )


  #### All infinite values -----------------------------------------------------
  x_all_inf <- rep_len(Inf, 1000L)

  testthat::test_that(
    paste0(
      "Transforming all-infinite values generates the correct results. ",
      "(", ii,
      "; method: ", parameter_set$method,
      "; robust: ", parameter_set$robust,
      "; shift: ", parameter_set$shift, ")"),
    {
      # Creating the transformer should throw an error.
      testthat::expect_error(
        do.call(
          power.transform::find_transformation_parameters,
          args=c(list("x"=x_all_inf),
                 parameter_set)),
        "x only contains NA or inf values.")
    }
  )


  #### Non-numerical (character) values ----------------------------------------
  x_all_char <- letters[round(x - min(x) + 1)]

  testthat::test_that(
    paste0(
      "Transforming non-numeric values generates the correct results. ",
      "(", ii,
      "; method: ", parameter_set$method,
      "; robust: ", parameter_set$robust,
      "; shift: ", parameter_set$shift, ")"),
    {
      # Creating the transformer should throw an error.
      testthat::expect_error(
        do.call(
          power.transform::find_transformation_parameters,
          args=c(list("x"=x_all_char),
                 parameter_set)),
        "x does not contain numeric values.")
    }
  )


  #### Categorical values ------------------------------------------------------
  x_categorical <- letters[round(x - min(x) + 1)]
  x_categorical <- factor(
    x_categorical,
    levels=sort(unique(x_categorical)))

  testthat::test_that(
    paste0(
      "Transforming categorical values generates the correct results. ",
      "(", ii,
      "; method: ", parameter_set$method,
      "; robust: ", parameter_set$robust,
      "; shift: ", parameter_set$shift, ")"),
    {
      # Creating the transformer should throw an error.
      testthat::expect_error(
        do.call(
          power.transform::find_transformation_parameters,
          args=c(list("x"=x_categorical),
                 parameter_set)),
        "x is categorical, and power transformations are not applicable.")
    }
  )


  #### Single value ------------------------------------------------------------
  x_single_value <- 1.0

  testthat::test_that(
    paste0(
      "Transforming single values generates the correct results. ",
      "(", ii,
      "; method: ", parameter_set$method,
      "; robust: ", parameter_set$robust,
      "; shift: ", parameter_set$shift, ")"),
    {
      # Creating the transformer should throw a warning, and produce a
      # transformationNone object instead.
      if(parameter_set$method == "none"){
        testthat::expect_no_warning(
          transformer <- do.call(
            power.transform::find_transformation_parameters,
            args=c(list("x"=x_single_value),
                   parameter_set)))

      } else {
        testthat::expect_warning(
          transformer <- do.call(
            power.transform::find_transformation_parameters,
            args=c(list("x"=x_single_value),
                   parameter_set)),
          "x contains three or fewer unique values")
      }

      testthat::expect_s4_class(transformer, "transformationNone")
    }
  )


  #### Few unique values (<= 3) ------------------------------------------------
  x_three_unique <- ceiling(stats::runif(1000L) * 3)

  testthat::test_that(
    paste0(
      "Transforming vector with 3 unique values generates the correct results. ",
      "(", ii,
      "; method: ", parameter_set$method,
      "; robust: ", parameter_set$robust,
      "; shift: ", parameter_set$shift, ")"),
    {
      # Creating the transformer should throw a warning, and produce a
      # transformationNone object instead.
      if(parameter_set$method == "none"){
        testthat::expect_no_warning(
          transformer <- do.call(
            power.transform::find_transformation_parameters,
            args=c(list("x"=x_three_unique),
                   parameter_set)))
      } else {
        testthat::expect_warning(
          transformer <- do.call(
            power.transform::find_transformation_parameters,
            args=c(list("x"=x_three_unique),
                   parameter_set)),
          "x contains three or fewer unique values")
      }

      testthat::expect_s4_class(transformer, "transformationNone")
    }
  )


  #### Few unique values (< 10) ------------------------------------------------
  x_few_unique <- ceiling(x - min(x) + 1E-8)

  testthat::test_that(
    paste0(
      "Transforming vector with fewer than 10 unique values generates the correct results. ",
      "(", ii,
      "; method: ", parameter_set$method,
      "; robust: ", parameter_set$robust,
      "; shift: ", parameter_set$shift, ")"),
    {
      # Creating the transformer should throw a warning, but otherwise function
      # normally.
      testthat::expect_warning(
        transformer <- do.call(
          power.transform::find_transformation_parameters,
          args=c(list("x"=x_few_unique),
                 parameter_set)),
        "x contains ten or fewer unique values")

      # Check lambda values.
      if(parameter_set$method == "box_cox"){

        if(parameter_set$robust){
          if(parameter_set$shift){
            # Very close to logarithmic transformation.
            testthat::expect_equal(transformer@lambda, 1.0, tolerance=0.2)
          } else {
            testthat::expect_equal(transformer@lambda, 0.6, tolerance=0.2)
          }

        } else {
          testthat::expect_equal(transformer@lambda, 1.0, tolerance=0.2)
        }

      } else if(parameter_set$method == "yeo_johnson"){

        if(parameter_set$robust){
          if(parameter_set$shift){
            testthat::expect_equal(transformer@lambda, 1.0, tolerance=0.2)
          } else {
            testthat::expect_equal(transformer@lambda, 0.6, tolerance=0.2)
          }

        } else {
          testthat::expect_equal(transformer@lambda, 1.0, tolerance=0.2)
        }
      }

      # Transform values.
      x_transformed <- power.transform::power_transform(
        x=x_few_unique,
        transformer = transformer)

      # Revert transformation.
      x_reverted <- power.transform::revert_power_transform(
        y=x_transformed,
        transformer = transformer)

      # Test that reverting the transform produces the original values.
      testthat::expect_equal(
        x_few_unique,
        x_reverted)
    }
  )
}
