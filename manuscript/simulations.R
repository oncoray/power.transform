.get_shifted_scaled_distribution_data <- function(manuscript_dir, main_manuscript) {
  # Set seed.
  set.seed(19L)

  if (main_manuscript) {
    file_name <- file.path(manuscript_dir, "shifted_scaled_distributions_plot_main.RDS")
    estimation_methods <- "mle"
  } else {
    file_name <- file.path(manuscript_dir, "shifted_scaled_distributions_plot_appendix.RDS")
    estimation_methods <- setdiff(
      power.transform:::..estimators_all(),
      power.transform:::..estimators_raymaekers_robust()
    )
  }

  if (!file.exists(file_name)) {
    # Normal distribution.
    x_normal <- power.transform::ragn(10000L, location = 0, scale = 1 / sqrt(2), alpha = 0.5, beta = 2)

    # Right skewed data
    x_right_skewed <- power.transform::ragn(10000L, location = 0, scale = 1 / sqrt(2), alpha = 0.2, beta = 2)

    # Left skewed data
    x_left_skewed <- power.transform::ragn(10000L, location = 0, scale = 1 / sqrt(2), alpha = 0.8, beta = 2)

    # generator ----------------------------------------------------------------
    generate_experiment_data <- coro::generator(
      function(
    x_normal,
    x_right_skewed,
    x_left_skewed,
    estimation_methods
      ) {
        shift_range <- 10^seq(from = 0, to = 6, by = 0.1)
        scale_range <- 10^seq(from = 0, to = 6, by = 0.1)

        for (distribution in c("normal", "right-skewed", "left-skewed")) {
          if (distribution == "normal") {
            x <- x_normal
          } else if (distribution == "right-skewed") {
            x <- x_right_skewed
          } else if (distribution == "left-skewed") {
            x <- x_left_skewed
          }

          for (method in c("box_cox", "yeo_johnson")) {
            for (invariant in c(FALSE, TRUE)) {
              for (robust in c(FALSE)) {
                for (estimation_method in estimation_methods) {

                  # Skip if robust, but not invariant.
                  if (robust && !invariant) next

                  if (!invariant && !robust) {
                    version <- "conventional"
                  } else if (!robust) {
                    version <- "invariant"
                  } else {
                    version <- "robust invariant"
                  }

                  # Iterate over shift range.
                  for (d in shift_range) {
                    yield(list(
                      "x" = x + d,
                      "d" = d,
                      "s" = NA_real_,
                      "distribution" = distribution,
                      "method" = method,
                      "invariant" = invariant,
                      "robust" = robust,
                      "version" = version,
                      "estimation_method" = estimation_method,
                      "data_type" = "shifted"
                    ))
                  }

                  # Iterate over scale range.
                  for (s in scale_range) {
                    yield(list(
                      "x" = x * s,
                      "d" = NA_real_,
                      "s" = s,
                      "distribution" = distribution,
                      "method" = method,
                      "invariant" = invariant,
                      "robust" = robust,
                      "version" = version,
                      "estimation_method" = estimation_method,
                      "data_type" = "scaled"
                    ))
                  }
                }
              }
            }
          }
        }
      }
    )

    .compute_lambda <- function(parameter_set) {

      if (parameter_set$version == "conventional") {
        # Conventional transformer

        transformer <- suppressWarnings(
          power.transform::find_transformation_parameters(
            x = parameter_set$x,
            method = parameter_set$method,
            robust = parameter_set$robust,
            invariant = parameter_set$invariant,
            estimation_method = parameter_set$estimation_method
          )
        )
      } else {
        transformer <- suppressWarnings(
          power.transform::find_transformation_parameters(
            x = parameter_set$x,
            method = parameter_set$method,
            robust = parameter_set$robust,
            invariant = parameter_set$invariant,
            estimation_method = parameter_set$estimation_method
          )
        )
      }


      return(
        data.table::data.table(
          "distribution" = parameter_set$distribution,
          "method" = parameter_set$method,
          "estimation_method" = parameter_set$estimation_method,
          "version" = parameter_set$version,
          "shift" = log10(parameter_set$d),
          "scale" = log10(parameter_set$s),
          "lambda" = transformer@lambda,
          "x_0" = transformer@shift,
          "s" = transformer@scale,
          "data_type" = parameter_set$data_type
        )
      )
    }

    # computations -------------------------------------------------------------

    # Generate all experiments.
    experiments <- coro::collect(generate_experiment_data(
      x_normal = x_normal,
      x_right_skewed = x_right_skewed,
      x_left_skewed = x_left_skewed,
      estimation_methods = estimation_methods
    ))

    # data <- lapply(
    #   X = experiments,
    #   FUN = .compute_lambda)

    # Start cluster
    cl <- parallel::makeCluster(16L)

    # Compute all data in parallel.
    data <- parallel::parLapply(
      cl = cl,
      X = experiments,
      fun = .compute_lambda
    )

    # Stop cluster.
    parallel::stopCluster(cl)

    # Combine data into a single table.
    data <- data.table::rbindlist(data)

    # Save to file.
    saveRDS(data, file_name)

  } else {
    data <- readRDS(file_name)
  }

  # Update distribution, method and version to factors.
  data$distribution <- factor(
    x = data$distribution,
    levels = c("normal", "right-skewed", "left-skewed")
  )
  data$method <- factor(
    x = data$method,
    levels = c("box_cox", "yeo_johnson"),
    labels = c("Box-Cox", "Yeo-Johnson")
  )
  data$version <- factor(
    x = data$version,
    levels = c("conventional", "invariant")
  )
  data$estimation_method <- factor(
    x = data$estimation_method,
    levels = c("mle", "anderson_darling", "cramer_von_mises", "jarque_bera", "dagostino"),
    labels = c("MLE", "Anderson-Darling", "Cramér-von Mises", "Jarque-Bera", "D'Agostino")
  )
  data$data_type <- factor(
    x = data$data_type,
    levels = c("shifted", "scaled")
  )

  return(data)
}



.get_shifted_outlier_plot_data <- function(manuscript_dir, offset = 0.0) {
  # Set seed.
  set.seed(19L)

  if (!file.exists(file.path(manuscript_dir, "shifted_outlier_plot.RDS"))) {
    # Normal distribution.
    x <- power.transform::ragn(1000L, location = offset, scale = 1 / sqrt(2), alpha = 0.5, beta = 2)

    # generator ----------------------------------------------------------------
    generate_experiment_data <- coro::generator(
      function(x, offset) {
        outlier_shift_range <- 10^seq(from = -1, to = 6, by = 0.1)

        for (d in outlier_shift_range) {
          for (method in c("box_cox", "yeo_johnson")) {
            yield(list(
              "x" = c(x, offset + d),
              "d" = d,
              "method" = method,
              "invariant" = FALSE,
              "robust" = FALSE
            ))
          }
        }
      }
    )

    .compute_lambda <- function(parameter_set) {
      # Create transformer object.
      transformer <- suppressWarnings(
        power.transform::find_transformation_parameters(
          x = parameter_set$x,
          method = parameter_set$method,
          robust = parameter_set$robust,
          invariant = parameter_set$invariant
        )
      )

      return(
        data.table::data.table(
          "method" = parameter_set$method,
          "d" = log10(parameter_set$d),
          "lambda" = transformer@lambda,
          "x_0" = transformer@shift
        )
      )
    }

    # computations -------------------------------------------------------------

    # Generate all experiments.
    experiments <- coro::collect(generate_experiment_data(x=x, offset = offset))

    data <- lapply(
      X = experiments,
      FUN = .compute_lambda
    )

    # # Start cluster
    # cl <- parallel::makeCluster(16L)
    #
    # # Compute all data in parallel.
    # data <- parallel::parLapply(
    #   cl = cl,
    #   X = experiments,
    #   fun = .compute_lambda
    # )
    #
    # # Stop cluster.
    # parallel::stopCluster(cl)

    # Combine data into a single table.
    data <- data.table::rbindlist(data)

    # Save to file.
    saveRDS(data, file.path(manuscript_dir, "shifted_outlier_plot.RDS"))

  } else {
    data <- readRDS(file.path(manuscript_dir, "shifted_outlier_plot.RDS"))
  }

  # Update method to factor.
  data$method <- factor(
    x = data$method,
    levels = c("box_cox", "yeo_johnson"),
    labels = c("Box-Cox", "Yeo-Johnson")
  )

  return(data)
}



.get_optimised_weighting_function_parameters <- function(manuscript_dir, side = "both") {
  # generator ------------------------------------------------------------------
  experiment_args <- coro::generator(
    function(
        manuscript_dir,
        side,
        transformation_methods = NULL,
        estimation_methods = NULL) {
      if (side == "both") {
        file_dir <- "robustness_comparison"
      } else if (side == "right") {
        file_dir <- "robustness_comparison_right"
      }

      if (is.null(transformation_methods)) {
        transformation_methods <- c("box_cox", "yeo_johnson")
      }

      if (is.null(estimation_methods)) {
        estimation_methods <- setdiff(power.transform:::..estimators_all(), power.transform:::..estimators_raymaekers_robust())
      }

      for (estimation_method in estimation_methods) {
        for (transformation_method in transformation_methods) {
          # Non-robust transformation.
          fun_args <- list(
            "method" = transformation_method,
            "estimation_method" = estimation_method,
            "robust" = FALSE,
            "shift" = TRUE
          )

          opt_args <- list()

          file_name <- file.path(
            manuscript_dir,
            file_dir,
            paste0(
              transformation_method, "_",
              estimation_method, "_",
              "non_robust.RDS"
            )
          )

          if (!file.exists(file_name)) {
            # Only yield arguments if the file does not exist.
            yield(list(
              "name" = "non-robust",
              "method" = transformation_method,
              "estimation_method" = estimation_method,
              "file_name" = file_name,
              "fun_args" = fun_args,
              "opt_args" = opt_args
            ))
          } else {
            yield(file_name)
          }
        }

        for (robustness_source in c("empirical_probability", "transformed", "residual")) {
          def_opt_limits <- list("k1" = c(0.0, 10.0), "k2" = c(0.0, 10.0))

          if (robustness_source == "empirical_probability") {
            def_opt_limits$k1[2] <- 1.0
            def_opt_limits$k2[2] <- 1.0

            def_opt_init <- list("k1" = 0.80, "k2" = 0.95)
          } else if (robustness_source == "transformed") {
            def_opt_init <- list("k1" = 1.28, "k2" = 1.96)
          } else if (robustness_source == "residual") {
            def_opt_init <- list("k1" = 0.50, "k2" = 1.0)
          }

          for (transformation_method in transformation_methods) {
            for (robustness_weighting_function in c("step", "triangle", "cosine")) {
              fun_args <- list(
                "method" = transformation_method,
                "estimation_method" = estimation_method,
                "robust" = TRUE,
                "shift" = TRUE,
                "weighting_function" = paste0(robustness_source, "_", robustness_weighting_function),
                "backup_use_default" = FALSE
              )

              opt_limits <- def_opt_limits
              opt_init <- def_opt_init

              # Step-weighting only has k1 as a parameter.
              if (robustness_weighting_function == "step") {
                opt_limits$k2 <- NULL
                opt_init$k2 <- NULL
              }

              file_name <- file.path(
                manuscript_dir,
                file_dir,
                paste0(
                  transformation_method, "_",
                  estimation_method, "_",
                  robustness_source, "_",
                  robustness_weighting_function, ".RDS"
                )
              )

              if (!file.exists(file_name)) {
                # Only yield arguments if the file does not exist.
                yield(list(
                  "name" = paste0(robustness_source, "-", robustness_weighting_function),
                  "method" = transformation_method,
                  "estimation_method" = estimation_method,
                  "file_name" = file_name,
                  "fun_args" = fun_args,
                  "opt_args" = list("initial" = opt_init, "limits" = opt_limits)
                ))
              } else {
                yield(file_name)
              }
            }
          }
        }
      }
    }
  )


  # Helper function for computing the target lambda value.
  .compute_target_lambda <- function(
      x,
      transformation_methods = NULL,
      estimation_methods = NULL) {
    if (is.null(transformation_methods)) {
      transformation_methods <- c("box_cox", "yeo_johnson")
    }

    if (is.null(estimation_methods)) {
      estimation_methods <- setdiff(power.transform:::..estimators_all(), power.transform:::..estimators_raymaekers_robust())
    }

    lambda_list <- list()
    ii <- 1L
    for (transformation_method in transformation_methods) {
      for (estimation_method in estimation_methods) {
        # Create transformer object.
        transformer <- suppressWarnings(
          power.transform::find_transformation_parameters(
            x = x,
            method = transformation_method,
            estimation_method = estimation_method,
            robust = FALSE,
            shift = TRUE
          )
        )

        lambda_list[[ii]] <- data.table::data.table(
          "method" = transformation_method,
          "estimation_method" = estimation_method,
          "lambda" = transformer@lambda
        )

        ii <- ii + 1L
      }
    }

    return(data.table::rbindlist(lambda_list))
  }


  # Helper function for population distributions with outliers.
  .populate_outliers <- function(
      x,
      n,
      alpha,
      beta,
      target_lambda,
      outlier_fraction,
      side) {
    # Set parameters.
    parameters <- list(
      "n" = n,
      "alpha" = alpha,
      "beta" = beta,
      "target_lambda" = target_lambda
    )

    # Add outliers.
    x <- lapply(
      outlier_fraction,
      function(k, x, parameters, side) {
        # Set parameters.
        parameters$k <- k

        # Compute interquartile range.
        interquartile_range <- stats::IQR(x)

        # Compute upper and lower quartiles.
        q_lower <- stats::quantile(x, probs = 0.25, names = FALSE)
        q_upper <- stats::quantile(x, probs = 0.75, names = FALSE)

        # Set data where the outliers will be copied into.
        x_outlier <- x

        if (k != 0.0) {
          n_draw <- ceiling(k * length(x))

          # Generate outlier values that are smaller than Q1 - 1.5 IQR or larger
          # than Q3 + 1.5 IQR.
          if (side == "both") {
            x_random <- stats::runif(n_draw, min = -2.0, max = 2.0)
          } else if (side == "right") {
            x_random <- stats::runif(n_draw, min = 0.0, max = 2.0)
          } else if (side == "left") {
            x_random <- stats::runif(n_draw, min = -2.0, max = 0.0)
          } else {
            stop(paste0("side was not recognised: ", side))
          }

          outlier <- numeric(n_draw)
          if (any(x_random < 0)) {
            outlier[x_random < 0] <- q_lower - 1.5 * interquartile_range + x_random[x_random < 0] * interquartile_range
          }

          if (any(x_random >= 0)) {
            outlier[x_random >= 0] <- q_upper + 1.5 * interquartile_range + x_random[x_random >= 0] * interquartile_range
          }

          # Randomly insert outlier values.
          x_outlier[sample(seq_along(x), size = n_draw, replace = FALSE)] <- outlier
        }

        return(list(
          "x" = x_outlier,
          "parameters" = parameters
        ))
      },
      x = x,
      parameters = parameters,
      side = side
    )

    return(x)
  }

  # Helper function for computing transformer parameters under optimisation
  # constraints.
  .compute_robust_lambda <- function(
      x,
      fun_args,
      opt_args) {
    # Custom parser.
    ..outlier_parser <- function(
        x,
        fun_args,
        opt_args) {
      # Prevent notes
      method <- estimation_method <- NULL

      # Create transformer.
      transformer <- suppressWarnings(do.call(
        power.transform::find_transformation_parameters,
        args = c(
          list("x" = x$x),
          list("weighting_function_parameters" = opt_args),
          fun_args
        )
      ))

      target_lambda <- x$parameters$target_lambda[method == fun_args$method & estimation_method == fun_args$estimation_method]$lambda

      parameter_data <- c(
        opt_args,
        list(
          "method" = fun_args$method,
          "estimation_method" = fun_args$estimation_method,
          "weight_method" = fun_args$weighting_function,
          "lambda" = transformer@lambda,
          "target_lambda" = target_lambda,
          "shift" = transformer@shift
        )
      )

      return(parameter_data)
    }

    parameter_data <- lapply(
      x,
      ..outlier_parser,
      fun_args = fun_args,
      opt_args = opt_args
    )

    return(parameter_data)
  }

  # Outer optimisation function used to optimise the weighting function parameters
  # over a dataset.
  .optimisation_outer <- function(
      experiment,
      data,
      verbose = FALSE) {
    set.seed(9)

    results_list <- list(
      "name" = experiment$name,
      "method" = experiment$method,
      "estimation_method" = experiment$estimation_method
    )

    # Run optimiser.
    if (length(experiment$opt_args) > 0) {
      x0 <- unlist(experiment$opt_args$initial)

      lower <- experiment$opt_args$limits$k1[1]
      upper <- experiment$opt_args$limits$k1[2]

      if (length(experiment$opt_args$limits) > 1) {
        lower <- c(lower, experiment$opt_args$limits$k2[1])
        upper <- c(upper, experiment$opt_args$limits$k2[2])
      }

      h <- nloptr::sbplx(
        x0 = x0,
        fn = .optimisation_inner,
        lower = lower,
        upper = upper,
        control = list("xtol_rel" = 1e-2, "maxeval" = 50),
        x = data,
        fun_args = experiment$fun_args,
        verbose = verbose
      )

      results_list <- c(results_list, list("k1" = h$par[1]))
      if (length(h$par) > 1) results_list <- c(results_list, list("k2" = h$par[2]))
      results_list <- c(results_list, list("value" = h$value))
    } else {
      value <- .optimisation_inner(
        opt_param = list(),
        x = data,
        fun_args = experiment$fun_args
      )

      results_list <- c(results_list, list("value" = value))
    }

    # Save to file.
    saveRDS(results_list, file = experiment$file_name)

    return(invisible(TRUE))
  }

  # Inner optimisation function that sets the loss.
  .optimisation_inner <- function(
      opt_param,
      x,
      fun_args,
      verbose = FALSE) {
    # Parse options.
    opt_args <- list()
    if (length(opt_param) >= 1) opt_args <- c(opt_args, list("k1" = opt_param[1]))
    if (length(opt_param) >= 2) opt_args <- c(opt_args, list("k2" = opt_param[2]))

    parameter_data <- .compute_robust_lambda(
      x = x,
      fun_args = fun_args,
      opt_args = opt_args
    )

    lambda_error <- sapply(
      parameter_data,
      function(x) (abs(x$lambda - x$target_lambda))
    )

    if (verbose) {
      message <- paste0("target: ", mean(lambda_error))
      if (length(opt_args) > 0) message <- c(message, paste0("k1: ", opt_param[1]))
      if (length(opt_args) > 1) message <- c(message, paste0("k2: ", opt_param[2]))
      cat(paste0(message, collapse = "; "))
      cat("\n")
    }

    return(mean(lambda_error))
  }

  # computations ---------------------------------------------------------------

  # Select experiments to conduct by filtering out experiments that are file
  # names.
  experiments <- coro::collect(experiment_args(
    manuscript_dir = manuscript_dir,
    side = side
  ))
  experiments <- experiments[!sapply(experiments, is.character)]

  if (length(experiments) > 0) {
    set.seed(95)

    # Generate table of asymmetric generalised normal distribution parameters.
    n_distributions <- 100L

    # Generate alpha, beta and n.
    n <- stats::runif(n = n_distributions, min = 2, max = 4)
    n <- ceiling(10^n)

    alpha <- stats::runif(n = n_distributions, min = 0.01, max = 0.99)
    beta <- stats::runif(n = n_distributions, min = 1.00, max = 5.00)

    # Generate corresponding distributions.
    x <- mapply(
      power.transform::ragn,
      n = n,
      alpha = alpha,
      beta = beta
    )

    # Compute lambda values without weighting.
    target_lambda <- lapply(x, .compute_target_lambda)

    # Add outliers, between 0.00 and 0.10.
    n_outliers <- 10
    outlier_fraction <- (seq_len(n_outliers) - 1)^2 / ((n_outliers - 1)^2 * 10)

    # Generate outliers in the data.
    x <- mapply(
      FUN = .populate_outliers,
      x = x,
      n = n,
      alpha = alpha,
      beta = beta,
      target_lambda = target_lambda,
      MoreArgs = list(
        "outlier_fraction" = outlier_fraction,
        "side" = side
      ),
      SIMPLIFY = FALSE,
      USE.NAMES = FALSE
    )

    x <- unlist(x, recursive = FALSE)

    # lapply(experiments, .optimisation_outer, data=x, verbose=TRUE)

    # Start cluster
    cl <- parallel::makeCluster(18L)

    parallel::clusterExport(
      cl = cl,
      varlist = c(".compute_robust_lambda", ".optimisation_inner"),
      envir = environment()
    )

    parallel::parLapplyLB(
      cl = cl,
      X = experiments,
      fun = .optimisation_outer,
      data = x,
      chunk.size = 1L
    )

    # Stop cluster.
    parallel::stopCluster(cl)
  }

  # Collect all files
  experiments <- coro::collect(experiment_args(
    manuscript_dir = manuscript_dir,
    side = side
  ))
  weighting_parameters <- lapply(experiments, readRDS)
  weighting_parameters <- lapply(weighting_parameters, data.table::as.data.table)
  weighting_parameters <- data.table::rbindlist(weighting_parameters, fill = TRUE)

  return(weighting_parameters)
}



.get_transformation_parameters <- function(manuscript_dir, side = "both") {
  if (side == "both") {
    file_name <- file.path(manuscript_dir, "robustness_comparison_marginal_parameter_data.RDS")
  } else if (side == "right") {
    file_name <- file.path(manuscript_dir, "robustness_comparison_marginal_parameter_data_right.RDS")
  }

  if (!file.exists(file_name)) {
    # generator ----------------------------------------------------------------
    experiment_args <- coro::generator(
      function(manuscript_dir,
               side,
               transformation_methods = NULL,
               estimation_methods = NULL) {
        if (side == "both") {
          file_dir <- "robustness_comparison"
        } else if (side == "right") {
          file_dir <- "robustness_comparison_right"
        }

        if (is.null(transformation_methods)) {
          transformation_methods <- c("box_cox", "yeo_johnson")
        }

        if (is.null(estimation_methods)) {
          estimation_methods <- setdiff(power.transform:::..estimators_all(), power.transform:::..estimators_raymaekers_robust())
        }

        for (estimation_method in estimation_methods) {
          for (transformation_method in transformation_methods) {
            # Non-robust transformation.
            fun_args <- list(
              "method" = transformation_method,
              "estimation_method" = estimation_method,
              "robust" = FALSE,
              "shift" = TRUE
            )

            file_name <- file.path(
              manuscript_dir,
              file_dir,
              paste0(
                transformation_method, "_",
                estimation_method, "_",
                "non_robust.RDS"
              )
            )

            # Only yield arguments if the file does not exist.
            yield(list(
              "name" = "non-robust",
              "method" = transformation_method,
              "estimation_method" = estimation_method,
              "file_name" = file_name,
              "fun_args" = fun_args
            ))
          }

          for (robustness_source in c("empirical_probability", "transformed", "residual")) {
            for (transformation_method in transformation_methods) {
              for (robustness_weighting_function in c("step", "triangle", "cosine")) {
                file_name <- file.path(
                  manuscript_dir,
                  file_dir,
                  paste0(
                    transformation_method, "_",
                    estimation_method, "_",
                    robustness_source, "_",
                    robustness_weighting_function, ".RDS"
                  )
                )

                optimisation_settings <- readRDS(file_name)
                weighting_function_parameters <- list("k1" = optimisation_settings$k1)
                if (!is.null(optimisation_settings$k2)) {
                  weighting_function_parameters <- c(
                    weighting_function_parameters,
                    list("k2" = optimisation_settings$k2)
                  )
                }

                fun_args <- list(
                  "method" = transformation_method,
                  "estimation_method" = estimation_method,
                  "robust" = TRUE,
                  "shift" = TRUE,
                  "weighting_function" = paste0(robustness_source, "_", robustness_weighting_function),
                  "weighting_function_parameters" = weighting_function_parameters
                )

                # Only yield arguments if the file does not exist.
                yield(list(
                  "name" = paste0(robustness_source, "-", robustness_weighting_function),
                  "method" = transformation_method,
                  "estimation_method" = estimation_method,
                  "file_name" = file_name,
                  "fun_args" = fun_args
                ))
              }
            }
          }
        }
      }
    )



    # Helper function for computing the target lambda value.
    .compute_target_lambda <- function(x,
                                       transformation_methods = NULL,
                                       estimation_methods = NULL) {
      if (is.null(transformation_methods)) {
        transformation_methods <- c("box_cox", "yeo_johnson")
      }

      if (is.null(estimation_methods)) {
        estimation_methods <- setdiff(power.transform:::..estimators_all(), power.transform:::..estimators_raymaekers_robust())
      }

      lambda_list <- list()
      ii <- 1L
      for (transformation_method in transformation_methods) {
        for (estimation_method in estimation_methods) {
          # Create transformer object.
          transformer <- suppressWarnings(
            power.transform::find_transformation_parameters(
              x = x,
              method = transformation_method,
              estimation_method = estimation_method,
              robust = FALSE,
              shift = TRUE
            )
          )

          lambda_list[[ii]] <- data.table::data.table(
            "method" = transformation_method,
            "estimation_method" = estimation_method,
            "lambda" = transformer@lambda
          )

          ii <- ii + 1L
        }
      }

      return(data.table::rbindlist(lambda_list))
    }


    # Helper function for population distributions with outliers.
    .populate_outliers <- function(x,
                                   n,
                                   alpha,
                                   beta,
                                   target_lambda,
                                   ii,
                                   outlier_fraction,
                                   side) {
      # Set parameters.
      parameters <- list(
        "n" = n,
        "alpha" = alpha,
        "beta" = beta,
        "target_lambda" = target_lambda,
        "ii" = ii
      )

      # Add outliers.
      x <- lapply(
        outlier_fraction,
        function(k, x, parameters, side) {
          # Set parameters.
          parameters$k <- k

          # Compute interquartile range.
          interquartile_range <- stats::IQR(x)

          # Compute upper and lower quartiles.
          q_lower <- stats::quantile(x, probs = 0.25, names = FALSE)
          q_upper <- stats::quantile(x, probs = 0.75, names = FALSE)

          # Set data where the outliers will be copied into.
          x_outlier <- x

          if (k != 0.0) {
            n_draw <- ceiling(k * length(x))

            # Generate outlier values that are smaller than Q1 - 1.5 IQR or larger
            # than Q3 + 1.5 IQR.
            if (side == "both") {
              x_random <- stats::runif(n_draw, min = -2.0, max = 2.0)
            } else if (side == "right") {
              x_random <- stats::runif(n_draw, min = 0.0, max = 2.0)
            } else if (side == "left") {
              x_random <- stats::runif(n_draw, min = -2.0, max = 0.0)
            } else {
              stop(paste0("side was not recognised: ", side))
            }

            outlier <- numeric(n_draw)
            if (any(x_random < 0)) {
              outlier[x_random < 0] <- q_lower - 1.5 * interquartile_range + x_random[x_random < 0] * interquartile_range
            }

            if (any(x_random >= 0)) {
              outlier[x_random >= 0] <- q_upper + 1.5 * interquartile_range + x_random[x_random >= 0] * interquartile_range
            }

            # Randomly insert outlier values.
            x_outlier[sample(seq_along(x), size = n_draw, replace = FALSE)] <- outlier
          }

          return(list(
            "x" = x_outlier,
            "parameters" = parameters
          ))
        },
        x = x,
        parameters = parameters,
        side = side
      )

      return(x)
    }



    # Helper function for computing transformer parameters under optimisation
    # constraints.
    .compute_robust_lambda <- function(experiment,
                                       data) {
      # Custom parser.
      ..outlier_parser <- function(x,
                                   fun_args) {
        # Create transformer.
        transformer <- suppressWarnings(do.call(
          power.transform::find_transformation_parameters,
          args = c(
            list("x" = x$x),
            fun_args
          )
        ))

        # Prevent warnings.
        method <- estimation_method <- NULL

        target_lambda <- x$parameters$target_lambda[method == fun_args$method & estimation_method == fun_args$estimation_method]$lambda

        parameter_data <- data.table::data.table(
          "method" = fun_args$method,
          "estimation_method" = fun_args$estimation_method,
          "weight_method" = fun_args$weighting_function,
          "lambda" = transformer@lambda,
          "target_lambda" = target_lambda,
          "shift" = transformer@shift,
          "k" = x$parameters$k,
          "ii" = x$parameters$ii
        )

        return(parameter_data)
      }

      parameter_data <- lapply(
        data,
        ..outlier_parser,
        fun_args = experiment$fun_args
      )

      parameter_data <- data.table::rbindlist(
        parameter_data,
        fill = TRUE
      )

      return(parameter_data)
    }

    # computations -------------------------------------------------------------
    set.seed(95)

    # Generate table of asymmetric generalised normal distribution parameters.
    n_distributions <- 100L

    # Generate alpha, beta and n.
    n <- stats::runif(n = n_distributions, min = 2, max = 4)
    n <- ceiling(10^n)

    alpha <- stats::runif(n = n_distributions, min = 0.01, max = 0.99)
    beta <- stats::runif(n = n_distributions, min = 1.00, max = 5.00)

    # Generate corresponding distributions.
    x <- mapply(
      power.transform::ragn,
      n = n,
      alpha = alpha,
      beta = beta
    )

    # Compute lambda values without weighting.
    target_lambda <- lapply(x, .compute_target_lambda)

    # Add outliers, between 0.00 and 0.10.
    n_outliers <- 10
    outlier_fraction <- (seq_len(n_outliers) - 1)^2 / ((n_outliers - 1)^2 * 10)

    # Generate outliers in the data.
    x <- mapply(
      FUN = .populate_outliers,
      x = x,
      n = n,
      alpha = alpha,
      beta = beta,
      target_lambda = target_lambda,
      ii = seq_along(x),
      MoreArgs = list(
        "outlier_fraction" = outlier_fraction,
        "side" = side
      ),
      SIMPLIFY = FALSE,
      USE.NAMES = FALSE
    )

    x <- unlist(x, recursive = FALSE)

    experiments <- coro::collect(experiment_args(
      manuscript_dir = manuscript_dir,
      side = side
    ))

    cl <- parallel::makeCluster(18L)

    data <- parallel::parLapplyLB(
      cl = cl,
      X = experiments,
      fun = .compute_robust_lambda,
      data = x,
      chunk.size = 1L
    )

    # Stop cluster.
    parallel::stopCluster(cl)

    # Merge to data.table.
    data <- data.table::rbindlist(
      data,
      fill = TRUE
    )

    saveRDS(data, file = file.path(manuscript_dir, "robustness_comparison_marginal_parameter_data.RDS"))
  } else {
    data <- readRDS(data, file = file.path(manuscript_dir, "robustness_comparison_marginal_parameter_data.RDS"))
  }

  data[is.na(weight_method), "weight_method" := "none"]

  data$method <- factor(
    x = data$method,
    levels = c("box_cox", "yeo_johnson"),
    labels = c("Box-Cox", "Yeo-Johnson")
  )
  data$weight_method <- factor(
    x = data$weight_method,
    levels = c(
      "none",
      "empirical_probability_step", "empirical_probability_triangle", "empirical_probability_cosine",
      "transformed_step", "transformed_triangle", "transformed_cosine",
      "residual_step", "residual_triangle", "residual_cosine"
    ),
    labels = c(
      "non-robust",
      "emp. prob. (step)", "emp. prob. (triangle)", "emp. prob. (tap. cos.)",
      "z-score (step)", "z-score (triangle)", "z-score (tap. cos.)",
      "residual (step)", "residual (triangle)", "residual (tap. cos.)"
    )
  )
  data$estimation_method <- factor(
    x = data$estimation_method,
    levels = c("mle", "anderson_darling", "cramer_von_mises", "jarque_bera", "dagostino"),
    labels = c("MLE", "Anderson-Darling", "Cramér-von Mises", "Jarque-Bera", "D'Agostino")
  )

  return(data)
}



.get_goodness_of_fit_data <- function(
    manuscript_dir,
    residual_fun = NULL,
    n_distributions = 10000L,
    parallel = TRUE) {
  # Non-standard evaluation in data.table.
  p <- outlier_id <- NULL

  external_fun <- FALSE
  if (!is.null(residual_fun)) external_fun <- TRUE

  # Helper function for population distributions with outliers.
  .populate_outliers <- function(
      x,
      n,
      alpha,
      beta,
      ii,
      outlier_fraction,
      side = "both") {
    # Set parameters.
    parameters <- list(
      "n" = n,
      "alpha" = alpha,
      "beta" = beta,
      "ii" = ii
    )

    # Add outliers.
    x <- mapply(
      function(k, jj, x, parameters, side) {
        # Set parameters.
        parameters$k <- k
        parameters$jj <- jj

        # Compute interquartile range.
        interquartile_range <- stats::IQR(x)

        # Compute upper and lower quartiles.
        q_lower <- stats::quantile(x, probs = 0.25, names = FALSE)
        q_upper <- stats::quantile(x, probs = 0.75, names = FALSE)

        # Set data where the outliers will be copied into.
        x_outlier <- x

        if (k != 0.0) {
          n_draw <- ceiling(k * length(x))

          # Generate outlier values that are smaller than Q1 - 1.5 IQR or larger
          # than Q3 + 1.5 IQR.
          if (side == "both") {
            x_random <- stats::runif(n_draw, min = -2.0, max = 2.0)
          } else if (side == "right") {
            x_random <- stats::runif(n_draw, min = 0.0, max = 2.0)
          } else if (side == "left") {
            x_random <- stats::runif(n_draw, min = -2.0, max = 0.0)
          } else {
            stop(paste0("side was not recognised: ", side))
          }

          outlier <- numeric(n_draw)
          if (any(x_random < 0)) {
            outlier[x_random < 0] <- q_lower - 1.5 * interquartile_range + x_random[x_random < 0] * interquartile_range
          }

          if (any(x_random >= 0)) {
            outlier[x_random >= 0] <- q_upper + 1.5 * interquartile_range + x_random[x_random >= 0] * interquartile_range
          }

          # Randomly insert outlier values.
          x_outlier[sample(seq_along(x), size = n_draw, replace = FALSE)] <- outlier
        }

        return(list(
          "x" = x_outlier,
          "parameters" = parameters
        ))
      },
      k = outlier_fraction,
      jj = seq_along(outlier_fraction),
      MoreArgs = list(
        "x" = x,
        "parameters" = parameters,
        "side" = side
      ),
      SIMPLIFY = FALSE,
      USE.NAMES = FALSE
    )

    return(x)
  }


  .compute_residuals <- function(x, n_sample = 200L) {
    # Set interpolation points.
    p_sample <- (seq_len(n_sample) - 1 / 3) / (n_sample + 1 / 3)

    # Determine residuals for Box-Cox transformations.
    transformer_bc <- suppressWarnings(power.transform::find_transformation_parameters(
      x = x$x,
      method = "box_cox",
      shift = TRUE,
      robust = TRUE,
      estimation_method = "mle",
      weighting_function = "empirical_probability_cosine"
    ))

    residual_data <- power.transform::get_residuals(
      x = x$x,
      transformer = transformer_bc
    )

    residual <- stats::approx(
      x = residual_data$p,
      y = abs(residual_data$residual),
      xout = p_sample,
      rule = 2,
      ties = max
    )$y

    residual_bc <- data.table::data.table(
      "distribution_id" = x$parameter$ii,
      "outlier_id" = x$parameters$jj,
      "p" = seq_along(p_sample),
      "residual" = residual,
      "method" = "box_cox"
    )

    # Determine residuals for Yeo-Johnson transformations.
    transformer_yj <- suppressWarnings(power.transform::find_transformation_parameters(
      x = x$x,
      method = "yeo_johnson",
      shift = TRUE,
      robust = TRUE,
      estimation_method = "mle",
      weighting_function = "empirical_probability_cosine"
    ))

    residual_data <- power.transform::get_residuals(
      x = x$x,
      transformer = transformer_yj
    )

    residual <- stats::approx(
      x = residual_data$p,
      y = abs(residual_data$residual),
      xout = p_sample,
      rule = 2,
      ties = max
    )$y

    residual_yj <- data.table::data.table(
      "distribution_id" = x$parameter$ii,
      "outlier_id" = x$parameters$jj,
      "p" = seq_along(p_sample),
      "residual" = residual,
      "method" = "yeo_johnson"
    )

    data <- data.table::rbindlist(list(residual_bc, residual_yj))

    data$method <- factor(
      x = data$method,
      levels = c("box_cox", "yeo_johnson"),
      labels = c("Box-Cox", "Yeo-Johnson")
    )

    return(data)
  }


  if (!file.exists(file.path(manuscript_dir, "residual_plot.RDS")) ||
    external_fun) {
    if (!external_fun) {
      residual_fun <- .compute_residuals
    }

    # computations -------------------------------------------------------------
    set.seed(95)

    # Generate alpha, beta and n. Unlike for determining weighing function
    # parameters, we assess smaller datasets, notably between 30 and 1000
    # samples.
    n <- stats::runif(n = n_distributions, min = 1.47, max = 3)
    n <- ceiling(10^n)

    alpha <- stats::runif(n = n_distributions, min = 0.01, max = 0.99)
    beta <- stats::runif(n = n_distributions, min = 1.00, max = 5.00)

    # Generate corresponding distributions.
    x <- mapply(
      power.transform::ragn,
      n = n,
      alpha = alpha,
      beta = beta
    )

    # Add outliers, between 0.00 and 0.10.
    n_outliers <- 10
    outlier_fraction <- (seq_len(n_outliers) - 1)^2 / ((n_outliers - 1)^2 * 10)

    # Generate outliers in the data.
    x <- mapply(
      FUN = .populate_outliers,
      x = x,
      n = n,
      alpha = alpha,
      beta = beta,
      ii = seq_along(x),
      MoreArgs = list(
        "outlier_fraction" = outlier_fraction
      ),
      SIMPLIFY = FALSE,
      USE.NAMES = FALSE
    )

    set.seed(95)

    x <- unlist(x, recursive = FALSE)

    if (parallel) {
      # Start cluster
      cl <- parallel::makeCluster(18L)

      # Compute all data in parallel.
      data <- parallel::parLapply(
        cl = cl,
        X = x,
        fun = residual_fun
      )

      # Stop cluster.
      parallel::stopCluster(cl)
    } else {
      data <- lapply(
        X = x,
        FUN = residual_fun
      )
    }

    data <- data.table::rbindlist(data)

    if (!external_fun) {
      saveRDS(
        object = data,
        file = file.path(manuscript_dir, "residual_plot.RDS")
      )
    }
  } else {
    data <- readRDS(file.path(manuscript_dir, "residual_plot.RDS"))
  }

  # Convert p to empirical probabilities.
  data[, "p" := (p - 1 / 3) / (200 + 1 / 3)]
  data[, "has_outliers" := outlier_id > 1L]

  return(data)
}



.get_test_statistics_data <- function(manuscript_dir, as_table = TRUE, reduce = FALSE) {
  # Non-standard evaluation in data.table.
  residual <- p <- n <- test_statistic <- method <- alpha <- NULL

  ..down_sample_data <- function(x, alpha, upper, lower, n) {
    # Find the x-values.
    alpha_out <- mapply(
      function(upper, lower, n) {
        # Linear placement of interpolation points. The starting point (upper)
        # is not included.
        x <- rev(seq(from = lower, to = upper, length.out = n + 1L))
        return(tail(x, n = n))
      },
      upper = upper,
      lower = lower,
      n = n,
      SIMPLIFY = FALSE,
      USE.NAMES = FALSE
    )

    alpha_out <- unlist(alpha_out, use.names = FALSE)

    # Find the corresponding
    x_out <- stats::spline(
      x = alpha,
      y = x,
      method = "hyman",
      xout = alpha_out
    )$y

    # Add start and end points.
    x_out <- c(min(x), x_out, max(x))
    alpha_out <- c(1.0, alpha_out, 0.0)

    return(list("alpha" = alpha_out, "test_statistic" = x_out))
  }

  # Compute test statistic values.
  data <- .get_goodness_of_fit_data(manuscript_dir = manuscript_dir)

  central_width <- 0.80

  significance_values <- c(0.50, 0.20, 0.10, 0.05, 0.02, 0.01, 0.001)

  p_lower <- 0.50 - central_width / 2
  p_upper <- 0.50 + central_width / 2

  # Clip empirical probabilities to the center.
  data <- data[p >= p_lower & p <= p_upper]

  # Compute mean absolute residual error per feature.
  data <- data[, list("test_statistic" = mean(residual)), by = c("distribution_id", "outlier_id", "method")]
  data <- data[, list("n" = .N), by = c("test_statistic", "method")][order(test_statistic, method)]
  data[, "alpha" := 1.0 - cumsum(n) / sum(n), by = "method"]

  if (reduce) {
    # Reduce complexity of the data to save storage size.

    # 5 points between 1.00 and 0.95
    # 10 points between 0.95 and 0.75
    # 10 points between 0.75 and 0.25
    # 10 points between 0.25 and 0.10
    # 10 points between 0.10 and 0.05
    # 10 points between 0.05 and 0.01
    # 10 points between 0.01 and 0.005
    # 10 points between 0.005 and 0.001
    # 10 points between 0.001 and 0.0001
    upper <- c(1.00, 0.95, 0.75, 0.25, 0.10, 0.05, 0.010, 0.005, 0.0010)
    lower <- c(0.95, 0.75, 0.25, 0.10, 0.05, 0.01, 0.005, 0.001, 0.0001)
    n_int <- c(5, 10, 10, 10, 10, 10, 10, 10, 10)

    data <- data[, ..down_sample_data(
      x = test_statistic,
      alpha = alpha,
      upper = upper,
      lower = lower,
      n = n_int
    ),
    by = "method"
    ]
  } else {
    data <- rbind(
      data.table::data.table(
        "test_statistic" = c(0.0, 0.0),
        "method" = c("Box-Cox", "Yeo-Johnson"),
        "n" = 0L,
        "alpha" = 1.0
      ),
      data
    )
  }

  # Use Yeo-Johnson for test statistics.
  data <- data[method == "Yeo-Johnson"]

  if (as_table) {
    # Find thresholds.
    data <- data[, list(
      "test_statistic" = stats::spline(
        x = alpha,
        y = test_statistic,
        method = "hyman",
        xout = significance_values
      )$y,
      "alpha" = significance_values
    ),
    by = "method"
    ]
  }

  return(data)
}


.get_test_statistics_data_appendix <- function(manuscript_dir) {
  # Prevent warnings due to non-standard evaluation.
  residual <- p <- n <- mare <- method <- method_tag <- NULL

  .compute_residuals <- function(x) {
    ..compute_residuals <- function(x,
                                    method,
                                    shift,
                                    robust,
                                    estimation_method,
                                    weighting_function,
                                    method_tag,
                                    n_sample = 200L) {
      # Set interpolation points.
      p_sample <- (seq_len(n_sample) - 1 / 3) / (n_sample + 1 / 3)

      transformer <- suppressWarnings(power.transform::find_transformation_parameters(
        x = x$x,
        method = method,
        shift = shift,
        robust = robust,
        estimation_method = estimation_method,
        weighting_function = weighting_function
      ))

      residual_data <- power.transform::get_residuals(
        x = x$x,
        transformer = transformer
      )

      residual <- stats::approx(
        x = residual_data$p,
        y = abs(residual_data$residual),
        xout = p_sample,
        rule = 2,
        ties = max
      )$y

      return(data.table::data.table(
        "distribution_id" = x$parameter$ii,
        "outlier_id" = x$parameters$jj,
        "p" = seq_along(p_sample),
        "residual" = residual,
        "method" = method,
        "method_tag" = method_tag
      ))
    }

    # Determine residuals for robust, shift-sensitive transformations.
    residual_bc_robust_shift_mle <- ..compute_residuals(
      x = x,
      method = "box_cox",
      shift = TRUE,
      robust = TRUE,
      estimation_method = "mle",
      weighting_function = "empirical_probability_cosine",
      method_tag = "robust_shift_mle"
    )

    residual_yj_robust_shift_mle <- ..compute_residuals(
      x = x,
      method = "yeo_johnson",
      shift = TRUE,
      robust = TRUE,
      estimation_method = "mle",
      weighting_function = "empirical_probability_cosine",
      method_tag = "robust_shift_mle"
    )

    # Determine residuals for normal transformations.
    if (x$parameters$jj == 1L) {
      residual_bc_mle <- ..compute_residuals(
        x = x,
        method = "box_cox",
        shift = FALSE,
        robust = FALSE,
        estimation_method = "mle",
        weighting_function = "none",
        method_tag = "conventional_mle"
      )

      residual_yj_mle <- ..compute_residuals(
        x = x,
        method = "yeo_johnson",
        shift = FALSE,
        robust = FALSE,
        estimation_method = "mle",
        weighting_function = "none",
        method_tag = "conventional_mle"
      )
    } else {
      residual_bc_mle <- residual_yj_mle <- NULL
    }

    # Determine residuals for shift-sensitive robust cramér-von Mises
    # transformations.
    residual_bc_robust_shift_cm <- ..compute_residuals(
      x = x,
      method = "box_cox",
      shift = TRUE,
      robust = TRUE,
      estimation_method = "cramer_von_mises",
      weighting_function = "empirical_probability_cosine",
      method_tag = "robust_shift_cm"
    )

    residual_yj_robust_shift_cm <- ..compute_residuals(
      x = x,
      method = "yeo_johnson",
      shift = TRUE,
      robust = TRUE,
      estimation_method = "cramer_von_mises",
      weighting_function = "empirical_probability_cosine",
      method_tag = "robust_shift_cm"
    )

    # Determine residuals for normal transformations.
    if (x$parameters$jj == 1L) {
      residual_bc_cm <- ..compute_residuals(
        x = x,
        method = "box_cox",
        shift = FALSE,
        robust = FALSE,
        estimation_method = "cramer_von_mises",
        weighting_function = "none",
        method_tag = "conventional_cm"
      )

      residual_yj_cm <- ..compute_residuals(
        x = x,
        method = "yeo_johnson",
        shift = FALSE,
        robust = FALSE,
        estimation_method = "cramer_von_mises",
        weighting_function = "none",
        method_tag = "conventional_cm"
      )
    } else {
      residual_bc_cm <- residual_yj_cm <- NULL
    }

    data <- data.table::rbindlist(list(
      residual_bc_robust_shift_mle,
      residual_yj_robust_shift_mle,
      residual_bc_mle,
      residual_yj_mle,
      residual_bc_robust_shift_cm,
      residual_yj_robust_shift_cm,
      residual_bc_cm,
      residual_yj_cm
    ))

    data$method <- factor(
      x = data$method,
      levels = c("box_cox", "yeo_johnson"),
      labels = c("Box-Cox", "Yeo-Johnson")
    )

    data$method_tag <- factor(
      x = data$method_tag,
      levels = c("robust_shift_mle", "conventional_mle", "robust_shift_cm", "conventional_cm"),
      labels = c("robust, shift sens. (MLE)", "conventional (MLE)", "robust, shift sens. (C-vM)", "conventional (C-vM)")
    )

    return(data)
  }

  # computations ---------------------------------------------------------------

  if (!file.exists(file.path(manuscript_dir, "residual_plot_appendix.RDS"))) {
    data <- .get_goodness_of_fit_data(
      manuscript_dir = manuscript_dir,
      residual_fun = .compute_residuals,
      n_distributions = 10000L,
      parallel = TRUE
    )

    central_width <- c(0.60, 0.70, 0.80, 0.90, 0.95, 1.00)

    mare_data <- list()

    for (ii in seq_along(central_width)) {
      p_lower <- 0.50 - central_width[ii] / 2
      p_upper <- 0.50 + central_width[ii] / 2

      x <- data.table::copy(data)

      # Clip empirical probabilities to centre.
      x <- x[p >= p_lower & p <= p_upper]

      # Compute mean absolute residual error per feature.
      x <- x[, list("mare" = mean(residual)), by = c("distribution_id", "outlier_id", "method", "method_tag")]
      x <- x[, list("n" = .N), by = c("mare", "method", "method_tag")][order(mare, method, method_tag)]
      x[, "rejected" := 1.0 - cumsum(n) / sum(n), by = c("method", "method_tag")]

      x[, "kappa" := central_width[ii]]

      mare_data[[ii]] <- x
    }

    data <- data.table::rbindlist(mare_data)

    saveRDS(
      object = data,
      file = file.path(manuscript_dir, "residual_plot_appendix.RDS")
    )
  } else {
    data <- readRDS(file.path(manuscript_dir, "residual_plot_appendix.RDS"))
  }

  return(data)
}



.get_ml_experiment_data <- function(manuscript_dir) {
  # Converts experiment results into something that can be interpreted.
  if (file.exists(file.path(manuscript_dir, "parsed_ml_results.RDS"))) {
    return(readRDS(file.path(manuscript_dir, "parsed_ml_results.RDS")))
  }

  if (!file.exists(file.path(manuscript_dir, "ml_exp_results_offset.RDS"))) {
    stop("Execute ml_exp.R to create the results and copy it to the manuscript folder.")
  }

  # non-standard evaluation
  experiment_parameters <- dataset_split <- value <- metric <- NULL

  results <- readRDS(file.path(manuscript_dir, "ml_exp_results_offset.RDS"))

  results[, "experiment_parameters" := sub(pattern = "_glm", replacement = "", x = experiment_parameters)]
  results[, "experiment_parameters" := sub(pattern = "_rf", replacement = "", x = experiment_parameters)]

  results <- results[dataset_split == "test"]

  # Set data difficulty.
  results[, "data_difficulty" := familiar.experiment::assess_difficulty(
    x = stats::median(value),
    metric = metric),
    by = c("dataset")
  ]

  # Convert dataset to category.
  results$dataset <- factor(
    x = results$dataset,
    levels = sort(unique(results$dataset))
  )

  # Rename experiment_parameters and convert to category.
  results$experiment_parameters <- factor(
    x = results$experiment_parameters,
    levels = c("config_no_transformation", "config_old_test", "config_new_test"),
    labels = c("none", "conventional", "proposed")
  )

  # Convert to ranks. Higher ranks are better results. Values are mapped to [0, 1].
  results[metric == "root_relative_squared_error_winsor", "value_rank" := (
    data.table::frank(-value, ties.method = "average") - 1.0) / (.N - 1),
    by = c("dataset")
  ]
  results[metric == "auc_roc", "value_rank" := (
    data.table::frank(value, ties.method = "average") - 1.0) / (.N - 1),
    by = c("dataset")
  ]

  # Map to [-1.0, 1.0] so that values are centered on 0.0.
  results[, "value_rank" := (value_rank - 0.5) * 2]

  data <- list(
    "n" = nrow(results),
    "n_dataset" = nlevels(results$dataset),
    "n_transformer" = nlevels(results$experiment_parameters),
    "n_learner" = nlevels(results$learner),
    "n_iteration" = data.table::uniqueN(results$iteration_id),
    "id_transformer" = as.numeric(results$experiment_parameters),
    "id_dataset" = as.numeric(results$dataset),
    "id_learner" = as.numeric(results$learner),
    "id_iteration" = results$iteration_id,
    "y" = results$value_rank
  )

  # TODO: Loop over problem difficulty. It may well be that easier problems do
  # not benefit as much from power transformations.

  model_data <- rstan::stan(
    file = file.path(manuscript_dir, "model.stan"),
    data = data,
    cores = 4,
    open_progress = FALSE
  )

  # results <- data.table::dcast(
  #   data = results,
  #   dataset + iteration_id + learner + dataset_split + metric + outcome_type ~ experiment_parameters,
  #   value.var = "value")
  results_test <- data.table::dcast(
    data = results,
    dataset + iteration_id + learner + data_difficulty ~ experiment_parameters,
    value.var = "value_rank"
  )


}
