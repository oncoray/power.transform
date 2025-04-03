.get_data_problematic_transformations <- function(manuscript_dir) {
  # Set seed.
  set.seed(19L)

  file_name <- file.path(manuscript_dir, "problematic_transformations_plot.RDS")

  if (!file.exists(file_name)) {
    x <- power.transform::ragn(
      10000L,
      location = 0.0,
      scale = 1.0 / sqrt(2.0),
      alpha = 0.5,
      beta = 2.0
    )

    # generator ----------------------------------------------------------------
    generate_experiment_data <- coro::generator(
      function(x) {
        shift_range <- 10^seq(from = 0, to = 6, by = 0.1)
        scale_range <- 10^seq(from = -6, to = 0, by = 0.1)
        outlier_range <- 10^seq(from = -1, to = 6, by = 0.1)

        for (method in c("box_cox", "yeo_johnson")) {
          # Iterate over shift range.
          for (d in shift_range) {
            coro::yield(list(
              "x" = x + d,
              "y" = d,
              "method" = method,
              "invariant" = FALSE,
              "robust" = FALSE,
              "estimation_method" = "mle",
              "data_type" = "shift"
            ))
          }

          # Iterate over scale range.
          for (s in scale_range) {
            coro::yield(list(
              "x" = x * s + 10.0,
              "y" = s,
              "method" = method,
              "invariant" = FALSE,
              "robust" = FALSE,
              "estimation_method" = "mle",
              "data_type" = "scale"
            ))
          }


          for (d in outlier_range) {
            coro::yield(list(
              "x" = c(x, d),
              "y" = d,
              "method" = method,
              "invariant" = FALSE,
              "robust" = FALSE,
              "estimation_method" = "mle",
              "data_type" = "outlier"
            ))
          }
        }
      }
    )

    .compute_lambda <- function(parameter_set) {

      transformer <- suppressWarnings(
        power.transform::find_transformation_parameters(
          x = parameter_set$x,
          method = parameter_set$method,
          robust = parameter_set$robust,
          invariant = parameter_set$invariant,
          estimation_method = parameter_set$estimation_method
        )
      )

      return(
        data.table::data.table(
          "method" = parameter_set$method,
          "y" = log10(parameter_set$y),
          "lambda" = transformer@lambda,
          "x_0" = transformer@shift,
          "s" = transformer@scale,
          "data_type" = parameter_set$data_type
        )
      )
    }

    experiments <- coro::collect(generate_experiment_data(x = x))

    data <- lapply(
      X = experiments,
      FUN = .compute_lambda
    )

    # Start cluster
    cl <- parallel::makeCluster(8L)

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

  # Update method and data type..
  data$method <- factor(
    x = data$method,
    levels = c("box_cox", "yeo_johnson"),
    labels = c("Box-Cox", "Yeo-Johnson")
  )
  data$data_type <- factor(
    x = data$data_type,
    levels = c("shift", "scale", "outlier")
  )

  return(data)
}



.get_shifted_scaled_distribution_data <- function(manuscript_dir, main_manuscript) {
  # Set seed.
  set.seed(19L)

  if (main_manuscript) {
    file_name <- file.path(manuscript_dir, "main_shifted_scaled_distributions_plot.RDS")
    estimation_methods <- "mle"
  } else {
    file_name <- file.path(manuscript_dir, "appendix_shifted_scaled_distributions_plot.RDS")
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
                    coro::yield(list(
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
                    coro::yield(list(
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

  if (!file.exists(file.path(manuscript_dir, "main_shifted_outlier_plot.RDS"))) {
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



.get_optimised_weighting_function_parameters <- function(
    manuscript_dir
) {

  # generator ------------------------------------------------------------------
  experiment_args <- coro::generator(
    function(
        manuscript_dir,
        transformation_method = NULL,
        estimation_methods = NULL
    ) {

      file_dir <- "robustness_parameters"

      if (is.null(estimation_methods)) {
        estimation_methods <- setdiff(
          power.transform:::..estimators_all(),
          power.transform:::..estimators_raymaekers_robust())
        # estimation_methods <- "mle"
      }

      for (estimation_method in estimation_methods) {

        # Non-robust transformation.
        fun_args <- list(
          "method" = transformation_method,
          "estimation_method" = estimation_method,
          "robust" = FALSE,
          "invariant" = TRUE
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

          for (robustness_weighting_function in c("step", "triangle", "cosine")) {
            fun_args <- list(
              "method" = transformation_method,
              "estimation_method" = estimation_method,
              "robust" = TRUE,
              "invariant" = TRUE,
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
  )


  # Helper function for population distributions with outliers.
  .populate_outliers <- function(
    x,
    n,
    alpha,
    beta,
    outlier_fraction,
    side
  ) {
    # Set parameters.
    parameters <- list(
      "n" = n,
      "alpha" = alpha,
      "beta" = beta
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
  .compute_robust_error <- function(
    x,
    fun_args,
    opt_args
  ) {
    # Custom parser.
    ..outlier_parser <- function(
      x,
      fun_args,
      opt_args
    ) {
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

      # Transform data
      transformed_data <- tryCatch(
        power.transform::power_transform(
          x = x$x,
          transformer = transformer
        ),
        error = identity
      )

      if (inherits(transformed_data, "simpleError")) {
        residual_error <- NA_real_
      } else {
        residual_data <- power.transform::get_residuals(
          x = x$x,
          transformer = transformer
        )

        residual_error <- mean(abs(residual_data[p_observed >= 0.1 & p_observed <= 0.90]$residual))
      }

      parameter_data <- c(
        opt_args,
        list(
          "method" = fun_args$method,
          "estimation_method" = fun_args$estimation_method,
          "weight_method" = fun_args$weighting_function,
          "lambda" = transformer@lambda,
          "shift" = transformer@shift,
          "scale" = transformer@scale,
          "residual_error" = residual_error
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

    parameter_data <- .compute_robust_error(
      x = x,
      fun_args = fun_args,
      opt_args = opt_args
    )

    errors <- sapply(parameter_data, function(x) (x$residual_error))

    if (verbose) {
      message <- paste0("target: ", sum(errors))
      if (length(opt_args) > 0) message <- c(message, paste0("k1: ", opt_param[1]))
      if (length(opt_args) > 1) message <- c(message, paste0("k2: ", opt_param[2]))
      cat(paste0(message, collapse = "; "))
      cat("\n")
    }

    return(sum(errors))
  }

  # computations ---------------------------------------------------------------

  # Select experiments to conduct by filtering out experiments that are file
  # names.
  experiments <- c(
    coro::collect(experiment_args(
      manuscript_dir = manuscript_dir,
      transformation_method = "yeo_johnson"
    )),
    coro::collect(experiment_args(
      manuscript_dir = manuscript_dir,
      transformation_method = "box_cox"
    ))
  )
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

    # Add outliers, between 0.00 and 0.10.
    outlier_fraction <- 0.10

    # Generate outliers in the data.
    x <- mapply(
      FUN = .populate_outliers,
      x = x,
      n = n,
      alpha = alpha,
      beta = beta,
      MoreArgs = list(
        "outlier_fraction" = outlier_fraction,
        "side" = "both"
      ),
      SIMPLIFY = FALSE,
      USE.NAMES = FALSE
    )

    x <- unlist(x, recursive = FALSE)

    # lapply(experiments, .optimisation_outer, data=x, verbose=TRUE)

    # Start cluster
    cl <- parallel::makeCluster(18L)
    on.exit(parallel::stopCluster(cl))

    parallel::clusterExport(
      cl = cl,
      varlist = c(".compute_robust_error", ".optimisation_inner"),
      envir = environment()
    )

    parallel::parLapplyLB(
      cl = cl,
      X = experiments,
      fun = .optimisation_outer,
      data = x,
      chunk.size = 1L
    )
  }

  # Collect all files
  experiments <- c(
    coro::collect(experiment_args(
      manuscript_dir = manuscript_dir,
      transformation_method = "yeo_johnson"
    )),
    coro::collect(experiment_args(
      manuscript_dir = manuscript_dir,
      transformation_method = "box_cox"
    ))
  )
  weighting_parameters <- lapply(experiments, readRDS)
  weighting_parameters <- lapply(weighting_parameters, data.table::as.data.table)
  weighting_parameters <- data.table::rbindlist(weighting_parameters, fill = TRUE)

  return(weighting_parameters)
}



.get_ml_experiment_data <- function(manuscript_dir, subset = "numeric") {
  # non-standard evaluation
  experiment_parameters <- dataset_split <- value <- value_rank <- metric <- NULL

  data <- readRDS(file.path(manuscript_dir, "ml_experiment", "results", "performance.RDS"))

  # Strip details from experiment_parameters
  data[, "experiment_parameters" := sub(pattern = "_glm", replacement = "", x = experiment_parameters)]
  data[, "experiment_parameters" := sub(pattern = "_rf", replacement = "", x = experiment_parameters)]
  data[, "experiment_parameters" := sub(pattern = "_xgboost", replacement = "", x = experiment_parameters)]
  data[, "experiment_parameters" := sub(pattern = "_lasso", replacement = "", x = experiment_parameters)]
  data[, "experiment_parameters" := sub(pattern = "config_", replacement = "", x = experiment_parameters)]

  # Detail transformation method.
  data[, "transformation_method" := "invariant_robust"]
  data[grepl("no_transformation", experiment_parameters, fixed = TRUE), "transformation_method" := "none"]
  data[grepl("conventional", experiment_parameters, fixed = TRUE), "transformation_method" := "conventional"]
  data[grepl("invariant_robust_gof", experiment_parameters, fixed = TRUE), "transformation_method" := "invariant_robust_gof"]
  data$transformation_method <- factor(
    data$transformation_method,
    levels = c("none", "conventional", "invariant_robust", "invariant_robust_gof")
  )

  # Detail normalisation method.
  data[, "normalisation_method" := "robust_standardisation"]
  data[grepl("no_normalisation", experiment_parameters, fixed = TRUE), "normalisation_method" := "none"]
  data$normalisation_method <- factor(
    data$normalisation_method,
    levels = c("none", "robust_standardisation")
  )

  # Select only test data.
  data <- data[dataset_split == "test"]

  # Set data difficulty.
  data[, "task_difficulty" := familiar.experiment::assess_difficulty(
    x = stats::median(value),
    metric = metric),
    by = c("dataset")
  ]

  # Convert to ranks. Higher ranks are better results. Values are mapped to [0,
  # 1].
  data[
    metric == "root_relative_squared_error_winsor",
    "value_rank" := data.table::frank(-value, ties.method = "average"),
    by = c("dataset")
  ]
  data[
    metric == "auc_roc",
    "value_rank" := data.table::frank(value, ties.method = "average"),
    by = c("dataset")
  ]
  data[, "value_rank" := (value_rank - min(value_rank)) / (max(value_rank) - min(value_rank)), by = c("dataset") ]

  # Select subset of data.
  # if (subset == "numeric") {
  #   feature_data <- .get_ml_experiment_feature_statistics(manuscript_dir = manuscript_dir)
  #   feature_data[, "dataset_real" := gsub(pattern = "_data.RDS", replacement = "", x = dataset, fixed = TRUE)]
  #   data <- data[dataset %in% feature_data$dataset_real]
  #
  # } else if (subset == "numeric_high_shift") {
  #   feature_data <- .get_ml_experiment_feature_statistics(manuscript_dir = manuscript_dir)
  #   feature_data <- feature_data[abs(mu) > 500.0]
  #   feature_data[, "dataset_real" := gsub(pattern = "_data.RDS", replacement = "", x = dataset, fixed = TRUE)]
  #   data <- data[dataset %in% feature_data$dataset_real]
  # }

  return(data)
}



.get_ml_experiment_feature_statistics <- function(manuscript_dir) {

  file_name <- file.path(manuscript_dir, "ml_experiment_feature_statistics.RDS")

  if (!file.exists(file_name)) {
    data_dir <- file.path(manuscript_dir, "ml_experiment", "data_cache")
    data_files <- list.files(path = data_dir, full.names = FALSE, pattern = ".RDS")

    feature_data <- lapply(
      data_files,
      function(data_file, data_dir) {
        # Read datasets, which are familiar DataObject objects.
        dataset <- readRDS(file.path(data_dir, data_file))
        feature_data <- lapply(
          familiar:::get_feature_columns(dataset),
          function(feature, dataset) {
            if (is.numeric(dataset@data[[feature]])) {
              x <- power.transform::huber_estimate(dataset@data[[feature]])
              return(data.table::data.table("feature" = feature, "mu" = x$mu, "sigma" = x$sigma))

            } else {
              return(NULL)
            }
          },
          dataset = dataset
        )
        feature_data <- data.table::rbindlist(feature_data)
        if (nrow(feature_data) > 0) {
          feature_data[, "dataset" := data_file]
        }
        return(feature_data)
      },
      data_dir = data_dir
    )

    feature_data <- data.table::rbindlist(feature_data)

    saveRDS(
      object = feature_data,
      file = file_name
    )
  } else {
    feature_data <- readRDS(file_name)
  }

  return(feature_data)
}



.assess_simulated_data <- function(
    manuscript_dir,
    method = "yeo_johnson",
    data_type = "clean",
    k = 0.10
) {
  file_name <- paste0(
    "simulated_data_",
    data_type,
    "_", method, ".RDS"
  )
  file_name <- file.path(manuscript_dir, file_name)

  set.seed(95L)

  if (file.exists(file_name)) return(readRDS(file_name))

  # Generate distributions to assess.
  n_distributions <- 10000L

  # Generate alpha, beta and n.
  n <- stats::runif(n = n_distributions, min = 1.0, max = 4.0)
  n <- ceiling(10^n)

  # Parameters for asymmetric generalised normal distributions.
  alpha <- stats::runif(n = n_distributions, min = 0.01, max = 0.99)
  beta <- stats::runif(n = n_distributions, min = 1.00, max = 5.00)

  # Parameters for normal distributed data that is inversely transformed.
  # Note that Yeo-Johnson is only directly invertible for lambda between 0 and
  # 2, and Box-Cox only if lambda > 0.
  lambda <- stats::runif(n = n_distributions, min = 0.00, max = 2.00)

  ..compute <- function(
    ii,
    n,
    alpha,
    beta,
    conf_lambda,
    data_type,
    method,
    k
  ) {
    set.seed(ii)
    if (data_type == "clean") {
      # Clean data is produced by first inverting the power transformation (with
      # known lambda) for normally distributed data. This means that the
      # transformation parameters can be exactly determined - aside from sampling
      # effects.

      x <- stats::rnorm(
        n = n,
        mean = 0.0,
        sd = 1.0
      )

      if (method == "box_cox" && any(x <= -1.0 / conf_lambda)) x <- x + 1E-4 - min(x) - 1.0 / conf_lambda

      # Perform inverse transform.
      x <- power.transform::revert_power_transform(
        y = x,
        transformer = power.transform::create_transformer_skeleton(
          method = method,
          lambda = conf_lambda
        )
      )

    } else if (data_type == "dirty") {
      # Dirty data are produced from asymmetric generalised normal distributions.
      # These are dirty in the sense that the parameters of the power
      # transformation cannot be directly derived from the distribution.
      x <- power.transform::ragn(
        n,
        location = 0,
        scale = 1 / sqrt(2),
        alpha = alpha,
        beta = beta
      )

      if (method == "box_cox" && any(x <= 0.0)) x <- x + 1E-4 - min(x)

    } else if (data_type == "dirty_shifted") {
      # Dirty data that is shifted and scaled to be far from 0.0.
      x <- power.transform::ragn(
        n,
        location = 100.0,
        scale = 0.001 * 1 / sqrt(2),
        alpha = alpha,
        beta = beta
      )

      if (method == "box_cox" && any(x <= 0.0)) x <- x + 1E-4 - min(x)

    } else {
      stop(paste0("data_type was not recognised: ", data_type))
    }

    # Ensure that x can be transformed, and is entirely positive


    # Compute upper and lower quartiles, and IQR.
    q_lower <- stats::quantile(x, probs = 0.25, names = FALSE)
    q_upper <- stats::quantile(x, probs = 0.75, names = FALSE)
    interquartile_range <- stats::IQR(x)

    # Set data where the outliers will be copied into.
    x_outlier <- x

    # Generate outlier values that are smaller than Q1 - 1.5 IQR or larger
    # than Q3 + 1.5 IQR.
    n_draw <- ceiling(k * length(x))
    # if (method == "box_cox") {
    #   # Only positive outliers for Box-Cox
    #   x_random <- stats::runif(n_draw, min = 0.0, max = 2.0)
    # } else {
    x_random <- stats::runif(n_draw, min = -2.0, max = 2.0)
    # }

    outlier <- numeric(n_draw)
    if (any(x_random < 0)) {
      outlier[x_random < 0] <- q_lower - 1.5 * interquartile_range + x_random[x_random < 0] * interquartile_range
    }

    if (any(x_random >= 0)) {
      outlier[x_random >= 0] <- q_upper + 1.5 * interquartile_range + x_random[x_random >= 0] * interquartile_range
    }

    # Randomly insert outlier values.
    x_outlier[sample(seq_along(x), size = n_draw, replace = FALSE)] <- outlier

    data_wo_outlier <- suppressWarnings(
      ..standardisation_before_transformation(
        x = x,
        data_name = data_type,
        method = method
      )
    )
    data_wo_outlier[, "outlier" := FALSE]
    data_w_outlier <- suppressWarnings(
      ..standardisation_before_transformation(
        x = x_outlier,
        data_name = data_type,
        method = method
      )
    )
    data_w_outlier[, "outlier" := TRUE]

    data <- rbind(data_wo_outlier, data_w_outlier)

    # Add data concerning the experiment.
    data[, "experiment_id" := ii]
    data[, "n_samples" := n]
    if (data_type == "clean") {
      data[, "conf_lambda" := conf_lambda]
    } else {
      data[, ":="("conf_alpha" = alpha, "conf_beta" = beta)]
    }

    return(data)
  }

  cl <- parallel::makeCluster(18L)
  on.exit(parallel::stopCluster(cl))

  parallel::clusterExport(cl = cl, "..standardisation_before_transformation")

  data <- parallel::clusterMap(
    cl = cl,
    fun = ..compute,
    ii = seq_len(n_distributions),
    n = n,
    alpha = alpha,
    beta = beta,
    conf_lambda = lambda,
    MoreArgs = list(
      "data_type" = data_type,
      "method" = method,
      "k" = k
    ),
    SIMPLIFY = FALSE,
    USE.NAMES = FALSE
  )
  data <- data.table::rbindlist(data)

  # Cache data.
  saveRDS(data, file_name)

  return(data)
}



.assess_standardisation_before_transformation_simulation <- function(
    lambda_limit = NULL,
    method = "yeo_johnson"
) {
  data <- list()
  set.seed(1L)

  n <- c(30L, 100L, 500L)
  lambda <- c(0.1, 1.0, 1.9)
  dataset_names <- character(length(n) * length(lambda))

  for (ii in seq_along(n)) {
    for (jj in seq_along(lambda)) {
      current_iteration <- (ii - 1L) * length(lambda) + jj

      set.seed(current_iteration)

      # Generate base, clean data.
      x <- stats::rnorm(
        n = n[ii],
        mean = 0.0,
        sd = 1.0
      )

      # Perform inverse transform.
      transformer <- power.transform::create_transformer_skeleton(
        method = method,
        lambda = lambda[jj]
      )
      x <- power.transform::revert_power_transform(
        y = x,
        transformer = transformer
      )

      dataset_names[current_iteration] <- paste0("n: ", n[ii], "; $\\lambda$:", lambda[jj])

      data[[current_iteration]] <- ..standardisation_before_transformation(
        x = x,
        data_name = dataset_names[current_iteration],
        lambda_limit = lambda_limit,
        method = method
      )
    }
  }

  data <- rbindlist(data, use.names = TRUE)
  data$dataset <- factor(data$dataset, levels = dataset_names)

  return(data)
}



.assess_standardisation_before_transformation <- function(
    lambda_limit = NULL,
    method = "yeo_johnson"
) {
  data <- list()

  # Age of lung cancer patients
  data[[1L]] <- ..standardisation_before_transformation(
    x = survival::lung$age,
    data_name = "age",
    lambda_limit = lambda_limit,
    method = method
  )

  # Penguin body mass
  x <- modeldata::penguins$body_mass_g
  x <- x[is.finite(x)]
  data[[2L]] <- ..standardisation_before_transformation(
    x = x,
    data_name = "penguin body mass",
    lambda_limit = lambda_limit,
    method = method
  )

  # Housing latitude
  data[[3L]] <- ..standardisation_before_transformation(
    x = modeldata::ames$Latitude,
    data_name = "latitude",
    lambda_limit = lambda_limit,
    method = method
  )

  # Fuel efficiency
  utils::data("TopGear", package = "robustHD", envir = environment())
  tg_data <- data.table::as.data.table(TopGear)
  x <- tg_data$MPG
  x <- x[!is.na(x)]
  data[[4L]] <- ..standardisation_before_transformation(
    x = x,
    data_name = "fuel efficiency",
    lambda_limit = lambda_limit,
    method = method
  )

  # Maximum arterial wall thickness
  data[[5L]] <- ..standardisation_before_transformation(
    x = modeldata::ischemic_stroke$max_max_wall_thickness,
    data_name = "arterial wall thickness",
    lambda_limit = lambda_limit,
    method = method
  )

  data <- rbindlist(data, use.names = TRUE)
  return(data)
}


..standardisation_before_transformation <- function(
    x,
    data_name,
    lambda_limit = c(-4.0, 6.0),
    method = "yeo_johnson"
) {

  ..parse_data <- function(
    transformation_method,
    x
  ) {

    if (endsWith(transformation_method, "(z-score norm.)")) {
      x <- (x - mean(x)) / stats::sd(x)

    } else if (endsWith(transformation_method, "(robust scaling)")) {
      x <- (x - stats::median(x)) / stats::IQR(x)
    }

    if (startsWith(transformation_method, "none")) {
      # No transformer
      transformer <- power.transform::find_transformation_parameters(
        x = x,
        method = "none"
      )

    } else if (startsWith(transformation_method, "conventional")) {
      # Conventional transformer
      transformer <- power.transform::find_transformation_parameters(
        x = x,
        method = method,
        robust = FALSE,
        invariant = FALSE,
        lambda = lambda_limit
      )

    } else if (startsWith(transformation_method, "Raymaekers-Rousseeuw")) {
      # Raymaekers robust transformer
      transformer <- power.transform::find_transformation_parameters(
        x = x,
        method = method,
        robust = TRUE,
        invariant = FALSE,
        estimation_method = "raymaekers_robust",
        lambda = lambda_limit
      )

    } else if (startsWith(transformation_method, "invariant")) {
      # Location- and scale-invariant transformer
      transformer <- power.transform::find_transformation_parameters(
        x = x,
        method = method,
        robust = FALSE,
        invariant = TRUE,
        lambda = lambda_limit
      )

    } else if (startsWith(transformation_method, "robust invariant")) {
      # Robust location- and scale-invariant transformer
      transformer <- power.transform::find_transformation_parameters(
        x = x,
        method = method,
        robust = TRUE,
        invariant = TRUE,
        lambda = lambda_limit
      )

    } else {
      stop(paste0("No known transformer_method: ", transformation_method))
    }

    # Transform data
    transformed_data <- power.transform::power_transform(
      x = x,
      transformer = transformer
    )

    residual_data <- power.transform::get_residuals(
      x = x,
      transformer = transformer
    )

    summary_data <- data.table::data.table(
      "residuals" = sum(abs(residual_data$residual)),
      "residuals_central" = sum(abs(residual_data[p_observed >= 0.1 & p_observed <= 0.9]$residual)),
      "mu" = mean(transformed_data),
      "sigma" = sd(transformed_data),
      "lambda" = power.transform::get_lambda(transformer),
      "shift" = power.transform::get_shift(transformer),
      "scale" = power.transform::get_scale(transformer),
      "p_value" = power.transform::assess_transformation(
        x = x,
        transformer = transformer,
        verbose = FALSE
      ),
      "transformation" = transformation_method
    )

    return(summary_data)
  }

  transformer_labels <- c(
    "none",
    "conventional", "conventional (z-score norm.)", "conventional (robust scaling)",
    "Raymaekers-Rousseeuw", "Raymaekers-Rousseeuw (z-score norm.)", "Raymaekers-Rousseeuw (robust scaling)",
    "invariant", "robust invariant"
  )

  data <- lapply(
    transformer_labels,
    ..parse_data,
    x = x
  )

  data <- data.table::rbindlist(data, use.names = TRUE)
  data$transformation <- factor(x = data$transformation, levels = transformer_labels)
  data[, "dataset" := data_name]

  return(data)
}




.get_test_statistics_data <- function(manuscript_dir, with_outliers = FALSE) {
  n_rep <- 30000L  # number of distributions drawn.
  k <- 0.10  # outlier fraction
  kappa <- c(0.60, 0.70, 0.80, 0.90, 0.95, 1.00)  # central portion of distribution.

  # From 5 to 10000 in equal steps (log 10 scale)
  n <- unique(floor(10^(seq(from = 0.75, to = 4.0, by = 0.0625))))

  if (with_outliers) {
    file_name <- file.path(manuscript_dir, "raw_ecn_statistics_w_outlier.RDS")
  } else {
    file_name <- file.path(manuscript_dir, "raw_ecn_statistics.RDS")
  }

  if (!file.exists(file_name)) {

    cl <- parallel::makeCluster(18L)
    on.exit(parallel::stopCluster(cl))

    data <- parallel::parLapplyLB(
      cl = cl,
      X = seq_along(n),
      fun = ..get_test_statistics_data,
      n = n,
      n_rep = n_rep,
      k = k,
      kappa = kappa,
      with_outliers = with_outliers
    )

    # Aggregate to single table.
    data <- data.table::rbindlist(data)
    data <- data[order(n, kappa, mean_residual_error)]

    # Store
    saveRDS(data, file_name)

  } else {
    data <- readRDS(file_name)
  }

  return(data)
}


..get_test_statistics_data <- function(ii, n, n_rep, k, kappa, with_outliers = FALSE) {
  set.seed(19L + ii)

  data <- list()
  # Repeat 1000 times.
  for (jj in seq_len(n_rep)) {
    x <- power.transform::ragn(floor(n[ii]), location = 0, scale = 1 / sqrt(2), alpha = 0.5, beta = 2)

    # Compute upper and lower quartiles, and IQR.
    q_lower <- stats::quantile(x, probs = 0.25, names = FALSE)
    q_upper <- stats::quantile(x, probs = 0.75, names = FALSE)
    interquartile_range <- stats::IQR(x)

    if (with_outliers) {
      # Set data where the outliers will be copied into.
      x_outlier <- x

      # Generate outlier values that are smaller than Q1 - 1.5 IQR or larger
      # than Q3 + 1.5 IQR.
      n_draw <- ceiling(k * length(x))
      x_random <- stats::runif(n_draw, min = -2.0, max = 2.0)

      outlier <- numeric(n_draw)
      if (any(x_random < 0)) {
        outlier[x_random < 0] <- q_lower - 1.5 * interquartile_range + x_random[x_random < 0] * interquartile_range
      }

      if (any(x_random >= 0)) {
        outlier[x_random >= 0] <- q_upper + 1.5 * interquartile_range + x_random[x_random >= 0] * interquartile_range
      }

      # Randomly insert outlier values.
      x_outlier[sample(seq_along(x), size = n_draw, replace = FALSE)] <- outlier
      x_outlier <- sort(x_outlier)

    } else {
      # Use data without outlier.
      x_outlier <- sort(x)
    }

    # Compute the expected z-score.
    z_expected <- power.transform:::compute_expected_z(x = x_outlier)

    # Compute M-estimates for locality and scale
    robust_estimates <- power.transform::huber_estimate(x_outlier, tol = 1E-3)

    # Avoid division by 0.0.
    if (robust_estimates$sigma == 0.0) robust_estimates$sigma <- 1.0

    # Compute the observed z-score.
    z_observed <- (x_outlier - robust_estimates$mu) / robust_estimates$sigma

    # Compute residuals.
    residual <- z_observed - z_expected
    residual_data <- data.table::data.table(
      "z_expected" = z_expected,
      "z_observed" = z_observed,
      "residual" = residual,
      "p_observed" = (seq_along(x) - 1 / 3) / (length(x) + 1 / 3)
    )

    # Compute mean residual error for all kappa.
    mean_residual_error <- numeric(length(kappa))
    for (kk in seq_along(kappa)) {
      p_lower <- 0.50 - kappa[kk] / 2
      p_upper <- 0.50 + kappa[kk] / 2

      mean_residual_error[kk] <- mean(abs(residual_data[p_observed >= p_lower & p_observed <= p_upper]$residual))
    }

    # Set data.
    data[[jj]] <- data.table::data.table(
      "n" = n[ii],
      "kappa" = kappa,
      "mean_residual_error" = mean_residual_error
    )
  }

  data <- data.table::rbindlist(data)
  return(data)
}



.get_test_statistic_lookup_table <- function(
    manuscript_dir,
    k = 0.80,
    for_manuscript = TRUE,
    with_outliers = FALSE
) {
  if (for_manuscript) {
    alpha_levels <- rev(1.0 - c(0.10, 0.20, 0.50, 0.80, 0.90, 0.95, 0.975, 0.99, 0.999))
    n <- c(5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000)
  } else {
    alpha_levels <- rev(1.0 - c(
      0.00, 0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55,
      0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95,
      0.96, 0.97, 0.975, 0.98, 0.9825, 0.985, 0.9875, 0.99,
      0.991, 0.992, 0.993, 0.994, 0.995, 0.996, 0.997, 0.998, 0.999, 0.9995
    ))
    n <- c(
      5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200,
      225, 250, 275, 300, 340, 380, 420, 460, 500,
      560, 620, 680, 740, 800, 900, 1000, 1100, 1200, 1500,
      2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)
  }

  data <- .get_test_statistics_data(manuscript_dir, with_outliers = with_outliers)
  data <- data[kappa == k]
  data[, "kappa" := NULL]

  # Compute alpha. Higher alpha yields higher tau (with n constant).
  data[, "alpha" := 1.0 - (seq_len(.N) - 1L) / (.N - 1L), by = c("n")]

  # Define all combinations.
  new_data <- data.table::as.data.table(expand.grid(list("alpha" = alpha_levels, "n" = n)))
  interp_list <- apply(new_data, function(ii) as.list(ii), MARGIN = 1L, simplify = FALSE)

  # Interpolate critical statistic values.
  interp_tau <- sapply(interp_list, power.transform:::.interpolate_2d, data = data)
  new_data[, "tau" := interp_tau]

  if (for_manuscript) {
    # Represent as 10^-2 for better interpretability.
    new_data[, "tau" := round(tau * 100.0, 2L)]
    return(data.table::dcast(data = new_data, formula = n ~ alpha, value.var = "tau"))

  } else {
    return(new_data)
  }
}

