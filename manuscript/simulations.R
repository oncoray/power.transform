.get_shifted_distribution_data <- function(manuscript_dir){

  # Set seed.
  set.seed(19L)

  if(!file.exists(file.path(manuscript_dir, "shifted_distributions_plot.RDS"))){

    # Normal distribution.
    x_normal <- power.transform::ragn(10000L, location=0, scale=1/sqrt(2), alpha=0.5, beta=2)

    # Right skewed data
    x_right_skewed <- power.transform::ragn(10000L, location=0, scale=1/sqrt(2), alpha=0.2, beta=2)

    # Left skewed data
    x_left_skewed <- power.transform::ragn(10000L, location=0, scale=1/sqrt(2), alpha=0.8, beta=2)


    # generator ----------------------------------------------------------------
    generate_experiment_data <- coro::generator(
      function(
    x_normal,
    x_right_skewed,
    x_left_skewed){

        shift_range <- 10^seq(from=0, to=6, by=0.1)

        for(d in shift_range){
          for(distribution in c("normal", "right-skewed", "left-skewed")){
            if(distribution == "normal"){
              x <- x_normal

            } else if(distribution == "right-skewed"){
              x <- x_right_skewed

            } else if(distribution == "left-skewed"){
              x <- x_left_skewed
            }

            for(method in c("box_cox", "yeo_johnson")){
              for(shift in c(FALSE, TRUE)){
                for(robust in c(FALSE)){
                  for(estimation_method in setdiff(power.transform:::..estimators_all(), power.transform:::..estimators_raymaekers_robust())){

                    # Skip if robust, but not shifted.
                    if(robust && !shift) next

                    if(!shift && !robust){
                      version <- "original"

                    } else if(!robust){
                      version <- "shift-sensitive"

                    } else {
                      version <- "robust shift-sensitive"
                    }

                    yield(list(
                      "x" = x + d,
                      "d" = d,
                      "distribution" = distribution,
                      "method" = method,
                      "shift" = shift,
                      "robust" = robust,
                      "version" = version,
                      "estimation_method" = estimation_method
                    ))
                  }
                }
              }
            }
          }
        }
      }
    )

    .compute_lambda <- function(parameter_set){
      # Create transformer object.
      if(parameter_set$version == "original"){
        transformer <- suppressWarnings(
          power.transform::find_transformation_parameters(
            x = parameter_set$x,
            method = parameter_set$method,
            robust = parameter_set$robust,
            shift = parameter_set$shift,
            estimation_method = parameter_set$estimation_method,
            lambda = NULL,
            optimiser_control = list("xtol_rel"=1e-5)))

      } else {
        transformer <- suppressWarnings(
          power.transform::find_transformation_parameters(
            x = parameter_set$x,
            method = parameter_set$method,
            robust = parameter_set$robust,
            shift = parameter_set$shift,
            estimation_method = parameter_set$estimation_method))
      }


      return(
        data.table::data.table(
          "distribution" = parameter_set$distribution,
          "method" = parameter_set$method,
          "estimation_method" = parameter_set$estimation_method,
          "version" = parameter_set$version,
          "d" = log10(parameter_set$d),
          "lambda" = transformer@lambda,
          "x_0" = transformer@shift
        ))
    }

    # computations -------------------------------------------------------------

    # Generate all experiments.
    experiments <- coro::collect(generate_experiment_data(
      x_normal = x_normal,
      x_right_skewed = x_right_skewed,
      x_left_skewed = x_left_skewed
    ))

    # data <- lapply(
    #   X = experiments,
    #   FUN = .compute_lambda)

    # Start cluster
    cl <- parallel::makeCluster(16L)

    # Compute all data in parallel.
    data <- parallel::parLapply(
      cl=cl,
      X=experiments,
      fun=.compute_lambda
    )

    # Stop cluster.
    parallel::stopCluster(cl)

    # Combine data into a single table.
    data <- data.table::rbindlist(data)

    # Save to file.
    saveRDS(data, file.path(manuscript_dir, "shifted_distributions_plot.RDS"))
  } else {
    data <- readRDS(file.path(manuscript_dir, "shifted_distributions_plot.RDS"))
  }

  # Update distribution, method and version to factors.
  data$distribution <- factor(
    x = data$distribution,
    levels = c("normal", "right-skewed", "left-skewed"))
  data$method <- factor(
    x = data$method,
    levels = c("box_cox", "yeo_johnson"),
    labels = c("Box-Cox", "Yeo-Johnson"))
  data$version <- factor(
    x = data$version,
    levels = c("original", "shift-sensitive"))
  data$estimation_method <- factor(
    x = data$estimation_method,
    levels = c("mle", "anderson_darling", "cramer_von_mises", "jarque_bera", "dagostino"),
    labels = c("MLE", "Anderson-Darling", "Cramér-von Mises", "Jarque-Bera", "D'Agostino"))

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

      if(is.null(transformation_methods)){
        transformation_methods <- c("box_cox", "yeo_johnson")
      }

      if(is.null(estimation_methods)){
        estimation_methods <- setdiff(power.transform:::..estimators_all(), power.transform:::..estimators_raymaekers_robust())
      }

      for(estimation_method in estimation_methods){

        for(transformation_method in transformation_methods){
          # Non-robust transformation.
          fun_args <- list(
            "method" = transformation_method,
            "estimation_method" = estimation_method,
            "robust" = FALSE,
            "shift" = TRUE)

          opt_args <- list()

          file_name <- file.path(
            manuscript_dir,
            file_dir,
            paste0(
              transformation_method, "_",
              estimation_method, "_",
              "non_robust.RDS"))

          if(!file.exists(file_name)){
            # Only yield arguments if the file does not exist.
            yield(list(
              "name" = "non-robust",
              "method" = transformation_method,
              "estimation_method" = estimation_method,
              "file_name" = file_name,
              "fun_args" = fun_args,
              "opt_args" = opt_args))
          } else {
            yield(file_name)
          }
        }

        for(robustness_source in c("empirical_probability", "transformed", "residual")){

          def_opt_limits <- list("k1"=c(0.0, 10.0), "k2"=c(0.0, 10.0))

          if(robustness_source == "empirical_probability"){
            def_opt_limits$k1[2] <- 1.0
            def_opt_limits$k2[2] <- 1.0

            def_opt_init <- list("k1" = 0.80, "k2" = 0.95)

          } else if(robustness_source == "transformed"){
            def_opt_init <- list("k1" = 1.28, "k2" = 1.96)

          } else if(robustness_source == "residual"){
            def_opt_init <- list("k1" = 0.50, "k2" = 1.0)
          }

          for(transformation_method in transformation_methods){
            for(robustness_weighting_function in c("step", "triangle", "cosine")){

              fun_args <- list(
                "method" = transformation_method,
                "estimation_method" = estimation_method,
                "robust" = TRUE,
                "shift" = TRUE,
                "weighting_function" = paste0(robustness_source, "_", robustness_weighting_function),
                "backup_use_default" = FALSE)

              opt_limits <- def_opt_limits
              opt_init <- def_opt_init

              # Step-weighting only has k1 as a parameter.
              if(robustness_weighting_function == "step"){
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
                  robustness_weighting_function, ".RDS"))

              if(!file.exists(file_name)){
                # Only yield arguments if the file does not exist.
                yield(list(
                  "name" = paste0(robustness_source, "-", robustness_weighting_function),
                  "method" = transformation_method,
                  "estimation_method" = estimation_method,
                  "file_name" = file_name,
                  "fun_args" = fun_args,
                  "opt_args" = list("initial" = opt_init, "limits" = opt_limits)))

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
    estimation_methods = NULL){

    if(is.null(transformation_methods)){
      transformation_methods <- c("box_cox", "yeo_johnson")
    }

    if(is.null(estimation_methods)){
      estimation_methods <- setdiff(power.transform:::..estimators_all(), power.transform:::..estimators_raymaekers_robust())
    }

    lambda_list <- list()
    ii <- 1L
    for(transformation_method in transformation_methods){
      for(estimation_method in estimation_methods){
        # Create transformer object.
        transformer <- suppressWarnings(
          power.transform::find_transformation_parameters(
            x = x,
            method = transformation_method,
            estimation_method = estimation_method,
            robust = FALSE,
            shift = TRUE))

        lambda_list[[ii]] <- data.table::data.table(
          "method" = transformation_method,
          "estimation_method" = estimation_method,
          "lambda" = transformer@lambda)

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
    side){

    # Set parameters.
    parameters <- list(
      "n" = n,
      "alpha" = alpha,
      "beta" = beta,
      "target_lambda" = target_lambda)

    # Add outliers.
    x <- lapply(
      outlier_fraction,
      function(k, x, parameters){
        # Set parameters.
        parameters$k <- k

        # Compute interquartile range.
        interquartile_range <- stats::IQR(x)

        # Compute upper and lower quartiles.
        q_lower <- stats::quantile(x, probs=0.25, names=FALSE)
        q_upper <- stats::quantile(x, probs=0.75, names=FALSE)

        # Set data where the outliers will be copied into.
        x_outlier <- x

        if(k != 0.0){
          n_draw <- ceiling(k * length(x))

          # Generate outlier values that are smaller than Q1 - 1.5 IQR or larger
          # than Q3 + 1.5 IQR.
          if (side == "both") {
            x_random <- stats::runif(n_draw, min=-2.0, max=2.0)
          } else if (side == "right") {
            x_random <- stats::runif(n_draw, min=0.0, max=2.0)
          } else if (side == "left") {
            x_random <- stats::runif(n_draw, min=-2.0, max=0.0)
          } else {
            stop(paste0("side was not recognised: ", side))
          }

          outlier <- numeric(n_draw)
          if(any(x_random < 0)){
            outlier[x_random < 0] <- q_lower - 1.5 * interquartile_range + x_random[x_random < 0] * interquartile_range
          }

          if(any(x_random >= 0)){
            outlier[x_random >= 0] <- q_upper + 1.5 * interquartile_range + x_random[x_random >= 0] * interquartile_range
          }

          # Randomly insert outlier values.
          x_outlier[sample(seq_along(x), size=n_draw, replace=FALSE)] <- outlier
        }

        return(list(
          "x" = x_outlier,
          "parameters" = parameters))
      },
      x = x,
      parameters)

    return(x)
  }

  # Helper function for computing transformer parameters under optimisation
  # constraints.
  .compute_robust_lambda <- function(
    x,
    fun_args,
    opt_args){

    # Custom parser.
    ..outlier_parser <- function(
    x,
    fun_args,
    opt_args){

      # Prevent notes
      method <- estimation_method <- NULL

      # Create transformer.
      transformer <- suppressWarnings(do.call(
        power.transform::find_transformation_parameters,
        args=c(
          list("x" = x$x),
          list("weighting_function_parameters" = opt_args),
          fun_args)))

      target_lambda <- x$parameters$target_lambda[method == fun_args$method & estimation_method == fun_args$estimation_method]$lambda

      parameter_data <- c(
        opt_args,
        list(
          "method" = fun_args$method,
          "estimation_method" = fun_args$estimation_method,
          "weight_method" = fun_args$weighting_function,
          "lambda" = transformer@lambda,
          "target_lambda" = target_lambda,
          "shift" = transformer@shift)
      )

      return(parameter_data)
    }

    parameter_data <- lapply(
      x,
      ..outlier_parser,
      fun_args = fun_args,
      opt_args = opt_args)

    return(parameter_data)
  }

  # Outer optimisation function used to optimise the weighting function parameters
  # over a dataset.
  .optimisation_outer <- function(
    experiment,
    data,
    verbose = FALSE){

    set.seed(9)

    results_list <- list(
      "name" = experiment$name,
      "method" = experiment$method,
      "estimation_method" = experiment$estimation_method)

    # Run optimiser.
    if(length(experiment$opt_args) > 0){

      x0 <- unlist(experiment$opt_args$initial)

      lower <- experiment$opt_args$limits$k1[1]
      upper <- experiment$opt_args$limits$k1[2]

      if(length(experiment$opt_args$limits) > 1){
        lower <- c(lower, experiment$opt_args$limits$k2[1])
        upper <- c(upper, experiment$opt_args$limits$k2[2])
      }

      h <- nloptr::sbplx(
        x0 = x0,
        fn = .optimisation_inner,
        lower = lower,
        upper = upper,
        control=list("xtol_rel"=1e-2, "maxeval"=50),
        x = data,
        fun_args = experiment$fun_args,
        verbose = verbose)

      results_list <- c(results_list, list("k1" = h$par[1]))
      if(length(h$par) > 1) results_list <- c(results_list, list("k2" = h$par[2]))
      results_list <- c(results_list, list("value" = h$value))

    } else {
      value <- .optimisation_inner(
        opt_param=list(),
        x = data,
        fun_args=experiment$fun_args)

      results_list <- c(results_list, list("value" = value))
    }

    # Save to file.
    saveRDS(results_list, file=experiment$file_name)

    return(invisible(TRUE))
  }

  # Inner optimisation function that sets the loss.
  .optimisation_inner <- function(
    opt_param,
    x,
    fun_args,
    verbose = FALSE){

    # Parse options.
    opt_args <- list()
    if(length(opt_param) >= 1) opt_args <- c(opt_args, list("k1"=opt_param[1]))
    if(length(opt_param) >= 2) opt_args <- c(opt_args, list("k2"=opt_param[2]))

    parameter_data <- .compute_robust_lambda(
      x = x,
      fun_args = fun_args,
      opt_args = opt_args)

    lambda_error <- sapply(
      parameter_data,
      function(x) (abs(x$lambda - x$target_lambda)))

    if(verbose){
      message <- paste0("target: ", mean(lambda_error))
      if(length(opt_args) > 0) message <- c(message, paste0("k1: ", opt_param[1]))
      if(length(opt_args) > 1) message <- c(message, paste0("k2: ", opt_param[2]))
      cat(paste0(message, collapse="; "))
      cat("\n")
    }

    return(mean(lambda_error))
  }

  # computations ---------------------------------------------------------------

  # Select experiments to conduct by filtering out experiments that are file
  # names.
  experiments <- coro::collect(experiment_args(
    manuscript_dir = manuscript_dir,
    side = side))
  experiments <- experiments[!sapply(experiments, is.character)]

  if(length(experiments) > 0){
    set.seed(95)

    # Generate table of asymmetric generalised normal distribution parameters.
    n_distributions <- 100L

    # Generate alpha, beta and n.
    n <- stats::runif(n=n_distributions, min=2, max=4)
    n <- ceiling(10^n)

    alpha <- stats::runif(n=n_distributions, min=0.01, max=0.99)
    beta <- stats::runif(n=n_distributions, min=1.00, max=5.00)

    # Generate corresponding distributions.
    x <- mapply(
      power.transform::ragn,
      n = n,
      alpha = alpha,
      beta = beta)

    # Compute lambda values without weighting.
    target_lambda <- lapply(x, .compute_target_lambda)

    # Add outliers, between 0.00 and 0.10.
    n_outliers <- 10
    outlier_fraction <- (seq_len(n_outliers) - 1)^2 /((n_outliers-1)^2 * 10)

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
        "side" = side),
      SIMPLIFY = FALSE,
      USE.NAMES = FALSE)

    x <- unlist(x, recursive = FALSE)

    # lapply(experiments, .optimisation_outer, data=x, verbose=TRUE)

    # Start cluster
    cl <- parallel::makeCluster(18L)

    parallel::clusterExport(
      cl=cl,
      varlist = c(".compute_robust_lambda", ".optimisation_inner"),
      envir = environment()
    )

    parallel::parLapplyLB(
      cl = cl,
      X = experiments,
      fun = .optimisation_outer,
      data = x,
      chunk.size = 1L)

    # Stop cluster.
    parallel::stopCluster(cl)
  }

  # Collect all files
  experiments <- coro::collect(experiment_args(
    manuscript_dir = manuscript_dir,
    side = side))
  weighting_parameters <- lapply(experiments, readRDS)
  weighting_parameters <- lapply(weighting_parameters, data.table::as.data.table)
  weighting_parameters <- data.table::rbindlist(weighting_parameters, fill = TRUE)

  return(weighting_parameters)
}



.get_transformation_parameters_two_sided <- function(manuscript_dir){

  if(!file.exists(file.path(manuscript_dir, "robustness_comparison_marginal_parameter_data.RDS"))){

    # generator ----------------------------------------------------------------
    experiment_args <- coro::generator(
      function(
    manuscript_dir,
    transformation_methods = NULL,
    estimation_methods = NULL){

        if(is.null(transformation_methods)){
          transformation_methods <- c("box_cox", "yeo_johnson")
        }

        if(is.null(estimation_methods)){
          estimation_methods <- setdiff(power.transform:::..estimators_all(), power.transform:::..estimators_raymaekers_robust())
        }

        for(estimation_method in estimation_methods){

          for(transformation_method in transformation_methods){
            # Non-robust transformation.
            fun_args <- list(
              "method" = transformation_method,
              "estimation_method" = estimation_method,
              "robust" = FALSE,
              "shift" = TRUE)

            file_name <- file.path(
              manuscript_dir,
              "robustness_comparison",
              paste0(
                transformation_method, "_",
                estimation_method, "_",
                "non_robust.RDS"))

            # Only yield arguments if the file does not exist.
            yield(list(
              "name" = "non-robust",
              "method" = transformation_method,
              "estimation_method" = estimation_method,
              "file_name" = file_name,
              "fun_args" = fun_args))
          }

          for(robustness_source in c("empirical_probability", "transformed", "residual")){
            for(transformation_method in transformation_methods){
              for(robustness_weighting_function in c("step", "triangle", "cosine")){

                file_name <- file.path(
                  manuscript_dir,
                  "robustness_comparison",
                  paste0(
                    transformation_method, "_",
                    estimation_method, "_",
                    robustness_source, "_",
                    robustness_weighting_function, ".RDS"))

                optimisation_settings <- readRDS(file_name)
                weighting_function_parameters <- list("k1"=optimisation_settings$k1)
                if(!is.null(optimisation_settings$k2)){
                  weighting_function_parameters <- c(
                    weighting_function_parameters,
                    list("k2"=optimisation_settings$k2))
                }

                fun_args <- list(
                  "method" = transformation_method,
                  "estimation_method" = estimation_method,
                  "robust" = TRUE,
                  "shift" = TRUE,
                  "weighting_function" = paste0(robustness_source, "_", robustness_weighting_function),
                  "weighting_function_parameters" = weighting_function_parameters)

                # Only yield arguments if the file does not exist.
                yield(list(
                  "name" = paste0(robustness_source, "-", robustness_weighting_function),
                  "method" = transformation_method,
                  "estimation_method" = estimation_method,
                  "file_name" = file_name,
                  "fun_args" = fun_args))
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
    estimation_methods = NULL){

      if(is.null(transformation_methods)){
        transformation_methods <- c("box_cox", "yeo_johnson")
      }

      if(is.null(estimation_methods)){
        estimation_methods <- setdiff(power.transform:::..estimators_all(), power.transform:::..estimators_raymaekers_robust())
      }

      lambda_list <- list()
      ii <- 1L
      for(transformation_method in transformation_methods){
        for(estimation_method in estimation_methods){
          # Create transformer object.
          transformer <- suppressWarnings(
            power.transform::find_transformation_parameters(
              x = x,
              method = transformation_method,
              estimation_method = estimation_method,
              robust = FALSE,
              shift = TRUE))

          lambda_list[[ii]] <- data.table::data.table(
            "method" = transformation_method,
            "estimation_method" = estimation_method,
            "lambda" = transformer@lambda)

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
    ii,
    outlier_fraction){

      # Set parameters.
      parameters <- list(
        "n" = n,
        "alpha" = alpha,
        "beta" = beta,
        "target_lambda" = target_lambda,
        "ii" = ii)

      # Add outliers.
      x <- lapply(
        outlier_fraction,
        function(k, x, parameters){
          # Set parameters.
          parameters$k <- k

          # Compute interquartile range.
          interquartile_range <- stats::IQR(x)

          # Compute upper and lower quartiles.
          q_lower <- stats::quantile(x, probs=0.25, names=FALSE)
          q_upper <- stats::quantile(x, probs=0.75, names=FALSE)

          # Set data where the outliers will be copied into.
          x_outlier <- x

          if(k != 0.0){
            n_draw <- ceiling(k * length(x))

            # Generate outlier values that are smaller than Q1 - 1.5 IQR or larger
            # than Q3 + 1.5 IQR.
            x_random <- stats::runif(n_draw, min=-2.0, max=2.0)
            outlier <- numeric(n_draw)
            if(any(x_random < 0)){
              outlier[x_random < 0] <- q_lower - 1.5 * interquartile_range + x_random[x_random < 0] * interquartile_range
            }

            if(any(x_random >= 0)){
              outlier[x_random >= 0] <- q_upper + 1.5 * interquartile_range + x_random[x_random >= 0] * interquartile_range
            }

            # Randomly insert outlier values.
            x_outlier[sample(seq_along(x), size=n_draw, replace=FALSE)] <- outlier
          }

          return(list(
            "x" = x_outlier,
            "parameters" = parameters))
        },
        x = x,
        parameters)

      return(x)
    }



    # Helper function for computing transformer parameters under optimisation
    # constraints.
    .compute_robust_lambda <- function(
    experiment,
    data){

      # Custom parser.
      ..outlier_parser <- function(
    x,
    fun_args){
        # Create transformer.
        transformer <- suppressWarnings(do.call(
          power.transform::find_transformation_parameters,
          args=c(
            list("x" = x$x),
            fun_args)))

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
          "ii" = x$parameters$ii)

        return(parameter_data)
      }

      parameter_data <- lapply(
        data,
        ..outlier_parser,
        fun_args = experiment$fun_args)

      parameter_data <- data.table::rbindlist(
        parameter_data,
        fill = TRUE)

      return(parameter_data)
    }

    # computations -------------------------------------------------------------
    set.seed(95)

    # Generate table of asymmetric generalised normal distribution parameters.
    n_distributions <- 100L

    # Generate alpha, beta and n.
    n <- stats::runif(n=n_distributions, min=2, max=4)
    n <- ceiling(10^n)

    alpha <- stats::runif(n=n_distributions, min=0.01, max=0.99)
    beta <- stats::runif(n=n_distributions, min=1.00, max=5.00)

    # Generate corresponding distributions.
    x <- mapply(
      power.transform::ragn,
      n = n,
      alpha = alpha,
      beta = beta)

    # Compute lambda values without weighting.
    target_lambda <- lapply(x, .compute_target_lambda)

    # Add outliers, between 0.00 and 0.10.
    n_outliers <- 10
    outlier_fraction <- (seq_len(n_outliers) - 1)^2 /((n_outliers-1)^2 * 10)

    # Generate outliers in the data.
    x <- mapply(
      FUN = .populate_outliers,
      x = x,
      n = n,
      alpha = alpha,
      beta = beta,
      target_lambda = target_lambda,
      ii = seq_along(x),
      MoreArgs=list("outlier_fraction"=outlier_fraction),
      SIMPLIFY=FALSE,
      USE.NAMES=FALSE)

    x <- unlist(x, recursive = FALSE)

    experiments <- coro::collect(experiment_args(manuscript_dir = manuscript_dir))

    cl <- parallel::makeCluster(18L)

    data <- parallel::parLapplyLB(
      cl = cl,
      X = experiments,
      fun = .compute_robust_lambda,
      data = x,
      chunk.size = 1L)

    # Stop cluster.
    parallel::stopCluster(cl)

    # Merge to data.table.
    data <- data.table::rbindlist(
      data,
      fill=TRUE)

    saveRDS(data, file=file.path(manuscript_dir, "robustness_comparison_marginal_parameter_data.RDS"))

  } else {
    data <- readRDS(data, file=file.path(manuscript_dir, "robustness_comparison_marginal_parameter_data.RDS"))
  }

  data[is.na(weight_method), "weight_method" := "none"]

  data$method <- factor(
    x = data$method,
    levels = c("box_cox", "yeo_johnson"),
    labels = c("Box-Cox", "Yeo-Johnson"))
  data$weight_method <- factor(
    x = data$weight_method,
    levels = c(
      "none",
      "empirical_probability_step", "empirical_probability_triangle", "empirical_probability_cosine",
      "transformed_step", "transformed_triangle", "transformed_cosine",
      "residual_step", "residual_triangle", "residual_cosine"),
    labels = c(
      "non-robust",
      "emp. prob. (step)", "emp. prob. (triangle)", "emp. prob. (tap. cos.)",
      "z-score (step)", "z-score (triangle)", "z-score (tap. cos.)",
      "residual (step)", "residual (triangle)", "residual (tap. cos.)"))
  data$estimation_method <- factor(
    x = data$estimation_method,
    levels = c("mle", "anderson_darling", "cramer_von_mises", "jarque_bera", "dagostino"),
    labels = c("MLE", "Anderson-Darling", "Cramér-von Mises", "Jarque-Bera", "D'Agostino"))

  return(data)
}
