..estimator_wrapper <- function(
    par,
    transformer,
    estimator,
    x,
    parameter_type,
    verbose = FALSE,
    ...){

  # Set lambda, shift, as provided by the optimisation algorithm.
  if(setequal(parameter_type, c("lambda", "shift"))){
    transformer@shift <- par[1]
    transformer@lambda <- par[2]

  } else if(parameter_type == "lambda" & length(par) == 1){
    transformer@lambda <- par[1]

  } else if(parameter_type == "shift" & length(par) == 1){
    transformer@shift <- par[1]

  } else {
    stop("DEV: Optimisation function was not specified correctly.")
  }

  # Compute target value.
  target <- .compute_objective(
    object = estimator,
    transformer = transformer,
    x = x,
    ...)

  if(verbose){
    cat(paste0(
      "target: ", target, "; ",
      paste0(mapply(function(type, value) (paste0(type, ": ", value)), type=parameter_type, value=par), collapse="; "),
      "\n"))
  }

  return(target)
}



# estimatorGeneric definition --------------------------------------------------
setClass(
  "estimatorGeneric",
  slots = list(
    "weighting_method" = "ANY"),
  prototype = list(
    "weighting_method" = NULL))

setClass(
  "estimatorEmpiricalDistributionFunction",
  contains = "estimatorGeneric")

setClass(
  "estimatorAndersonDarling",
  contains = "estimatorEmpiricalDistributionFunction")

setClass(
  "estimatorCramervonMises",
  contains = "estimatorEmpiricalDistributionFunction")

setClass(
  "estimatorKolmogorovSmirnov",
  contains = "estimatorEmpiricalDistributionFunction")

setClass(
  "estimatorSkewnessKurtosis",
  contains = "estimatorGeneric")

setClass(
  "estimatorJarqueBera",
  contains = "estimatorSkewnessKurtosis")

setClass(
  "estimatorDAgostino",
  contains = "estimatorSkewnessKurtosis")



# .compute_objective (generic) -------------------------------------------------
setGeneric(
  ".compute_objective",
  function(object, ...) standardGeneric(".compute_objective"))



# .compute_objective (general) -------------------------------------------------
setMethod(
  ".compute_objective",
  signature(object = "estimatorGeneric"),
  function(object, ...){

    # Estimators should be specified for each estimator - an error is thrown to
    # ensure that this general method is never used.
    stop(paste0("DEV: missing .compute_objective method for the ", paste_s(class(object)), " class."))
  }
)



# .optimise_transformation_parameters (generic) --------------------------------
setGeneric(
  ".optimise_transformation_parameters",
  function(object, ...) standardGeneric(".optimise_transformation_parameters"))



# ..optimise_transformation_parameters (general) -------------------------------
setMethod(
  ".optimise_transformation_parameters",
  signature(object = "estimatorGeneric"),
  function(
    object,
    transformer,
    optimiser,
    optimisation_parameters,
    x,
    optimiser_control = list("xtol_rel"=1e-3),
    ...){

    # Check that we do not inadvertently pass problems that do not require
    # optimisation to the optimiser.
    if(is.null(optimisation_parameters)) stop("DEV: optimisation_parameters cannot be empty.")
    if(!is(transformer, "transformationPowerTransform")) stop("DEV: transformer should be a valid power transformation object.")
    if(!is.numeric(x)) stop("DEV: x should be numeric.")

    # Check fall-back option.
    if(!is_package_installed("nloptr") & optimiser %in% c("direct-l", "subplex", "nelder-mead")){
      warning(paste0(
        "The nloptr package is required to optimise power transformation parameters using the ",
        optimiser, " algoritm. stats::optim is used as a fallback option."))

      optimiser <- "optim-nelder-mead"
    }

    if(optimiser == "direct-l"){
      # DIRECT-L algorithm
      #
      # D. R. Jones, C. D. Perttunen, and B. E. Stuckmann, “Lipschitzian
      # optimization without the lipschitz constant,” J. Optimization Theory and
      # Applications, vol. 79, p. 157 (1993).
      #
      # J. M. Gablonsky and C. T. Kelley, “A locally-biased form of the DIRECT
      # algorithm," J. Global Optimization, vol. 21 (1), p. 27-37 (2001).

      results <- tryCatch(
        nloptr::directL(
          fn = ..estimator_wrapper,
          lower = optimisation_parameters$lower,
          upper = optimisation_parameters$upper,
          control = optimiser_control,
          transformer = transformer,
          estimator = object,
          x = x,
          parameter_type = optimisation_parameters$parameter_type),
        error = identity)

    } else if(optimiser == "subplex"){
      # SUBPLEX algorithm
      #
      # T. Rowan, “Functional Stability Analysis of Numerical
      # Algorithms”, Ph.D. thesis, Department of Computer Sciences, University
      # of Texas at Austin, 1990.

      results <- tryCatch(
        nloptr::sbplx(
          x0 = optimisation_parameters$initial,
          fn = ..estimator_wrapper,
          lower = optimisation_parameters$lower,
          upper = optimisation_parameters$upper,
          control = optimiser_control,
          transformer = transformer,
          estimator = object,
          x = x,
          parameter_type = optimisation_parameters$parameter_type),
        error = identity)

    } else if(optimiser == "nelder-mead"){
      # Nelder-Mead simplex algorithm
      #
      # J. A. Nelder and R. Mead, “A simplex method for function minimization,”
      # The Computer Journal 7, p. 308-313 (1965).
      #
      # M. J. Box, “A new method of constrained optimization and a comparison
      # with other methods,” Computer J. 8 (1), 42-52 (1965).

      results <- tryCatch(
        nloptr::neldermead(
          x0 = optimisation_parameters$initial,
          fn = ..estimator_wrapper,
          lower = optimisation_parameters$lower,
          upper = optimisation_parameters$upper,
          control = optimiser_control,
          transformer = transformer,
          estimator = object,
          x = x,
          parameter_type = optimisation_parameters$parameter_type),
        error = identity)

    } else if(optimiser == "optim-nelder-mead"){
      # Fall-back optimiser in case nloptr is not available. The Nelder-Mead
      # algorithm in stats::optim does not yield results as consistent as the
      # nloptr optimisers.
      results <- stats::optim(
        par = optimisation_parameters$initial,
        fn = ..estimator_wrapper,
        gr = NULL,
        transformer = transformer,
        estimator = object,
        x = x,
        parameter_type = optimisation_parameters$parameter_type,
        control = list(
          "abstol" = 1E-5,
          "reltol" = 1E-5))

    } else {
      stop(paste0("Optimiser not recognised: ", optimiser))
    }

    # Check for known errors. If an unknown error is encountered, raise this
    # error to the user.
    if(inherits(results, "error")){
      if(!results$message %in% c(
        "objective in x0 returns NA")){
        stop(results)
      }
    }

    # Initialise a parameter list to return to the calling process.
    parameter_list <- lapply(
      optimisation_parameters$parameter_type,
      function(x) (NA_real_))
    names(parameter_list) <- optimisation_parameters$parameter_type

    # Insert optimal parameter values into the parameter list, if any.
    if(!inherits(results, "error")){
      for(ii in seq_along(optimisation_parameters$parameter_type)){
        parameter_value <- results$par[ii]
        if(is.finite(parameter_value)){
          parameter_list[[optimisation_parameters$parameter_type[ii]]] <- parameter_value
        }
      }
    }

    return(parameter_list)
  }
)



# ..get_available_weighting_functions (generic) --------------------------------
setGeneric(
  "..get_available_weighting_functions",
  function(transformer, estimator, ...) standardGeneric("..get_available_weighting_functions"))



# ..get_available_weighting_functions (general) --------------------------------
setMethod(
  "..get_available_weighting_functions",
  signature(transformer = "transformationPowerTransform", estimator = "estimatorGeneric"),
  function(transformer, estimator, ...){

    # If transformation methods is not robust, do not use weights.
    if(!transformer@robust) return("none")

    return(..weighting_functions_all())
  }
)



# ..get_default_weighting_function (generic) -----------------------------------
setGeneric(
  "..get_default_weighting_function",
  function(transformer, estimator, ...) standardGeneric("..get_default_weighting_function"))



# ..get_default_weighting_function (general) -----------------------------------
setMethod(
  "..get_default_weighting_function",
  signature(transformer = "transformationPowerTransform", estimator = "estimatorGeneric"),
  function(transformer, estimator, ...){

    if(transformer@robust){
      return("empirical_probability_cosine")

    } else {
      return("none")
    }
  }
)



# .get_weights (generic) -------------------------------------------------------
setGeneric(
  ".get_weights",
  function(object, ...) standardGeneric(".get_weights"))



# .get_weights (general) -------------------------------------------------------
setMethod(
  ".get_weights",
  signature(object = "estimatorGeneric"),
  function(object, transformer, x, ...){
    return(.get_weights(
      object = object@weighting_method,
      transformer = transformer,
      x = x))
  }
)



.set_estimator <- function(
  transformer,
  estimation_method,
  weighting_function,
  weighting_function_parameters){

  # Check estimation method.
  available_estimators <- ..get_available_estimators(object=transformer)
  if(!estimation_method %in% available_estimators){
    stop(paste0(
      "The desired estimator ", estimation_method, " is not available for ",
      ifelse(transformer@robust, "robust ", ""),
      "optimisation for ", paste_s(class(transformer)), " objects."))
  }

  # Create estimator instance.
  if(estimation_method %in% ..estimators_mle()){
    estimator <- methods::new("estimatorMaximumLikelihoodEstimation")

  } else if(estimation_method %in% ..estimators_raymaekers_robust()){
    estimator <- methods::new("estimatorRaymaekersRobust")

  } else {
    stop(paste0("DEV: unknown estimation_method: ", estimation_method))
  }

  estimator@weighting_method <- .set_weighting_function(
    transformer = transformer,
    estimator = estimator,
    weighting_function = weighting_function,
    weighting_function_parameters = weighting_function_parameters)

  return(estimator)
}



..estimators_all <- function(){
  # Create a list with all estimators implemented in power.transform. Update
  # when more estimators are added.

  return(c(
    ..estimators_mle(),
    ..estimators_raymaekers_robust()
  ))
}
