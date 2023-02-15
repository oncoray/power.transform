#' @include WeightingFunctions.R

# ..get_default_weighting_parameters (none) ------------------------------------
setMethod(
  "..get_default_weighting_parameters",
  signature(object = "weightingMethodNone"),
  function(object, transformer, estimator, ...){

    # No weights required.
    return(NULL)
  }
)



# ..get_default_weighting_parameters (emp prob, step, general) -----------------
setMethod(
  "..get_default_weighting_parameters",
  signature(object = "weightingMethodEmpiricalProbabilityStep"),
  function(object, transformer, estimator, ...){

    values <- list("k1" = 0.90)

    if(is(transformer, "transformationBoxCox")){
      if(is(estimator, "estimatorMaximumLikelihoodEstimation")){
        values <- list("k1" = 0.84)
      }

    } else if(is(transformer, "transformationYeoJohnson")){
      if(is(estimator, "estimatorMaximumLikelihoodEstimation")){
        values <- list("k1" = 0.92)
      }
    }

    return(values)
  }
)



# ..get_default_weighting_parameters (emp prob, triangle, general) -------------
setMethod(
  "..get_default_weighting_parameters",
  signature(object = "weightingMethodEmpiricalProbabilityTriangle"),
  function(object, transformer, estimator, ...){

    values <- list("k1" = 0.85, "k2" = 0.95)

    if(is(transformer, "transformationBoxCox")){
      if(is(estimator, "estimatorMaximumLikelihoodEstimation")){
        values <- list("k1" = 0.78, "k2" = 0.99)
      }

    } else if(is(transformer, "transformationYeoJohnson")){
      if(is(estimator, "estimatorMaximumLikelihoodEstimation")){
        values <- list("k1" = 0.88, "k2" = 0.93)
      }
    }

    return(values)
  }
)



# ..get_default_weighting_parameters (emp prob, cosine, general) ---------------
setMethod(
  "..get_default_weighting_parameters",
  signature(object = "weightingMethodEmpiricalProbabilityCosine"),
  function(object, transformer, estimator, ...){

    values <- list("k1" = 0.85, "k2" = 0.95)

    if(is(transformer, "transformationBoxCox")){
      if(is(estimator, "estimatorMaximumLikelihoodEstimation")){
        values <- list("k1" = 0.76, "k2" = 0.98)
      }

    } else if(is(transformer, "transformationYeoJohnson")){
      if(is(estimator, "estimatorMaximumLikelihoodEstimation")){
        values <- list("k1" = 0.89, "k2" = 0.93)
      }
    }

    return(values)
  }
)



# ..get_default_weighting_parameters (transformed, step, general) --------------
setMethod(
  "..get_default_weighting_parameters",
  signature(object = "weightingMethodTransformedStep"),
  function(object, transformer, estimator, ...){

    values <- list("k1" = 1.96)

    if(is(transformer, "transformationBoxCox")){
      if(is(estimator, "estimatorMaximumLikelihoodEstimation")){
        values <- list("k1" = 1.09)
      }

    } else if(is(transformer, "transformationYeoJohnson")){
      if(is(estimator, "estimatorMaximumLikelihoodEstimation")){
        values <- list("k1" = 1.04)
      }
    }

    return(values)
  }
)



# ..get_default_weighting_parameters (transformed, triangle, general) ----------
setMethod(
  "..get_default_weighting_parameters",
  signature(object = "weightingMethodTransformedTriangle"),
  function(object, transformer, estimator, ...){

    values <- list("k1" = 0.50, "k2" = 8.00)

    if(is(transformer, "transformationBoxCox")){
      if(is(estimator, "estimatorMaximumLikelihoodEstimation")){
        values <- list("k1" = 0.26, "k2" = 5.47)
      }

    } else if(is(transformer, "transformationYeoJohnson")){
      if(is(estimator, "estimatorMaximumLikelihoodEstimation")){
        values <- list("k1" = 0.12, "k2" = 9.91)
      }
    }

    return(values)
  }
)



# ..get_default_weighting_parameters (transformed, cosine, general) ------------
setMethod(
  "..get_default_weighting_parameters",
  signature(object = "weightingMethodTransformedCosine"),
  function(object, transformer, estimator, ...){

    values <- list("k1" = 0.50, "k2" = 8.00)

    if(is(transformer, "transformationBoxCox")){
      if(is(estimator, "estimatorMaximumLikelihoodEstimation")){
        values <- list("k1" = 0.03, "k2" = 6.91)
      }

    } else if(is(transformer, "transformationYeoJohnson")){
      if(is(estimator, "estimatorMaximumLikelihoodEstimation")){
        values <- list("k1" = 0.23, "k2" = 5.63)
      }
    }

    return(values)
  }
)



# ..get_default_weighting_parameters (residual, step, general) -----------------
setMethod(
  "..get_default_weighting_parameters",
  signature(object = "weightingMethodResidualStep"),
  function(object, transformer, estimator, ...){

    values <- list("k1" = 2.00)

    if(is(transformer, "transformationBoxCox")){
      if(is(estimator, "estimatorMaximumLikelihoodEstimation")){
        values <- list("k1" = 1.56)
      }

    } else if(is(transformer, "transformationYeoJohnson")){
      if(is(estimator, "estimatorMaximumLikelihoodEstimation")){
        values <- list("k1" = 10.00)
      }
    }

    return(values)
  }
)



# ..get_default_weighting_parameters (residual, triangle, general) -------------
setMethod(
  "..get_default_weighting_parameters",
  signature(object = "weightingMethodResidualTriangle"),
  function(object, transformer, estimator, ...){

    values <- list("k1" = 0.50, "k2" = 8.00)

    if(is(transformer, "transformationBoxCox")){
      if(is(estimator, "estimatorMaximumLikelihoodEstimation")){
        values <- list("k1" = 1.49, "k2" = 1.51)
      }

    } else if(is(transformer, "transformationYeoJohnson")){
      if(is(estimator, "estimatorMaximumLikelihoodEstimation")){
        values <- list("k1" = 9.96, "k2" = 10.00)
      }
    }

    return(values)
  }
)



# ..get_default_weighting_parameters (residual, cosine, general) ---------------
setMethod(
  "..get_default_weighting_parameters",
  signature(object = "weightingMethodResidualCosine"),
  function(object, transformer, estimator, ...){

    values <- list("k1" = 0.50, "k2" = 8.00)

    if(is(transformer, "transformationBoxCox")){
      if(is(estimator, "estimatorMaximumLikelihoodEstimation")){
        values <- list("k1" = 1.32, "k2" = 1.65)
      }

    } else if(is(transformer, "transformationYeoJohnson")){
      if(is(estimator, "estimatorMaximumLikelihoodEstimation")){
        values <- list("k1" = 9.97, "k2" = 10.00)
      }
    }

    return(values)
  }
)
