#' @include ParameterEstimators.R
NULL



# estimatorEmpiricalDistributionFunction definition ----------------------------
setClass(
  "estimatorEmpiricalDistributionFunction",
  contains = "estimatorGeneric")

# estimatorAndersonDarling definition ------------------------------------------
setClass(
  "estimatorAndersonDarling",
  contains = "estimatorEmpiricalDistributionFunction")

# estimatorCramervonMises definition -------------------------------------------
setClass(
  "estimatorCramervonMises",
  contains = "estimatorEmpiricalDistributionFunction")

# estimatorKolmogorovSmirnov definition -----------------------------------------
setClass(
  "estimatorKolmogorovSmirnov",
  contains = "estimatorEmpiricalDistributionFunction")



# .compute_objective (general EDF) ---------------------------------------------
setMethod(
  ".compute_objective",
  signature(object="estimatorEmpiricalDistributionFunction"),
  function(object, transformer, x, ...){

    # Prevent NOTE due to non-standard evaluation in data.table.
    p_empirical <- NULL

    # Make sure that x is sorted.
    if(is.unsorted(x)) x <- sort(x)

    # Transform feature values
    y <- ..transform(
      object = object,
      x = x)

    if(any(!is.finite(y))) return(NA_real_)

    # Compute (merged) empirical probabilities. Average empirical probabilities
    # for when x has multiple values. Though this necessitates using the
    # data.table package, this is by far the fastest implementation.
    p_empirical <- (seq_along(x) - 1/3) / (length(x) + 1/3)

    # Set up a data.table.
    data <- data.table::data.table("y"=y, "p"=p_empirical)
    data[, "p_group":=mean(p_empirical), by="y"]
    p_empirical <- data$p_group

    # Determine mu and sigma for putative normal distribution from the data.
    # This is necessary to compute the probabilities expected by theory.
    if(transformer@robust){
      # Compute M-estimates for locality and scale
      estimates <- huber_estimate(y, tol=1E-3)

    } else {
      estimates <- list("mu" = mean(y), "sigma" = stats::sd(y))
    }

    # Check problematic values.
    if(!is.finite(estimates$sigma)) return(NA_real_)
    if(estimates$sigma == 0.0) return(NA_real_)

    # Compute expected probabilities according to the cumulative density
    # function of the normal distribution parametrised by mu and sigma.
    p_expected <- stats::pnorm(
      q = y,
      mean = estimates$mu,
      sd = estimates$sigma)

    return(list(
      "p_empirical" = p_empirical,
      "p_expected" = p_expected))
  }
)



# .compute_objective (Anderson-Darling) ----------------------------------------
setMethod(
  ".compute_objective",
  signature(object="estimatorAndersonDarling"),
  function(object, transformer, x, ...){
    # Based on the Anderson-Darling test statistic.

    # Get general EDF data first.
    p <- callNextMethod()

    # Get weights from the weighting function..
    w <- .get_weights(
      object = object,
      transformer = transformer,
      x = x)

    if(!all(is.finite(w))) return(NA_real_)

    sum_w <- sum(w)
    if(sum_w == 0.0) return(NA_real_)

    # Get weights for Anderson-Darling.
    w_ad <- 1.0 / (p$p_expected * (1.0 - p$p_expected))

    # Compute (weighted) statistic.
    t <- w * w_ad * (p$p_empirical - p$p_expected)^2

    # Normalise by weights. This is to prevent the optimiser from cheating by
    # finding fitting parameters where most of the weights will be (close to)
    # zero.
    t <- t / sum_w

    # Statistic should be minimised for better fits.
    return(t)
  }
)



# .compute_objective (Cramér-von Mises) ----------------------------------------
setMethod(
  ".compute_objective",
  signature(object="estimatorCramervonMises"),
  function(object, transformer, x, ...){
    # Based on the Cramér-von Mises test statistic.

    # Get general EDF data first.
    p <- callNextMethod()

    # Get weights from the weighting function..
    w <- .get_weights(
      object = object,
      transformer = transformer,
      x = x)

    if(!all(is.finite(w))) return(NA_real_)

    sum_w <- sum(w)
    if(sum_w == 0.0) return(NA_real_)

    # Compute (weighted) statistic.
    t <- w * (p$p_empirical - p$p_expected)^2

    # Normalise by weights. This is to prevent the optimiser from cheating by
    # finding fitting parameters where most of the weights will be (close to)
    # zero.
    t <- t / sum_w

    # Statistic should be minimised for better fits.
    return(t)
  }
)



# .compute_objective (Kolmogorov-Smirnov) --------------------------------------
setMethod(
  ".compute_objective",
  signature(object="estimatorKolmogorovSmirnov"),
  function(object, transformer, x, ...){
    # Based on the Kolmogorov-Smirnov test statistic.

    # Get general EDF data first.
    p <- callNextMethod()

    # Get weights from the weighting function..
    w <- .get_weights(
      object = object,
      transformer = transformer,
      x = x)

    if(!all(is.finite(w))) return(NA_real_)
    if(all(w == 0.0)) return(NA_real_)

    # Compute absolute error between empirical and theoretic probabilities.
    t <- abs(p$p_empirical - p$p_expected)

    # Find maximum value. The weights are only used for finding the maximum
    # value, and unlike other EDF-based tests we do not normalise by the sum of
    # weights.
    t <- t[which.max(w * t)]

    # Statistic should be minimised for better fits.
    return(t)
  }
)



..estimators_anderson_darling <- function(){
  return(c("anderson_darling"))
}

..estimators_mises_von_cramer <- function(){
  return(c("mises_von_cramer"))
}

..estimators_kolmogorov_smirnov <- function(){
  return(c("kolmogorov_smirnov"))
}
