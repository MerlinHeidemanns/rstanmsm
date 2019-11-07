
coef.stan_msm <- function(object, digits = 1, kind = NULL){
  if (is.null(kind)){
     x <- cbind(object$coefs_median$alpha, object$ses$alpha)
     x <- rbind(x, cbind(object$coefs_median$beta, object$ses$beta))
     x <- rbind(x, cbind(object$coefs_median$gamma, object$ses$gamma))
     x <- rbind(x, cbind(object$coefs_median$delta, object$ses$delta))
     x <- rbind(x, cbind(object$coefs_median$lambda, object$ses$lambda))
  }
  if (!is.null(kind)){
    options <- c("alpha", "beta", "gamma", "delta", "lambda")
    if (is.element(kind, options)){
      x <- cbind(object$coefs_median[[kind]], object$ses[[kind]])
    } else {
      stop(cat("Please choose one of", paste0(options, sep = ","), "."))
    }
  }
  colnames(x) <- c("Median", "MAD_SD")
  return(round(x, digits = digits))
}

#' print.stan_msm
#'
#' @param object A stan_msm fit object.
#' @param digits The number of digits to round to.
#' @param detail Print output.
#'
#' A S3 function to print a summary output of a stan_msm object.

print.stan_msm <- function(object, digits = 1, detail = TRUE, ...){
  if (detail){
      cat(object$stan_function)
      cat("\n Continuous formula:      ", object$formula_continuous)
      cat("\n Observations:            ", length(object$data$n) )
    if (length(unique(object$data$n)) > 1){
      cat("\n Units:                   ",  length(unique(object$data$n)) )
      cat("\n min(T):                  ",  length(unique(object$data$n)) )
      cat("\n max(T):                  ",  length(unique(object$data$n)) )
      cat("\n N(states)                ",  object$data$K)
    }
    cat("\n ----------------------------\n")
    cat(" State-constant predictors\n ")
    print(coef(object, digits = digits, kind = "alpha"))
    cat("\n State-varying predictors\n ")
    print(coef(object, digits = digits, kind = "beta"))
  }
  if (detail) {
    cat("\n------\n")
    cat("* For help interpreting the printed output see ?print.stan_msm\n")
  }
  invisible(x)
}
