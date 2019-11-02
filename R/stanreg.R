
# Create a stanreg object
#
# @param object A list provided by one of the \code{stan_*} modeling functions.
# @return A stanreg object.

stanreg <- function(object) {
  opt <- object$algorithm == "optimizing"
  stanfit <- object$stanfit
  family <- object$family
  y <- object$y
  x <- object$x
  parsed_formula <- object$parsed_formula
  names_list <- list(alpha = parsed_formula$d,
                     beta = parsed_formula$e,
                     gamma = parsed_formula$a,
                     delta = parsed_formula$b,
                     lambda = parsed_formula$c)

  nvars <- list(nvars_a = ncol(x$x_a), nvars_b = ncol(x$x_b),
                nvars_c = ncol(x$x_c), nvars_d = ncol(x$x_d),
                nvars_e = ncol(x$x_e))
  nobs <- NROW(y)

  stan_summary <- make_stan_summary(stanfit)
  end_coef <- grep("logalpha\\[.+\\]", rownames(stan_summary))[1] - 1 # first non parameter, position preceding it

  # get subset of interest
  coefs <- stan_summary[1:end_coef, select_median(object$algorithm)]
  stanmat <- as.matrix(stanfit)[, 1:end_coef, drop = FALSE]
  colnames(stanmat) <- names(coefs)
  ses <- apply(stanmat, 2L, mad)


  coefs <- split_naming(x = coefs, names_list = names_list, N = N, K = K)
  ses <- split_naming(x = ses, names_list = names_list, N = N, K = K)
  rownames(stan_summary) <- naming_fun(x = rownames(stan_summary), para_names = ses)

  covmat <- cov(stanmat)
  # rownames(covmat) <- colnames(covmat) <- rownames(stan_summary)[1:nrow(covmat)]
  if (object$algorithm == "sampling")
    check_rhats(stan_summary[, "Rhat"])

  # linear predictor, fitted values
  fit_extract <- rstan::extract(stanfit)
  mu <- fit_extract$y_hat
  res <- y - mu

  # coefs

  out <- nlist(
    coefs_median = coefs,
    ses = ses,
    fitted.values = mu,
    residuals = res,
    covmat,
    y,
    x = x,
    model = object$model,
    data = object$data,
    family,
    formula_discrete = object$formula_discrete,
    formula_continuous = object$formula_continuous,
    #prior.info = attr(stanfit, "prior.info"),
    algorithm = object$algorithm,
    stan_summary = stan_summary,
    stanfit = if (opt) stanfit$stanfit else stanfit,
    rstan_version = packageVersion("rstan"),
    call = object$call,
    stan_function = object$stan_function
  )

  return(structure(out, class = c("stanreg_msm", "stan_msm")))
}