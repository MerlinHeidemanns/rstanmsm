
# Create a stanreg object
#
# @param object A list provided by one of the \code{stan_*} modeling functions.
# @return A stanreg object.

stanreg <- function(object) {
  opt <- object$algorithm == "optimizing"
  stanfit <- object$stanfit
  family <- object$family
  data_lst <- object$data[["data_lst"]]
  names_lst <-object$data[["names_lst"]]
  y <- data_lst$y
  N <- data_lst$N
  K <- data_lst$K
  T <- data_lst$T

  stan_summary <- make_stan_summary(stanfit)
  end_coef <- grep("logalpha\\[.+\\]", rownames(stan_summary))[1] - 1 # first non parameter, position preceding it

  # get subset of interest
  coefs <- stan_summary[1:end_coef, select_median(object$algorithm)]
  stanmat <- as.matrix(stanfit)[, 1:end_coef, drop = FALSE]
  colnames(stanmat) <- names(coefs)
  ses <- apply(stanmat, 2L, mad)


  coefs <- split_naming(x = coefs, names_list = names_lst, N = N, K = K)
  ses <- split_naming(x = ses, names_list = names_lst, N = N, K = K)
  rownames(stan_summary) <- naming_fun(x = rownames(stan_summary), para_names = ses)

  covmat <- cov(stanmat)
  # rownames(covmat) <- colnames(covmat) <- rownames(stan_summary)[1:nrow(covmat)]
  if (object$algorithm == "sampling")
    check_rhats(stan_summary[, "Rhat"])

  # linear predictor, fitted values
  fit_extract <- rstan::extract(stanfit)
  yrep <- fit_extract$yrep
  res <- y - yrep

  # coefs

  out <- list(
    coefs_median = coefs,
    ses = ses,
    yrep = yrep,
    residuals = res,
    covmat,
    model = object$model,
    data = data_lst,
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
