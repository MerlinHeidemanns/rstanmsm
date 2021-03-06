% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stan_msm.R
\name{stan_msm}
\alias{stan_msm}
\title{stan_msm}
\usage{
stan_msm(formula_discrete = NULL, formula_continuous,
  family = gaussian(), data = data, n = n, t = t, K = NULL,
  shared_TP = TRUE, shared_S = FALSE, state_varying_continuous = c(),
  state_varying_discrete = list(), order_continuous = c(),
  na.action = NULL, ... = ..., prior = normal(),
  prior_intercept = normal(), prior_aux = exponential(),
  prior_PD = FALSE, algorithm = c("sampling", "optimizing"),
  init_prior = FALSE, adapt_delta = NULL)
}
\arguments{
\item{formula_discrete}{The formula of the discrete Markov process}

\item{formula_continuous}{The formula of the continuous process}
}
\description{
stan_msm
}
