#' stan_msm
#'
#'
#'
#' @param formula_discrete The formula of the discrete Markov process
#' @param formula_continuous The formula of the continuous process
#' @order_continuous A character vector containing the names of the predictors that are supposed to be ordered.
#' It can include "Intercept" to indicate that the intercept is supposed to be ordered.
#' @export

stan_msm <- function(formula_discrete = NULL, formula_continuous, family = gaussian(),
                     data = data, n = n, t = t, K = NULL,
                     shared_TP = TRUE, shared_S = FALSE,
                     state_varying_continuous = c(), state_varying_discrete   = list(),
                     order_continuous = c(), na.action = NULL,
                     ... = ...,
                     prior = normal(),
                     prior_intercept = normal(),
                     prior_aux = exponential(),
                     prior_PD = FALSE,
                     algorithm = c("sampling", "optimizing"),
                     init_prior = FALSE,
                     adapt_delta = NULL) {


    algorithm <- match.arg(algorithm)
    family <- validate_family(family)
    check_tp_s(shared_TP = shared_TP, shared_S = shared_S, n = n) # check for combinations that are not allowed.

    call <- match.call(expand.dots = TRUE)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula_discrete", "formula_continuous"),
               table = names(mf), nomatch = 0L)
    mf <- mf[c(1L, m)]

    # parse formula
    formula_discrete <- deparse_call_formula(mf$formula_discrete)     # extract the discrete formula as a string
    formula_continuous <- deparse_call_formula(mf$formula_continuous) # extract the continuous formula as a string

    # check data for inclusion
    check_data(data = data, order_continuous = order_continuous,
               formula_continuous = formula_continuous, formula_discrete = formula_discrete,
               n_var = n, t_var = t,
               state_varying_continuous = state_varying_continuous, state_varying_discrete = state_varying_discrete)

    parsed_data_names <- data_parse(formula_continuous = formula_continuous,
                                    formula_discrete = formula_discrete,
                                    state_varying_continuous = state_varying_continuous,
                                    state_varying_discrete = state_varying_discrete,
                                    data = data, n_var = n, t_var = t, K = K)

    stanfit <- stan_msm.fit(data = parsed_data_names[["data_lst"]], K = K, shared_TP = shared_TP, shared_S = shared_S,
                            order_continuous = order_continuous,
                            family = family, init.prior = init_prior,
                            algorithm = algorithm, iter = 1000, chains = 1)


    fit <- list(stanfit = stanfit, algorithm = algorithm, family = family,
                 data = parsed_data_names, stan_function = "stan_msm",
                 formula_discrete = formula_discrete, formula_continuous = formula_continuous,
                 order_continuous = order_continuous, shared_TP = shared_TP, shared_S = shared_S,
                 call = call)

    out <- stanreg(fit)
    return(out)
}



