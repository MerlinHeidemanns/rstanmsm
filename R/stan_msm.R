#' stan_msm
#'
#'
#'
#' @param formula_discrete The formula of the discrete Markov process
#' @param formula_continuous The formula of the continuous process
#' @param order_continuous A character vector containing the names of the predictors that are supposed to be ordered.
#' It can include "Intercept" to indicate that the intercept is supposed to be ordered.
#' @param data A dataframe containing the variables references in `formula_discrete` and `formula_continuous`.
#' @param family The distribution of the outcome.
#' @param j State-process association
#' @param q Transition-process association
#'
#' @export

stan_msm <- function(formula_discrete = NULL, formula_continuous, family = gaussian(),
                     data = data, n = NULL, t = NULL, j = NULL, q = NULL, K = NULL, state_sigma = FALSE,
                     state_varying_continuous = c(), state_varying_discrete = list(state = NULL, state_state = NULL),
                     order_continuous = c(), na.action = NULL,
                     ...,
                     prior = c(alpha = normal(), beta = normal(), gamma = normal(),
                               delta = normal(), eta = normal(), sigma()),
                     algorithm = c("sampling", "optimizing"),
                     init_prior = FALSE,
                     adapt_delta = NULL) {

    algorithm <- match.arg(algorithm)
    family <- validate_family(family)
    #check_tp_s(shared_TP = shared_TP, shared_S = shared_S, n = n) # check for combinations that are not allowed.

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
               n = n, t = t,
               state_varying_continuous = state_varying_continuous, state_varying_discrete = state_varying_discrete)

    # parse
    parsed_data_names <- data_parse(formula_continuous = formula_continuous,
                                    formula_discrete = formula_discrete,
                                    state_varying_continuous = state_varying_continuous,
                                    state_varying_discrete = state_varying_discrete,
                                    data = data, n = n, t = t, j = j, q = q, K = K)

    # priors
    priors <- prior_mat(prior = prior, K = K, y = parsed_data_names[["data_lst"]]$y)

    stanfit <- stan_msm.fit(data = parsed_data_names[["data_lst"]], K = K,
                            state_sigma = state_sigma,
                            order_continuous = order_continuous,
                            family = family, init_prior = init_prior,
                            id_miss = parsed_data_names[["id_miss"]],
                            priors = priors,
                            algorithm = algorithm, ...)

    fit <- list(stanfit = stanfit, algorithm = algorithm, family = family,
                 data = parsed_data_names, stan_function = "stan_msm",
                 formula_discrete = formula_discrete, formula_continuous = formula_continuous,
                 order_continuous = order_continuous, state_sigma = state_sigma,
                 call = call)

    out <- stanreg(fit)
    return(out)
}



