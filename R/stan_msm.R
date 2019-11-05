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
                     data = data, n = n, t = t, K = NULL, shared_TP = TRUE, order_continuous = c(),
                     na.action = NULL,
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

    call <- match.call(expand.dots = TRUE)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula_discrete", "formula_continuous"),
               table = names(mf), nomatch = 0L)
    mf <- mf[c(1L, m)]
    mf$data <- data

    # parse formula
    formula_discrete <- deparse_call_formula(mf$formula_discrete)     # extract the discrete formula as a string
    formula_continuous <- deparse_call_formula(mf$formula_continuous) # extract the continuous formula as a string
    parsed_formula <- formula_parse(formula_discrete = formula_discrete,
                                    formula_continuous = formula_continuous) # parse the formula

    # check data for inclusion
    check_data(data = data, order_continuous = order_continuous, parsed_formula = parsed_formula, n_var = n, t_var = t)


    data <- data_split(data = data, tvtp = FALSE, parsed_formula = parsed_formula, n = n, t = t)
    x_d <- data$d
    x_e <- data$e
    y <- data$y
    n <- data$n
    t <- data$t
    has_intercept <- parsed_formula$has_intercept
    N <- max(data$n)

    stanfit <- stan_msm.fit(x_e = x_e, x_d = x_d, y = y, n = n, t = t, K = K, has_intercept = has_intercept,
                            shared_TP = shared_TP, order_continuous = order_continuous,
                            formula = parsed_formula, family = family, init.prior = init_prior,
                            algorithm = algorithm, iter = 1000, chains = 1)


    fit <- list(stanfit = stanfit, algorithm = algorithm, family = family,
                 data = data, y = y, x = list(x_d = x_d, x_e = x_e), N = N, K = K,
                 stan_function = "stan_msm", model = mf, parsed_formula = parsed_formula,
                 formula_discrete = formula_discrete, formula_continuous = formula_continuous,
                 order_continuous = order_continuous, shared_TP = shared_TP,
                 call = call)

    out <- stanreg(fit)
    return(out)
}



