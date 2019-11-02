#'
#'
#'
#'
#' @param formula_discrete The formula of the discrete Markov process
#' @param formula_continuous The formula of the continuous process

stan_msm <- function(formula_discrete = NULL, formula_continuous, family = gaussian(),
                     data = data, n = n, t = t, K = NULL, na.action = NULL,
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

    formula_discrete <- deparse_call_formula(mf$formula_discrete)
    formula_continuous <- deparse_call_formula(mf$formula_continuous)
    parsed_formula <- formula_parse(formula_discrete = formula_discrete,
                                    formula_continuous = formula_continuous)

    data <- data_split(data = data, tvtp = FALSE, parsed_formula = parsed_formula, n = n, t = t)
    x_d <- data$d
    x_e <- data$e
    y <- data$y
    n <- data$n
    t <- data$t
    has_intercept <- parsed_formula$has_intercept

    stanfit <- stan_msm.fit(x_e = x_e, x_d = x_d, y = y, n = n, t = t, K = K, has_intercept = has_intercept,
                            formula = parsed_formula, family = family, init.prior = init_prior,
                            algorithm = algorithm, iter = 1000, chains = 1)


    fit <- nlist(stanfit, algorithm, family, formula_discrete, formula_continuous,
                 data, y = y, x = list(x_d = x_d, x_e = x_e),
                 stan_function = "stan_msm",
                 model = mf, parsed_formula = parsed_formula, formula_discrete = formula_discrete, formula_continuous = formula_continuous,
                 call)

    out <- stanreg(fit)
    return(out)
}



