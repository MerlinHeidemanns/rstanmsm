#' stan_msm.fit
#'
#' @param data Data is a list including matrixes for the sorted predictors and indicator vectors for n and t
#' @param priors A list object containing the prior for the Dirichlet and all other parameters.
#' @param order_continuous text
#' @param state_sigma text


stan_msm.fit <- function(data = data, state_sigma = FALSE, order_continuous = order_continuous,
                         priors = priors, family = gaussian(), init_prior = FALSE,
                         algorithm = c("optimizing", "sampling"), ...){

    # family
    family <- validate_family(family)
    supported_families <- c("gaussian")
    fam <- which(pmatch(supported_families, family$family, nomatch = 0L) == 1L)
    famname <- supported_families[fam]
    is_gaussian <- is.gaussian(famname)

    # standata
    standata <- prepare_standata(data. = data, priors. = priors,
                             order_continuous. = order_continuous,
                             state_sigma. = state_sigma)

    # stanfit
    stanfit <- stanmodels$msm_constant_continuous

    # parameters to exclude
    pars <- pars_include(pp1 = standata$pp1, pp2 = standata$pp2, pp3 = standata$pp3,
                         Mx_d = standata$Mx_d, Mx_e = standata$Mx_e)

    # initialization
    init <- "random"

    # execution
    if (algorithm == "optimizing") {
      optimizing_args <- list(...)
      if (is.null(optimizing_args$draws)) optimizing_args$draws <- 1000L
      optimizing_args$object <- stanfit
      optimizing_args$data <- standata
      optimizing_args$init <- init
      if (is.null(optimizing_args$tol_rel_grad)) optimizing_args$tol_rel_grad <- 10000L
      out <- do.call(optimizing, args = optimizing_args)
      out$stanfit <- suppressMessages(sampling(stanfit, data = standata, chains = 0))
    } else {
      if (algorithm == "sampling") {
          sampling_args <- set_sampling_args(
            object = stanfit,
            pars = pars,
            include = FALSE,
            user_dots = list(...),
            user_adapt_delta = adapt_delta,
            data = standata,
            show_messages = FALSE)
          stanfit <- do.call(rstan::sampling, sampling_args)
      }
      check <- try(check_stanfit(stanfit))
      if (!isTRUE(check)) return(standata)

      return(structure(stanfit))
    }
}
