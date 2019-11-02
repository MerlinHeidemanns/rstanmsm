



stan_msm.fit <- function(x_e, x_d, y = y, n, t, K = 2, has_intercept = c(0, 0),
                         formula = parsed_formula, family = gaussian(), init.prior = TRUE,
                         algorithm = c("optimizing", "sampling"), ... = ...){

    # Coerce to matrixes
    x_e <- as.matrix(x_e)
    x_d <- as.matrix(x_d)

    # NT
    N <- max(n)
    T <- max(t)
    NT <- N * T

    # family
    family <- validate_family(family)
    supported_families <- c("gaussian")
    fam <- which(pmatch(supported_families, family$family, nomatch = 0L) == 1L)
    famname <- supported_families[fam]
    is_gaussian <- is.gaussian(famname)

    # start.stop
    start.stop <- matrix(NA, ncol = 2, nrow = N)
    for (nn in 1:N){
      start.stop[nn, 1] <- min(seq(1:NT)[n == nn])
      start.stop[nn, 2] <- max(seq(1:NT)[n == nn])
    }

    # names
    names_list = list(
      alpha = formula$d,
      beta = formula$e
    )

    # dimensions
    Mx_d <- ncol(x_d)
    Mx_e <- ncol(x_e)

    # standata
    standata <- list(
       N = N,
       T = T,
       NT = NT,
       startstop = start.stop,
       K = K,
       Mx_e = Mx_e,
       Mx_d = Mx_d,
       x_d = x_d,
       x_e = x_e,
       y = y,
       has_intercept = has_intercept
     )

    # stanfit
    stanfit <- m # EXCHANGE

    # initialization
    if (init.prior) {
      initf2 <- function(chain_id = 1) {
        list(
          pi1 = matrix(rep(1/K, N * K), ncol = K, nrow = N),
          lambda = rdirichlet(K, rep(1, K)),
          phi = matrix(sort(rnorm(K, 0.5, 0.25), decreasing = FALSE), ncol = K, nrow = 1),
          mu = matrix(sort(rnorm(K, 0, 5)), ncol = K, nrow = 1),
          alpha = rnorm(Mx_d_, 0, 1),
          beta = matrix(rnorm(K * Mx_e_, 0, 1), ncol = K, nrow = Mx_e_),
          sigma = 1)
        }
      init_ll <- lapply(1:n_chains, function(id) initf2(chain_id = id))
      init = init_ll
    } else {
      init = "random"
    }

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
            #prior = prior,
            user_dots = list(...),
            user_adapt_delta = adapt_delta,
            data = standata,
            show_messages = FALSE)
          stanfit <- do.call(rstan::sampling, sampling_args)
      }
      check <- try(check_stanfit(stanfit))
      if (!isTRUE(check)) return(standata)

      # naming
      #stanfit@sim$fnames_oi <- naming_fun(N = N,
      #                                    T = T,
      #                                    K = K,
      #                                    names = stanfit@sim$fnames_oi,
      #                                    names_list = names_list)

      return(structure(stanfit))
    }
}
