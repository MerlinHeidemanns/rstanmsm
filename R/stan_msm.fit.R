



stan_msm.fit <- function(x_e, x_d, y = y, n, t, K = 2, has_intercept = rep(0, 5),
                         shared_TP = TRUE, order_continuous,
                         formula = parsed_formula, family = gaussian(), init.prior = TRUE,
                         algorithm = c("optimizing", "sampling"), ... = ...){

    # NT
    N <- max(n)
    T <- max(t)
    NT <- N * T

    # Coerce to matrixes
    x_e <- as.matrix(x_e)
    x_d <- as.matrix(x_d)

    # dimensions
    Mx_d <- ncol(x_d)
    Mx_e <- ncol(x_e)

    if (Mx_e == 0){
      x_e <- matrix(0, ncol = 0, nrow = NT)
    }
    if (Mx_d == 0){
      x_d <- matrix(0, ncol = 0, nrow = NT)
    }

    # order vector
    order_x_e <- create_order_vector(formula, order_continuous)

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
    if (has_intercept[1] == 1) names_list$alpha <- c("Intercept", names_list$alpha)
    if (has_intercept[2] == 1) names_list$beta <- c("Intercept", names_list$beta)


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
       has_intercept = has_intercept,
       shared_TP = shared_TP,
       order_x_e = order_x_e
     )

    # stanfit
    stanfit <- stanmodels$msm_constant_continuous_v2

    # parameters to exclude
    pars <- pars_include(Mx_d = Mx_d, Mx_e = Mx_e)

    # initialization
    N_ <- if (shared_TP) 1 else N
    Mx_d_ <- if (has_intercept[1] == 1) Mx_d + 1 else if (Mx_d == 0) 1 else Mx_d
    Mx_e_ <- if (has_intercept[2] == 1) Mx_e + 1 else if (Mx_e == 0) 1 else Mx_e

    if (init.prior) {
      initf2 <- function(chain_id = 1) {
        list(
          pi1 = matrix(rep(1/K, N * K), ncol = K, nrow = N),
          A = array(rdirichlet(K * N_, rep(1, K)), dim = c(N_, K, K)),
          phi = matrix(sort(rnorm(K, 0.5, 0.25), decreasing = FALSE), ncol = K, nrow = 1),
          mu = matrix(sort(rnorm(K, 0, 5)), ncol = K, nrow = 1),
          alpha = array(rnorm(Mx_d_, 0, 1)),
          beta = matrix(rnorm(K * Mx_e_, 0, 1), ncol = K, nrow = Mx_e_),
          sigma = 1)
        }
      init_ll <- lapply(1:n_chains, function(id) initf2(chain_id = id))
      init = init_ll
    } else {
      init = "random"
    }

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
