#' stan_msm.fit
#'
#' @param data Data is a list including matrixes for the sorted predictors and indicator vectors for n and t



stan_msm.fit <- function(data, K = 2, shared_TP = TRUE, shared_S = FALSE, order_continuous,
                          family = gaussian(), init.prior = TRUE,
                          algorithm = c("optimizing", "sampling"), ... = ...){

    # NT
    N <- data$N
    n <- data$n
    T <- data$T
    NT <- N * T

    # model specifications
    K <- data$K
    has_intercept <- data$has_intercept

    # Coerce to matrixes
    x_d <- data$x_d
    x_e <- data$x_e

    # dimensions
    Mx_d <- ncol(x_d)
    Mx_e <- ncol(x_e)

    # fix if NULL
    if (Mx_d == 0) x_d <- matrix(0, ncol = 0, nrow = NT)
    if (Mx_e == 0) x_e <- matrix(0, ncol = 0, nrow = NT)

    # order vector
    order_x_e <- create_order_vector(data, order_continuous)

    # output
    y <- data$y

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
       shared_S = shared_S,
       order_x_e = order_x_e
     )

    # stanfit
    stanfit <- stanmodels$msm_constant_continuous

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
