#' stan_msm.fit
#'
#' @param data Data is a list including matrixes for the sorted predictors and indicator vectors for n and t
#' @param K The number of states
#' @param shared_TP
#' @param shared_state
#' @param state_sigma


stan_msm.fit <- function(data, K = 2, shared_TP = TRUE, shared_state = FALSE, state_sigma = FALSE, order_continuous,
                          family = gaussian(), init_prior = TRUE,
                          algorithm = c("optimizing", "sampling"), ... = ...){

    # NT
    N <- data$N
    n <- data$n
    T <- data$T
    NT <- data$NT

    # model specifications
    K <- data$K
    has_intercept <- data$has_intercept
    id_miss <- data$id_miss

    # Coerce to matrixes
    x_d <- data$x_d
    x_e <- data$x_e
    z <- cbind(data$x_a, data$x_b, data$x_c)

    # dimensions
    Mx_d <- ncol(x_d)
    Mx_e <- ncol(x_e)
    Mz <- ncol(z)

    # fix if NULL
    if (Mx_d == 0) x_d <- matrix(0, ncol = 0, nrow = NT)
    if (Mx_e == 0) x_e <- matrix(0, ncol = 0, nrow = NT)
    if (Mx_z == 0) z <- matrix(0, ncol = 0, nrow = NT)

    # pp


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

    # state process
    NS <- data$J
    NTP <- data$Q
    id_tp <- data$id_tp

    # slicer
    slicer_T <- slicer_time(data$j)
    start_stop <- start_stop_slicer(j = data$j, t = data$t)

    # standata
    standata <- list(
      N = N, # units
      T = T, # max of time points
      NT = NT, # N of observations
      NS = NS, # N of state processes
      NTP = NTP, # N of transition probabilities
      id_tp = id_tp,      # [NS] which state process belongs to which tp matrix
      slicer_T = slicer_T,   # min(N, T):max(N, T)
      startstop = start_stop,  # slicer for s units
      K = K,         # N of states
      has_intercept = has_intercept, # [5] 1: alpha, 2: beta; 3: gamma; 4: delta; 5: eta
      Mz = Mz, # N of tp predictors
      Mx_d = Mx_d, # N of fixed parameters of continuous process
      Mx_e = Mx_e, # N of continuous parameters of continuous process
      pp1 = pp1, pp2 = pp2, pp3 = pp3, # N of general, state, and state-state specific predictors
      pp_lambda = , # [3, Mz] which are varying at which level for tps
      pp_gamma = , # [pp2] 0/1 of varying at state
      pp_eta = , # [pp2]   0/1 of varying at state-state
      z = z,      # [NTP * T, Mz] matrix of predictors of discrete process
      x_d = x_d,    # [NT, Mx_d] matrix of fixed predictors of continuous process
      x_e = x_e,    # [NT, Mx_e] matrix of varying predictors of continuous process
      y = y,      # [NT] vector of output
      state_sigma = state_sigma, # 0: general,    1: state-specific
      tvtp = ,
      order_x_e = order_x_e, # [Mx_e + has_intercept[2]] 0: unordered, 1: ordered
      A_prior = , # [K] A_prior;
      priors = ,  #[7,4];     // 1: Kind, 2: mean, 3: sd, 4: df; 1: normal, 2: cauchy, 3: student-t
      id_miss = id_miss  # 1: missing, 0: present / at least one observation
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
