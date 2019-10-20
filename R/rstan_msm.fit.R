# rstanmsm
#
# rstan_msm.fit
# Call from rstan_msm


#' rstan_msm_fit: Call for rstan_msm
#'
#'
#' @export


rstan_msm.fit <- function(x_var, x_sha,
                         z_var, z_sha,
                         y = y, n, t, K = 2, family = gaussian(),
                         nchains = 4, niter = 2000, warmup = floor(niter/2),
                         thin = 1, controls = NULL,
                         init.prior = FALSE){

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
    names_x_var <- colnames(x_var)
    names_x_sha <- colnames(x_sha)
    names_z_var <- colnames(z_var)
    names_z_sha <- colnames(z_sha)

    # dimensions
    Md_var <- ncol(z_var)
    Md_sha <- ncol(z_sha)
    Mc_var <- ncol(x_var)
    Mc_sha <- ncol(x_sha)

    # data
    standata <- list(
       N = N,
       T = T,
       NT = NT,
       startstop = start.stop,
       K = K,
       Md_var = Md_var,
       Md_sha = Md_sha,
       Mc_var = Mc_var,
       Mc_sha = Mc_sha,
       z_sha = z_sha,
       z_var = z_var,
       x_sha = x_sha,
       x_var = x_var,
       y = y
     )

    # initialization
    if (init.prior) {
      mu_prior_mean <- gamma_prior_mean <- lambda_prior_mean <- zeta_prior_mean <- beta_prior_mean <- 0
      phi_prior_mean <- 0.5; phi_prior_sd <- 0.25
      mu_prior_sd <- 4; gamma_prior_sd <- 0.5; lambda_prior_sd <- 0.5; zeta_prior_sd <- 4; beta_prior_sd <- 4

      mu_init <- sort(rnorm(K, mu_prior_mean, mu_prior_sd), decreasing = FALSE)
      phi_init <- sort(rnorm(K, phi_prior_mean, phi_prior_sd), decreasing = FALSE)
      gamma_init <- rnorm(Md_sha, gamma_prior_mean, gamma_prior_sd)
      lambda_init <- matrix(rnorm(K * Md_var, lambda_prior_mean, lambda_prior_sd), ncol = K, nrow = Md_var)
      zeta_init <- rnorm(Mc_sha, zeta_prior_mean, zeta_prior_sd)
      beta_init <- matrix(rnorm(K * Mc_var, beta_prior_mean, beta_prior_sd), ncol = K, nrow = Mc_var)
      pi1_init <- matrix(rep(0.5, times = N * K), ncol = K, nrow = N)

      initf2 <- function(chain_id = 1) {
      list(pi1 = pi1_init,
           gamma = gamma_init,
           lambda = lambda_init,
           mu = mu_init,
           phi = phi_init,
           zeta = zeta_init,
           beta = beta_init,
           sigma = 1)
      }
      init_ll <- lapply(1:n_chains, function(id) initf2(chain_id = id))
      init = init_ll
    } else {
      init = "random"
    }

  fit <- rstan::sampling(stanmodels$rstan_msm_fit, data = standata,
                         chains = nchains, warmup = warmup, iter = 2000,
                         init = init)
  return(fit)
}
