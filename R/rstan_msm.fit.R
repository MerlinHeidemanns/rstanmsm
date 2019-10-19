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
                         y = y, n, t, K = 2, family = gaussian()){

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
  fit <- rstan::sampling(stanmodels$rstan_msm_fit, data = standata)
  return(fit)
}
