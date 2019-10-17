###################################
# rstanmsm
# Data Generating Processes
# Merlin Heidemanns
###################################

#' DGP for MSM-TVTP balanced panel
#'
#' Generate data for a MSM-TVTP balanced panel
#'
#' \code{dgp_nt} simulates data to test the model. It outputs the observable data,
#' the parameters used to generate the data, and summary statistics. The summary statistics
#' include the number of switches between states and the time spent in each state.
#' @param N The number of separate time series.
#' @param T The length of each time series.
#' @param gamma A length 2 vector of intercepts for the TVTP.
#' @param lambda A length 2 vector of coefficients on the predictor on z for the TVTP.
#' @param mu A length 2 vector of intercepts for the linear models.
#' @param phi A length 2 vector of AR1 coefficients. Each should be between 0 and 1.
#' @return dgp_nt() returns three lists. The first lists contains the observed data.
#' The second list contains the parameters.
#' The third list contains summary statistics.
#' \describe{
#'   \item{switches}{The total number of switches of all discrete processes.}
#'   \item{t.1}{The total time spent in state 1.}
#'   \item{t.2}{The total time spent in state 1.}
#' }
#' @keywords DGP
#' @examples
#' dgp_nt(N = 2, T = 300, gamma = c(0.25, 0.5), lambda = c(0.75, 0.5), mu = c(2, 6), phi = c(0.2, 0.8))


dgp_nt_k2 <- function(N, T, Mx_var, Mx_sha, Mz_var, Mz_sha,
                      Mx_var_int = FALSE, Mx_sha_int = FALSE,
                      Mz_var_int = FALSE, Mz_sha_int = FALSE,
                      para.var = FALSE,
                      gamma = NULL, lambda = NULL, zeta = NULL, beta = NULL,
                      fit = FALSE){
  K <- 2

  # parameters
  mu <- sort(rnorm(K, 0, 3))
  phi <- rnorm(K, 0, 0.5)
  if (para.var){
    gamma <- matrix(rnorm(Mz_sha , 0, 1), ncol = Mz_sha, nrow = 1)
    lambda <- matrix(rnorm(Mz_var * K, 0, 1), ncol = Mz_var, nrow = K)
    zeta <- matrix(rnorm(Mx_sha, 0, 1), ncol = Mx_var, nrow = 1)
    beta <- matrix(rnorm(Mz_var * K, 0, 1), ncol = Mx_sha, nrow = K)
    sigma <- 1
  } else if (para.var == FALSE){
    if (is.null(gamma)){
      stop("Assign gamma[Mz_sha] parameters or set para.var to TRUE.")
    }
    if (is.null(lambda)){
      stop("Assign lambda[K, Mz_var] parameters or set para.var to TRUE.")
    }
    if (is.null(zeta)){
      stop("Assign zeta[Mx_sha] parameters or set para.var to TRUE.")
    }
    if (is.null(beta)){
      stop("Assign beta[K, Mz_var] parameters or set para.var to TRUE.")
    }
  }


  ## slicer
  start.stop <- matrix(NA, ncol = 2, nrow = N)
  for (n in 1:N){
    start.stop[n, ] <- c(1 + (n - 1) * T, (n - 1) * T + T)
  }

  # predictors
  z_var <- matrix(rbinom(N * T * Mz_var, 1, 0.5), ncol = Mz_var, nrow = N * T)
  z_sha <- matrix(rbinom(N * T * Mz_sha, 1, 0.5), ncol = Mz_sha, nrow = N * T)
  x_var <- matrix(rbinom(N * T * Mx_var, 1, 0.5), ncol = Mx_var, nrow = N * T)
  x_sha <- matrix(rbinom(N * T * Mx_sha, 1, 0.5), ncol = Mx_sha, nrow = N * T)

  # add intercepts
  if (Mx_var_int) {x_var <- cbind(1, x_var)}
  if (Mx_sha_int) {x_sha <- cbind(1, x_sha)}
  if (Mz_var_int) {z_var <- cbind(1, z_var)}
  if (Mz_sha_int) {z_sha <- cbind(1, z_sha)}

  # tvtp
  A <- array(NA, dim = c(2, 2, T, N))
  for (n in 1:N){
    for (t in 1:T){
      A[1,1,t, n] <- pnorm(t(z_sha[start.stop[n, 1] + t - 1, ]) %*% gamma +
                           t(z_var[start.stop[n, 1] + t - 1, ]) %*% lambda[1, ])
      A[1,2,t, n] <- 1- A[1,1,t, n]
      A[2,2,t, n] <- pnorm(t(z_sha[start.stop[n, 1] +  t - 1]) %*% gamma +
                           t(z_var[start.stop[n, 1] + t - 1, ]) %*% lambda[2, ])
      A[2,1,t, n] <- 1 - A[2,2,t, n]
    }
  }

  # states
  s <- matrix(NA, nrow = T, ncol = N)
  s[1, ] <- rbinom(N, 1, 0.5) + 1
  for (n in 1:N){
    for (t in 2:T){
      s[t, n] <- rbinom(1, 1, A[s[t - 1, n], 2, t, n]) + 1
    }
  }

  # outcomes
  y <- rep(NA, times = N * T)
  for (n in 1:N){
    # t == 1
    y[start.stop[n, 1]] <- rnorm(1, mu[s[1, n]], sigma)

    # t >= 2
    for (t in 2:T){
      y[start.stop[n, 1] + t - 1] <- rnorm(1, mu[s[t, n]] +
                                              phi[s[t, n]] * y[start.stop[n, 1] + t - 2],
                                           sigma)
    }
  }

  data <- list(N = N,
              T = T,
              y = y,
              z = z)
  para <- list(pi1 = c(0.5, 0.5),
               gamma = gamma,
               mu = mu,
               lambda = lambda,
               phi = phi)
  other <- list(switches = sum(abs(s[2:T] - s[1:T-1])),
                t.1 = length(s[s == 1]),
                t.2 = length(s[s == 2]))
  out <- list(data, para, other)
  return(out)
}

