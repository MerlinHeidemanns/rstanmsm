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


dgp_nt_istates <- function(N, T, gamma, mu, lambda, phi){
  K <- 2

  # tvtp
  z <- matrix(rbinom(N * T, 1, 0.5), ncol = N, nrow = T)
  A <- array(NA, dim = c(2, 2, T, N))
  for (t in 1:T){
    A[1,1,t, ] <- pnorm(gamma[1] + lambda[1] * z[t, ])
    A[1,2,t, ] <- 1- A[1,1,t, ]
    A[2,2,t, ] <- pnorm(gamma[2] + lambda[2] * z[t, ])
    A[2,1,t, ] <- 1 - A[2,2,t, ]
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
  y <- matrix(NA, nrow = T, ncol = N)

  # t == 1
  y[1, ] <- rnorm(N, mu[s[1,]], 1)

  # t >= 2
  for (t in 2:T){
    y[t, ] <- rnorm(N, mu[s[t, ]] +  phi[s[t, ]] * y[t - 1, ], 1)
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

