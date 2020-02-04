#' msm_dgp
#'
#' @param N Number of units
#' @param T Number of average timepoints
#' @param K Number of states
#' @param S Number of maximum state processes
#' @param N_tp Number of transition probabilities
#' @param p Probability that an observation is missing.

msm_dgp <- function(N, T, K, S, N_tp, p){

  ## sample state processes
  ns <- sample(1:S, N, replace = TRUE)
  ns <- match(ns, unique(ns))
  S <- max(ns)

  ## assign state processes
  id_tp <- sample(1:N_tp, S, replace = TRUE)
  id_tp <- match(id_tp, unique(id_tp))
  NTP <- max(id_tp)

  ## sample length of time interval
  Tn <- rpois(N, T)
  NT <- sum(Tn)
  nts <- matrix(NA, ncol = 3, nrow = NT)
  for (n in 1:N){
    if (n == 1){
      nts[1:Tn[n], 1] <- n
      nts[1:Tn[n], 2] <- seq(1, Tn[n])
      nts[1:Tn[n], 3] <- ns[n]
    } else {
      nts[(sum(Tn[1:n - 1]) + 1):(sum(Tn[1:n - 1]) + Tn[n]), 1] <- n
      nts[(sum(Tn[1:n - 1]) + 1):(sum(Tn[1:n - 1]) + Tn[n]), 2] <- seq(1, Tn[n]) + rnbinom(1, 1, 0.7)
      nts[(sum(Tn[1:n - 1]) + 1):(sum(Tn[1:n - 1]) + Tn[n]), 3] <- ns[n]
    }
  }
  nts <- as_tibble(nts)
  nts <- rename(data, n = V1, t = V2, s = V3)
  nts <- mutate(data, q = id_tp[s])
  lenT <- group_by(nts, s)
  lenT <- summarize(nts, min.t = min(t), max.t = max(t))
  lenT <- mutate(nts, lenT = max.t - min.t)
  t.seq <- matrix(NA, ncol = S, nrow = max(lenT$max.t))

  ## coefficients
  mu <- sort(rnorm(K, 0, 5))
  sigma <- abs(rnorm(1, 0, 1))

  ## transition probabilities
  A <- array(NA, dim = c(NTP, K, K))
  for (n in 1:NTP){
    A[n, , ] <- rdirichlet(K, rep(1, K))
  }

  ## state processes
  s <- matrix(NA, ncol = S, nrow = max(lenT$max.t))
  seq.K <- seq(1, K, 1)
  for (i in 1:S){
    if (i %in% lenT$s){
      min.t <- as.integer(lenT[lenT$s == i, "min.t"])
      max.t <- as.integer(lenT[lenT$s == i, "max.t"])
      s[min.t, i] <- colSums(rmultinom(1, 1, prob = rep(1/K, K)) * seq.K)
      for (iter in (min.t + 1):max.t){
        s[iter, i] <- colSums(rmultinom(1, 1, prob = A[id_tp[i], s[iter - 1, i], ]) * seq.K)
      }
    }
  }

  # outcomes
  y <- rep(NA, times = NT)
  for (n in 1:N){
    idx <- which(nts$n == n)
    for (idx_id in idx){
      s_id <-
      if (idx_id == idx[1]){
        y[idx_id] <- rnorm(1, mu[s[nts$t[idx_id], nts$s[idx_id]]], sigma)
      } else {
        y[idx_id] <- rnorm(1, mu[s[nts$t[idx_id], nts$s[idx_id]]], sigma)
      }
    }
  }


  nts <- add_column(nts, y = y)
  nts <- dplyr::arrange(nts, q, s, t)

  start.stop.t <- matrix(NA, ncol = 2 * S, nrow = max(lenT$max.t))
  for (s.id in 1:S){
    for (t in 1:nrow(start.stop.t)){
      range <- which(nts$t == t & nts$s == s.id)
      if (length(range) != 0){
        start.stop.t[t, (s.id * 2) - 1] <- min(range)
        start.stop.t[t, (s.id * 2)] <- max(range)
      }
    }
  }
  slicer_T <- matrix(NA, ncol = 2, nrow = S)
  for (s.id in 1:S){
    range <- which(!is.na(s[, s.id]))
    slicer_T[s.id, ] <- c(min(range), max(range))
  }

  start.stop.t[is.na(start.stop.t)] <- 1e5


  pp1 <- pp2 <- pp3 <- 0
  Z <- 0
  z <- matrix(0, nrow = NTP * max(lenT$max.t), ncol = 0)
  pp_lambda <- matrix(0, nrow = 3, ncol = 0)
  pp_gamma <- array(0, dim = c(0))
  pp_eta <- array(0, dim = c(0))

  id_miss <- rbinom(NT * S, size = 1, p)

  # out
  stan_data <- list(N = N,
                    NS = S,
                    T = max(lenT$max.t),
                    NT = NT,
                    NTP = NTP,
                    id_tp = id_tp,
                    slicer_T = slicer_T,
                    startstop = start.stop.t,
                    pp1 = pp1,
                    pp2 = pp2,
                    pp3 = pp3,
                    Mz = Z,
                    z = z,
                    pp_lambda = pp_lambda,
                    pp_gamma = array(pp2, dim = c(pp2)),
                    pp_eta = array(pp3, dim = c(pp2)),
                    K = K,
                    Mx_d = 0,
                    Mx_e = 0,
                    z = z,
                    x_d = matrix(0, ncol = 0, nrow = NT),
                    x_e = matrix(0, ncol = 0, nrow = NT),
                    y = nts$y,
                    has_intercept = c(0, 1, 0, 0, 0),
                    shared_TP = 0,
                    state_sigma = 0,
                    tvtp = 0,
                    order_x_e = array(1, dim = 1),
                    id_miss = id_miss
  )

  data <- rename(nts, j = "s")
  data[which(id_miss == 1), "y"] <- NA

  parameters <- list(
    s = s,
    nts = nts,
    mu = mu,
    sigma = sigma,
    A = A
  )
  return(list(stan_data = stan_data, data = data, parameters = parameters))
}


