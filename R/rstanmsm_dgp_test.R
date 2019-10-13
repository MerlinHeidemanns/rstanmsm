###################################
# rstanmsm
# Test DGP
# Merlin Heidemanns
###################################



dgp_nt_sim <- function(N.sim = 1, N.seq = 1, T.seq = 300, n_chains = 2, n_iter = 2000){

  # Start time
  orig_time <- Sys.time()
  file <- file("progress.txt")

  # Output container
  out_par <- matrix(ncol = 12, nrow = 0)
  colnames(out_par) <- c("sim", "N", "T", "name", "mean", "sd", "q025", "q50", "q975", "switches", "t.1", "t.2")

  # Total number of simulations
  total.sim <- length(T.seq) * N.sim * length(N.seq)

  # parameter
  gamma <- rnorm(2, 0, 0.75); mu <- rnorm(2, 0, 4); lambda <- rnorm(2, 0.5, 0.25); phi <- rnorm(2, 0.5, 0.25)
  n_parameters <- length(gamma) + length(mu) + 2 + length(lambda) + length(phi)

  ## Initializing, chains, iter
  initf2 <- function(chain_id = 1) {
  # cat("chain_id =", chain_id, "\n")
  list(mu = sort(rnorm(2, 0, 4)), gamma = rnorm(2, 0, 0.75), pi1 = c(0.5, 0.5),
       lambda = rnorm(2, 0.5, 0.25), phi = sort(rnorm(2, 0.5, 0.25)), sigma = 1)
  }
  init_ll <- lapply(1:n_chains, function(id) initf2(chain_id = id))


  # Simulation
  sim_counter <- 0
  for (n1 in N.seq){
    for (t1 in T.seq){
      for (i in 1:N.sim){
        sim_counter <- sim_counter + 1
        data <- dgp_nt(N = n1, T = t1, gamma = gamma, mu = mu, lambda = lambda, phi = phi)
        fit <- stan("MSM_TVTP_sharedS_varyingNT.stan", data = data[[1]], chains = n_chains,
                    iter = n_iter,
                    init = init_ll)
        fit.ext <- as.matrix(fit)
        for (ii in 1:n_parameters){
          i_mean <- mean(fit.ext[, ii])
          i_sd <- sd(fit.ext[, ii])
          i_quan <- quantile(fit.ext[, ii], probs = c(0.025, 0.50, 0.975), names = F)
          out_par <- rbind(out_par, c(sim_counter, n1, t1, colnames(fit.ext)[ii], i_mean, i_sd, i_quan, data[[3]]$switches, data[[3]]$t.1, data[[3]]$t.2 ))
        }
        print(paste0("I've done ", sim_counter, " of ", total.sim, " simulations so far. ","The last size of T was ",
                            t1 , ". The last size of N was ", n1, ". The current date and time are ", as.character(Sys.time()), " and I started at ", as.character(orig_time), "."))
      }
      save(out_par, file = "out_tmp.RData")
    }
  }


  ## Clean
  out_par <- as.data.frame(out_par)
  i_seq <- c(1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12)
  for (i in i_seq){
    out_par[, i] <- as.numeric(as.character(out_par[, i]))
  }
  return(out_par)
}

## save output
save(out_par, file = "out_final.RData")

