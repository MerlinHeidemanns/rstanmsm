#' prior_mat
#'
#' @param prior: A vector for the priors
#'
#' @description If

prior_mat <- function(prior = prior, K = K, outcome = outcome){
  if (is.null(prior)){
    y_sd <- sd(y)
    priors <- matrix(NA, ncol = 4, nrow = 7)
    # discrete predictors
    priors[1:3,] <- cbind(rep(1,3), rep(0,3), rep(1, 3), 0)
    # linear predictors
    priors[4:6,] <- cbind(rep(1,3), rep(0,3), rep(y_sd, 3), 0)
    # variance
    priors[7,] <- c(1, 0, y_sd, 0)
    A_prior <- rep(1, K)
  }
}
