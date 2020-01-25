#' prior_mat
#'
#' @param prior: A list of the priors.
#' @param K: The number of states
#' @param y: The outcome
#'
#' @details If priors are null, then everything becomes normal.
#'
#' @return Returns a list with the transition probability matrix Dirichlet parameters and a matrix of the priors for the rest.

prior_mat <- function(prior = prior, K = K, y = y){
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
  return(list(A_prior = A_prior, priors = priors))
}
