###################################
# rstanmsm
# Main function
# Merlin Heidemanns
###################################

#'
#'
#'
#'
#' @param sd Indicates whether standard deviation of the distribution is shared or state-specific.



rstan_msm <- function(formula = c() , identical = FALSE, K = FALSE, data = FALSE, na.action = "PLACEHOLDER",
                      sd = c("shared")){

  call <- match.call(expand.dots = TRUE)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula"), table = names(mf), nomatch = 0L)
  mf <- mf[c(1L, m)]
  mf$data <- data


  if(data_validate(data) != TRUE){
    stop("Not all variables were found in the dataframe supplied.")
  }
  data_generate
  #if (K | data){
  #  stop("Please provide the number of states you intend to estimate.")
  #}
  #if (length(formula) != K & identical == FALSE){
  #  stop("Please provide a ")
  #}
}
