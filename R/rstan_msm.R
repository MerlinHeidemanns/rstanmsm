###################################
# rstanmsm
# Main function
# Merlin Heidemanns
###################################

#'
#'
#'
#'
#'



rstan_msm <- function(formula = c() , identical = FALSE, K = FALSE, data = FALSE){
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
