



coef.stan_msm <- function(object){
  coef_discrete <- c()
  coef_continuous <- c()
  coef_continuous <- rbind(fit$coefs_median$alpha, fit$coefs_median$beta)
  return(list(coef_discrete = coef_discrete, coef_continuous = coef_continuous))
}



#
# print.stanreg_msm <- function(x, digits = 1, detail = TRUE, ...) {
#   if (detail) {
#     cat(x$stan_function)
#     cat("\n Family:      ", family_plus_link(x))
#     cat("\n Formula:     ", formula_string(formula(x)))
#     cat("\n Observations:", nobs(x))
#     if (isTRUE(x$stan_function %in%
#                c("stan_msm"))) {
#       cat("\n Predictors:  ", length(coef(x)$))
#       cat("\n Time-varying Predictors:   ", )
#     } else if (isTRUE(x$stan_function %in%
#                c("stanmsm_tvtp")))
#
#     cat("\n------\n")
#   }
#
#   if (detail) {
#     cat("\n------\n")
#     cat("* For help interpreting the printed output see ?print.stanreg\n")
#     cat("* For info on the priors used see ?prior_summary.stanreg\n")
#   }
#   invisible(x)
# }
