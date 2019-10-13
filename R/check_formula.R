###################################
# rstanmsm
# Check formula
# Merlin Heidemanns
###################################

# Don't expose

check.formula <- function(formula, K, data = FALSE){
  tmp.form <- formula
  tmp.form <- strsplit(tmp.form, split = "[^a-z1-9A-Z_\\.]+")  #splits on anything but characters usable for variables
  if (data == TRUE){

  }
}
