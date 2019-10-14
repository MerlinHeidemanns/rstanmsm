###################################
# rstanmsm
# Validate data
# Merlin Heidemanns
###################################

# Only reports whether the variables are there.
#

data_validate <- function(formula, data){
  tmp.form <- formula  # assign to tmp objects
  tmp.form <- strsplit(tmp.form, split = "[^a-z1-9A-Z_\\.]+")  #split apart
  tmp.var.con <- c()  # empty container
  for (i in tmp.form){
    tmp.var.con <- c(tmp.var.con, i) #  store all in container
  }
  tmp.var.con <- unique(tmp.var.con)  # unique elements
  all.in.data <- ifelse(all(tmp.var.con %in% colnames(data)), TRUE, FALSE)  # test whether in data
  return(all.in.data)
}
