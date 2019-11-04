
#' formula_parse
#'
#' Parses the formulas into their components and stores them in a list

formula_parse <- function(formula_discrete = NULL, formula_continuous = formula_continuous) {

  int_pattern <- "(^1$)|(^1\\|2$)|(^1\\|3$)"

  # out list
  out.lst <- list(y = NULL,
                   a = NULL,
                   b = NULL,
                   c = NULL,
                   d = NULL,
                   e = NULL,
                   has_intercept = rep(0, 5), # a OR b and c OR d OR e
                   all.var = NULL)

  # formula continuous
  tmp.c <- formula_continuous
  tmp.c <- gsub("\\s", "", tmp.c)
  tmp.c <- unique(strsplit(tmp.c, split = "~|\\+" )[[1]])
  tmp.c.d <- c()
  tmp.c.e <- c()
  for (i in tmp.c[2:length(tmp.c)]){
    if (grepl(int_pattern, i)){
      if (grepl("\\|2", i)){
        out.lst$has_intercept[2] <- 1
      } else {
        out.lst$has_intercept[1] <- 1
      }
    } else if (grepl("\\|2", i)){
      tmp.var <- i
      tmp.var <- gsub("\\|2", "", tmp.var)
      tmp.c.e <- c(tmp.c.e, tmp.var)
    } else {
      tmp.c.d <- c(tmp.c.d, i)
    }
  }
  out.lst$all.var <- c(out.lst$all.var, tmp.c.d, tmp.c.e)
  out.lst$y <- tmp.c[1]; out.lst$d <- tmp.c.d; out.lst$e <- tmp.c.e
  # formula discrete
  if (!is.null(formula_discrete)){
    tmp.d <- formula_discrete
    tmp.d <- gsub("\\s", "", tmp.d)
    tmp.d <- unique(strsplit(tmp.d, split = "~|\\+" )[[1]])
    tmp.d.a <- NULL
    tmp.d.b <- NULL
    tmp.d.c <- NULL
    for (i in tmp.d){
      if (grepl(int_pattern, i)){
        if (grepl("\\|3", i)){
          out.lst$has_intercept[5] <- 1
        } else if (grepl("\\|2", i)){
          out.lst$has_intercept[4] <- 1
        } else {
          out.lst$has_intercept[3] <- 1
        }
      } else if (grepl("\\|3", i)){
        tmp.var <- i
        tmp.var <- gsub("\\|3", "", tmp.var)
        tmp.d.c <- c(tmp.d.c, tmp.var)
      } else if (grepl("\\|2", i)){
        tmp.var <- i
        tmp.var <- gsub("\\|2", "", tmp.var)
        tmp.d.b <- c(tmp.d.b, tmp.var)
      } else {
        tmp.d.a <- c(tmp.d.a, i)
      }
    }
    out.lst$a = tmp.d.a
    out.lst$b = tmp.d.b
    out.lst$c = tmp.d.c
    out.lst$all.var <- c(out.lst$all.var, tmp.d.a, tmp.d.b, tmp.d.c)
  }
  return(out.lst)
}



# data_check
#
# @data A dataframe object
# @par All predictors that are to be included in the model.
#
# Stops execution if a variable is not found in the dataframe.

data_check <- function(data = data, par = par) {
  dta.par <- colnames(data)
  for (i in par){
    if (!(i %in% dta.par)){
      stop(cat(i, "not found in supplied dataframe."))
    }
  }
}

# Data split
#
# @data A dataframe object.
# @tvtp TRUE/FALSE Indicator for whether the function includes time-varying transition probabilities.
# @parsed_formula A parsed formula, i.e. a list object supplied by formula_parse

data_split <- function(data = data, tvtp = FALSE, parsed_formula = parsed_formula, n = NULL, t = NULL) {

  # initialize output list
  out.lst <- list(a = NULL, b = NULL, c = NULL, d = NULL, e = NULL, y = NULL, n = NULL, t = NULL)
  if (!is.data.frame(data)){
    data <- as.data.frame(data)
  }
  names_data <- colnames(data)
  if (is.element(n, names_data) & is.element(t, names_data)){
    data <- data[order(data$n, data$t),]
  } else {
    stop("Please supply indexes for n and t.")
  }
  if (tvtp) {
    out.lst$a <- data[, parsed_formula$a]
    out.lst$b <- data[, parsed_formula$b]
    out.lst$c <- data[, parsed_formula$c]
  } else {
    out.lst$d <- data[, parsed_formula$d]
    out.lst$e <- data[, parsed_formula$e]
  }

  out.lst$y <- data[ , parsed_formula$y]
  out.lst$n <- data[ , n]
  out.lst$t <- data[ , t]

  return(out.lst)
}




#' check_order
#'
#' @param data Included data frame
#' @param order_continuous Vector of parameter names declared to be continuous

check_data <- function(data, order_continuous, parsed_formula, n_var, t_var){
  names_dta <- colnames(data)
  # check inclusion
  .check_inclusion(names_dta, parsed_formula = parsed_formula)
  # check for n and t
  .check_nt(names_dta = names_dta, n_var = n_var, t_var = t_var)
  # check for order predictor
  .check_order(names_dta = parsed_formula$e, order_continuous = order_continous, has_intercept = parsed_formula$has_intercept)
}


.check_order <- function(names_dta, order_continuous, has_intercept){
  cnd1 <- !(any(grepl(order_continuous, names_data) == TRUE))
  cnd2 <- !(grepl("Intercept", order_continuous) & (has_intercept[2] == 1))
  if (cnd1 & cnd2) {
    stop("Predictors that are supposed to be ordered are not indicated as varying by state.")
  }
}

.check_nt <- function(names_dta, n_var, t_var){
  if (!(is.element(n_var, names_dta) & is.element(t_var, names_dta))){
    stop("Please supply indexes for n and t.")
  }
}

.check_inclusion <- function(names_dta, parsed_formula){
  for (i in parsed_formula$all.var){
    if (!is.element(i, names_dta)) stop(paste0("The predictor ", i, " has not been found in the data."))
  }
}

#' create_order_vector
#'
#' @param formula The formula list
#' @param order_continuous A character vector indicating which parameters are ordered.
#'
#' This function creates a 0/1 vector indicating whether a particular parameter that is varying across states is ordered.

create_order_vector <- function(formula, order_continuous){
  tmp <- c("Intercept", order_continuous)
  out <- rep(0, length(tmp))
  out[match(formula$e, tmp)] <- 1
  return(out)
}


split_coef <- function(x, formula){
  # Initialize list
  out.lst <- list(coefs_a = NULL, coefs_b = NULL, coefs_c = NULL, coefs_d = NULL, coefs_e = NULL,
                  A = NULL, pi1 = NULL)

  # subset parameter vector on names
  out.lst$A <- x[grepl("A\\[.+\\]",names(x))]
  out.lst$pi1 <- x[grepl("pi1\\[.+\\]",names(x))]
  out.lst$coefs_a <- x[formula$a]
  out.lst$coefs_b <- x[formula$b]
  out.lst$coefs_c <- x[formula$c]
  out.lst$coefs_d <- x[formula$d]
  out.lst$coefs_e <- x[formula$e]

  # return
  return(out.lst)
}

#' state_name
#'
#'
#' @param x A character vector
#' @param K The number of states

naming_state <- function(names, K) {
  tmp <- c()
  for (j in names){
    for (k in seq(1, K, 1)){
      tmp <- c(tmp, paste0("S", k, "_", j))
    }
  }
  return(tmp)
}

#' split_naming
#'
#'
#' @param x A named vector of coefficients or standard errors
#'

split_naming <- function(x, names_list, N, K, shared_TP = TRUE){

  out.lst <- list(init_prob = NULL, tp = NULL, intercept = NULL, AR1 = NULL, alpha = NULL,
                  beta = NULL, gamma = NULL, delta = NULL, lambda = NULL)
  init_prob <- x[grepl("pi1\\[.+", names(x))]
  tp <- x[grepl("A\\[[0-9]+,[0-9]+,[0-9]+\\]", names(x))]

  # Intercept
  mu <- x[grepl("mu\\[.+\\]", names(x))]  # intercept
  names(mu) <- naming_state("Intercept", K)

  # AR1
  phi <- x[grepl("phi\\[.+\\]", names(x))]# ar1
  names(phi) <- naming_state("AR1", K)

  # alphs
  alpha <- x[grepl("^alpha", names(x))]
  if (length(alpha) == 0){ alpha <- NULL } else { names(alpha) <- names_list[["alpha"]]}

  # beta
  beta <- x[grepl("beta", names(x))]
  if (length(beta) == 0){beta <- NULL} else {names(beta) <- naming_state(names_list[["beta"]], K)}

  # out.lst para
  out.lst$init_prob <- init_prob
  out.lst$tp <- tp
  out.lst$intercept <- mu
  out.lst$AR1 <- phi
  out.lst$alpha <- if (is.null(alpha)) NULL else alpha
  out.lst$beta <- if (is.null(beta)) NULL else beta

  # return
  return(out.lst)
}


#' naming_fun
#'
#' @param x is a vector of names
#' @param para_names is a list of objects with names

naming_fun <- function(x, para_names){
  res.start <- grep("logalpha\\[.+", x)[1]
  res.end   <- length(x)
  res.para <- x[res.start:res.end]
  name_list <- c(names(para_names$init_prob), names(para_names$tp), names(para_names$intercept),
                 names(para_names$AR1), names(para_names$alpha), names(para_names$beta), "sigma" , res.para)
  return(name_list)
}






#' default_stan_control
#'
#' Sets controls to default unless otherwise specified.


default_stan_control <- function (adapt_delta = NULL, max_treedepth = 15L) {
    adapt_delta = NULL
    if (is.null(adapt_delta)){
        adapt_delta <- 0.95
    }
    nlist(adapt_delta, max_treedepth)
}


#' pars_include
#'
#' Set parameters to include in output.

pars_include <- function(Mx_a = 0, Mx_b = 0, Mx_c = 0, Mx_d = 0, Mx_e = 0){
  pars <- c()
  if (Mx_a == 0) pars <- c(pars, "gamma")
  if (Mx_b == 0) pars <- c(pars, "delta")
  if (Mx_c == 0) pars <- c(pars, "lambda")
  if (Mx_d == 0) pars <- c(pars, "alpha")
  if (Mx_e == 0) pars <- c(pars, "beta")
  return(pars)
}

#' set_sampling_args
#' From rstanarm (November 1st, 2019)
#'
#' Set the sampling arguments

set_sampling_args <- function (object, user_dots = list(), user_adapt_delta = NULL, ...) {
    user_adapt_delta <- NULL
    args <- list(object = object, ...)
    unms <- names(user_dots)
    for (j in seq_along(user_dots)) {
        args[[unms[j]]] <- user_dots[[j]]
    }
    defaults <- default_stan_control(adapt_delta = user_adapt_delta)
    if (!"control" %in% unms) {
        args$control <- defaults
    }
    else {
        if (!is.null(user_adapt_delta)) {
            args$control$adapt_delta <- user_adapt_delta
        }
        else {
            args$control$adapt_delta <- defaults$adapt_delta
        }
        if (is.null(args$control$max_treedepth)) {
            args$control$max_treedepth <- defaults$max_treedepth
        }
    }
    args$save_warmup <- FALSE
    return(args)
}


# Get the correct column name to use for selecting the median
# From rstanarm (November 1st, 2019)
#
# @param algorithm String naming the estimation algorithm (probably
#   \code{fit$algorithm}).
# @return Either \code{"50%"} or \code{"Median"} depending on \code{algorithm}.
select_median <- function(algorithm) {
  switch(algorithm,
         sampling = "50%",
         meanfield = "50%",
         fullrank = "50%",
         optimizing = "Median",
         stop("Bug found (incorrect algorithm name passed to select_median)",
              call. = FALSE))
}


# Issue warning if high rhat values
# From rstanarm (November 1st, 2019)
#
# @param rhats Vector of rhat values.
# @param threshold Threshold value. If any rhat values are above threshold a
#   warning is issued.
check_rhats <- function(rhats, threshold = 1.1, check_lp = FALSE) {
  if (!check_lp)
    rhats <- rhats[!names(rhats) %in% c("lp__", "log-posterior")]

  if (any(rhats > threshold, na.rm = TRUE))
    warning("Markov chains did not converge! Do not analyze results!",
            call. = FALSE, noBreaks. = TRUE)
}


#' check_stanfit
#'
#'
#' Check that a stanfit object (or list returned by rstan::optimizing) is valid

check_stanfit <- function(x) {
  if (is.list(x)) {
    if (!all(c("par", "value") %in% names(x)))
      stop("Invalid object produced please report bug")
  }
  else {
    stopifnot(is(x, "stanfit"))
    if (x@mode != 0)
      stop("Invalid stanfit object produced please report bug")
  }
  return(TRUE)
}

# Wrapper for rstan::summary
# @param stanfit A stanfit object created using rstan::sampling or rstan::vb
# @return A matrix of summary stats
make_stan_summary <- function(stanfit) {
  levs <- c(0.5, 0.8, 0.95)
  qq <- (1 - levs) / 2
  probs <- sort(c(0.5, qq, 1 - qq))
  rstan::summary(stanfit, probs = probs, digits = 10)$summary
}


# Test family
#
# @param x A character vector (probably x = family(fit)$family)
is.gaussian <- function(x) x == "gaussian"


# Check family argument
#
# @param f The \code{family} argument specified by user (or the default).
# @return If no error is thrown, then either \code{f} itself is returned (if
#   already a family) or the family object created from \code{f} is returned (if
#   \code{f} is a string or function).
validate_family <- function(f) {
  if (is.character(f))
    f <- get(f, mode = "function", envir = parent.frame(2))
  if (is.function(f))
    f <- f()
  if (!is(f, "family"))
    stop("'family' must be a family.", call. = FALSE)

  return(f)
}



#' deparse_call_formula
#'
#'
#' @param x A language object retrieved from call

deparse_call_formula <- function(x){
  out <- deparse(x)
  if (out == "NULL"){
    out <- NULL
    return(out)
  } else {
    return(out)
  }
}







# return.coef
#
# Return coefficients from Stan fit
#
# @param fit_extract: A list object from rstan::extract(fit)
# @param names_list: A list object referencing the variable names associated with the greek letters used in the Stan code

return.coef <- function(fit_extract = fit_extract, names_list = names_list){
  # initialize output mat
  coefficients <- matrix(NA, ncol = 3, nrow = 0)
  colnames(coefficients) <- c("state", "median", "MAD_sd")
  names_par <- c("gamma", "lambda" , "beta", "alpha", "mu", "phi")

  # extract coef of interest
  for (i in names_par){
    if (!is.null(fit_extract[[i]])){
      tmp <- fit_extract[[i]]
      tmp.names <- names_list[[i]]
      dim2 <- dim(tmp)[2]
      dim3 <- dim(tmp)[3]
      if (!is.na(dim3)){
        for (j in 1:dim2){
          for (q in 1:dim3){
            tmp.mat <- matrix(c(q, mean(tmp[ , j, q]), MAD_sd(tmp[ , j, q])), ncol = 3, nrow = 1)
            rownames(tmp.mat) <- paste0("S", q , "_", tmp.names[j], sep = "")
            coefficients <- rbind(coefficients, tmp.mat)
          }
        }
      } else {
        for (j in 1:dim2){
          tmp.mat <- matrix(c(0, mean(tmp[, j]), MAD_sd(tmp[, j])), ncol = 3, nrow = 1)
          rownames(tmp.mat) <- paste0("Sall_", tmp.names[j], sep = "")
          coefficients <- rbind(coefficients, tmp.mat)
        }
      }
    }
  }
  # return
  return(coefficients)
}


# MAD_sd
#
# @param x A numeric vector.
MAD_sd <- function(x){
  tilde_x <- median(x)
  x <- abs(x - tilde_x)
  out <- 1.482 * median(x)
  return(out)
}
