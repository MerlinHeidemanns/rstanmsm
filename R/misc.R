




#' check_tp_s
#'
#' @param shared_TP TRUE/FALSE indicating whether transition probabilities are being shared.
#' @param shared_S TRUE/FALSE indicating whether the discrete process is being shared
#'
#' This function checks for the case that the user indicates that transition probabilities are individual but the state process is shared
#' which is not possible. Returns an error.

check_tp_s <- function(shared_TP = NULL, shared_S = NULL, n = NULL){
  if (isTRUE(shared_S) & isFALSE(shared_TP) & isTRUE(length(unique(n)) > 1)){
    stop("It is not possible to specify the discrete process as shared and the transition probabilities to be unit-specific for N = 1.")
  }
}

#' data_parse
#'
#'

data_parse <- function(formula_continuous, formula_discrete, state_varying_continuous, state_varying_discrete, data, n_var, t_var, K){
  data_lst <- list(y = NULL, x_a = NULL, x_b = NULL, x_c = NULL, x_d = NULL, x_e = NULL, has_intercept = NULL, n_var = NULL, t_var = NULL, N = NULL)
  names_lst <- list(y = NULL, x_a = NULL, x_b = NULL, x_c = NULL, x_d = NULL, x_e = NULL)
  form_con <- if (!is.null(formula_continuous)) .split_formula(formula_continuous) else NULL
  form_dis <- if (!is.null(formula_discrete)) .split_formula(formula_discrete) else NULL

  # Intercept
  data_lst$has_intercept <- .has_intercept(form_con = form_con, form_dis = form_dis, state_varying_continuous = state_varying_continuous, state_varying_discrete = state_varying_discrete)

  # y
  data_lst$y <- data[, form_con[1]]
  names_lst$y <- form_con[1]
  form_con <- form_con[2:length(form_con)]

  # create names
  names_x_d <- form_con[!grepl(paste0(paste0("^" ,state_varying_continuous, "$"), collapse = "|"), form_con)]
  names_x_e <- state_varying_continuous
  names_x_a <- form_dis[!grepl(paste0(paste0("^" ,unlist(state_varying_discrete), "$") , collapse = "|"), form_dis)]
  names_x_b <- state_varying_discrete[["state"]]
  names_x_c <- state_varying_discrete[["state_state"]]

  # create predictor matrixes
  data_lst$x_d <- .create_predmat(names_x_d, data)
  data_lst$x_e <- .create_predmat(names_x_e, data)
  data_lst$x_a <- .create_predmat(names_x_a, data)
  data_lst$x_b <- .create_predmat(names_x_b, data)
  data_lst$x_c <- .create_predmat(names_x_c, data)

  # names
  names_lst$x_d <- .add_intercept_name(data_lst$has_intercept[1], colnames(data_lst$x_d))
  names_lst$x_e <- .add_intercept_name(data_lst$has_intercept[2], colnames(data_lst$x_e))
  names_lst$x_a <- .add_intercept_name(data_lst$has_intercept[3], colnames(data_lst$x_a))
  names_lst$x_b <- .add_intercept_name(data_lst$has_intercept[4], colnames(data_lst$x_b))
  names_lst$x_c <- .add_intercept_name(data_lst$has_intercept[5], colnames(data_lst$x_c))

  # n and t
  data_lst$t_var <- data[, t_var]
  data_lst[["n_var"]] <- .n_var_adj(n_var = n_var, data)
  data_lst$T <- max(unique(data_lst$t_var))
  data_lst$N <- length(unique(data_lst[["n_var"]]))

  # K
  data_lst$K <- K

  return(list(data_lst = data_lst, names_lst = names_lst))
}

.add_intercept_name <- function(has_intercept, names){
  out <- names
  if (has_intercept == 1) out <- c("Intercept", out)
  return(out)
}

.n_var_adj <- function(n_var = NULL, data){
  if (is.null(n_var)){
    n_var <- rep(1, nrow(data))
  } else {
    n <- data[, n_var]
    unique.n <- unique(n)
    n_var <- match(n, unique.n)
    names(n_var) <- as.character(n)
  }
  return(n_var)
}

.factor_mat <- function(x, name){
  tmp.lvl <- levels(x)
  tmp.mat <- matrix(0, nrow = length(x), ncol = length(tmp.lvl))
  for (j in 1:length(tmp.lvl)){
    tmp.mat[, j][x == tmp.lvl[j]] <- 1
  }
  colnames(tmp.mat) <- paste0(name, "_", tmp.lvl)
  if (all(rowSums(tmp.mat) == 1)){
    tmp.mat <- tmp.mat[, 2:length(tmp.lvl)]
  }
  x <- tmp.mat
  if (is.null(dim(x))){
    x <- as.matrix(x)
    colnames(x) <- name
  }
  return(x)
}

.interaction <- function(x, data){
  var <- strsplit(x, ":")[[1]] # split interaction into constituent terms
  var.lst <- list() # container for interaction terms
  for (i in 1:2){
    var.name <- var[i] # name of the variable
    var.tmp <- data[, var.name] # variable
    if (is.factor(var.tmp)){
      var.lst[[i]] <- .factor_mat(var.tmp, var.name) # create factor
    } else {
      var.lst[[i]] <- as.matrix(var.tmp, ncol = 1) # if not factor then make matrix
      colnames(var.lst[[i]]) <- var.name # name matrix
    }
  }
  out <- c() # initialize output matrix
  lng.lst <- length(var.lst) # number of interaction terms
  pred <- rep(NA, lng.lst) # initialize empty vector to hold N(cols)
  for (q in 1:lng.lst){
    pred[q] <- if (is.null(ncol(var.lst[[q]]))) 1 else ncol(var.lst[[q]]) # assign N of columns for cycling
  }
  for (i in 1:pred[1]){   # cycle through first matrix
    for (j in 1:pred[2]){ # cycle through second matrix
      tmp <- as.matrix(var.lst[[1]][, i] * var.lst[[2]][, j]) # create interaction term
      colnames(tmp) <- paste0(colnames(var.lst[[1]])[i], ":", colnames(var.lst[[2]])[j]) # rename
      out <- cbind(out, tmp) # add to output matrix
    }
  }
  return(out)
}

.create_predmat <- function(names, data){
  out <- c()
  names <- names[!names == "Intercept"]
  for (i in names){
    tmp <- i
    if (grepl(":", tmp)){
      tmp.var <- .interaction(x = tmp, data = data)
    } else if (is.factor(data[, tmp])){
      tmp.var <- data[, tmp]
      tmp.var <- .factor_mat(tmp.var, tmp)
    } else {
      tmp.var <- as.matrix(data[, tmp])
      colnames(tmp.var) <- i
    }
    out <- cbind(out, tmp.var)
  }
  return(out)
}

.has_intercept <- function(form_con, form_dis, state_varying_continuous, state_varying_discrete){
  has_intercept <- rep(0, 5)
  if (is.element("Intercept", form_con)){
    if (is.element("Intercept", state_varying_continuous)){
      has_intercept[2] <- 1
    } else {
      has_intercept[1] <- 1
    }
  }
  if (is.element("Intercept", form_dis)){
    if (is.element("Intercept", state_varying_discrete[["state"]])){
      has_intercept[4] <- 1
    } else if (is.element("Intercept", state_varying_discrete[["state_staet"]])){
      has_intercept[5] <- 1
    } else {
      has_intercept[3] <- 1
    }
  }
  return(has_intercept)
}


.split_formula <- function(x){
  if (!is.null(x)){
    x <- gsub("\\s", "", x)
    x <- unique(strsplit(x, split = "~|\\+" )[[1]])
    x[x == "1"] <- "Intercept"
    out <- c()
    for (i in x){
      tmp <- i
      if (grepl("#", i)){
        tmp <- c(strsplit(i, "#")[[1]], gsub("#", ":", i))
      }
      out <- c(out, tmp)
    }
  } else {
    out <- NULL
  }
  return(out)
}

# data <- as.data.frame(matrix(t(rmultinom(10, 1, rep(1/3, 3))) %*% c(1, 2, 3))); colnames(data) <- "x1"
# x <- "as.factor(x1)"
# .adj_type <- function(x, data){
#   # factor
#   for (i in x){
#     factor <- if ()
#     if (grepl("as\\.factor(.+)", i)){
#       tmp <- gsub("as\\.factor((.+))", "", i)
#       data[, tmp] <- as.factor(data[, tmp])
#     }
#   }
#   return(data)
# }
# .adj_type(x = x, data = data)

#' formula_parse
#'
#' Parses the formulas into their components and stores them in a list




.sort_predictor <- function(x){
  out.lst <- list(base = NULL, state = NULL, state.state = NULL, all.var = NULL)
  if (grepl("\\*", x)){
    tmp <- strsplit(x, "\\*")[[1]]; tmp <- sort(tmp); tmp <- c(tmp, paste(tmp[1], tmp[2], sep = ":"))
  } else if (grepl(":", x)){
    tmp <- strsplit(x, ":")[[1]]; tmp <- sort(tmp); tmp <- paste(tmp[1], tmp[2], sep = ":")
  } else {
    tmp <- x
  }
  for (j in tmp){
    if (grepl("\\|3", j)){
      out.lst$state.state <- c(out.lst$state.state, gsub("\\|3|\\|2", "", j))
    } else if (grepl("\\|2", j)){
      out.lst$state <- c(out.lst$state, gsub("\\|2", "", j))
    } else {
      out.lst$base <- c(out.lst$base, j)
    }
    if (!grepl("\\*", j)){
      if (grepl(":", j)){
        out.lst$all.var <- c(out.lst$all.var, gsub("\\|2|\\|3", "", strsplit(j, ":")[[1]]))
      } else {
      out.lst$all.var <- c(out.lst$all.var, gsub("\\|2|\\|3", "", j))
      }
    }
  }
  out.lst$all.var <- unique(out.lst$all.var)
  return(out.lst)
}

formula_parse <- function(formula_discrete = NULL, formula_continuous = formula_continuous) {

  int_pattern <- "(^1$)|(^1\\|2$)|(^1\\|3$)"

  # out list
  out.lst <- list(y = NULL, a = NULL, b = NULL, c = NULL, d = NULL, e = NULL,
                   has_intercept = rep(0, 5), all.var = NULL)

  # containers
  tmp.d.a <- tmp.d.b <- tmp.d.c <- tmp.c.d <- tmp.c.e <- tmp.all.var <- c()
  has_intercept <- rep(0, 5)

  # formula continuous
  tmp.c <- .split_formula(formula_continuous)
  for (i in tmp.c[2:length(tmp.c)]){
    if (grepl(int_pattern, i)){
      if (grepl("\\|2", i)) {has_intercept[2] <- 1} else {has_intercept[1] <- 1}
    } else {
      out <- .sort_predictor(i)
      tmp.c.d <- c(tmp.c.d, out$base)
      tmp.c.e <- c(tmp.c.e, out$state)
      tmp.all.var <- c(tmp.all.var, out$all.var)
    }
  }

  # formula discrete
  if (!is.null(formula_discrete)){
    tmp.d <- .split_formula(formula_discrete)
    for (i in tmp.d){
      if (grepl(int_pattern, i)){
        if (grepl("\\|3", i)){
          has_intercept[5] <- 1
        } else if (grepl("\\|2", i)){
          has_intercept[4] <- 1
        } else {
          has_intercept[3] <- 1
        }
      } else {
        out <- .sort_predictor(i)
        tmp.d.a <- c(tmp.d.a, out$base)
        tmp.d.b <- c(tmp.d.b, out$state)
        tmp.d.c <- c(tmp.d.c, out$state.state)
        tmp.all.var <- c(tmp.all.var, out$all.var)
      }
    }
  }
  out.lst$a <- unique(tmp.d.a)
  out.lst$b <- unique(tmp.d.b)
  out.lst$c <- unique(tmp.d.c)
  out.lst$d <- unique(tmp.c.d)
  out.lst$e <- unique(tmp.c.e)
  out.lst$y <- tmp.c[1]
  out.lst$has_intercept <- has_intercept
  out.lst$all.var <- unique(c(tmp.all.var, tmp.c[1]))
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

.create_interactions <- function(data, formula){
  seq.types <- c("a", "b", "c", "d", "e")
  for (i in seq.types){
    int <- formula[[i]][grepl(":", formula[[i]])]
    for (j in int){
      tmp <- strsplit(j, split = ":")[[1]]
      tmp <- data[, tmp[1]] * data[, tmp[2]]
      data <- cbind(data, tmp)
      colnames(data)[ncol(data)] <- j
    }
  }
  return(data)
}



data_split <- function(data = data, tvtp = FALSE, parsed_formula = parsed_formula, n = NULL, t = NULL, K = NULL) {

  # initialize output list
  out.lst <- list(a = NULL, b = NULL, c = NULL, d = NULL, e = NULL, y = NULL, n = NULL, t = NULL, K = K, N = NULL)
  if (!is.data.frame(data)){
    data <- as.data.frame(data)
  }
  names_data <- colnames(data)
  data <- data[order(data$n, data$t),]

  # create interactions
  data <- .create_interactions(data = data, formula = parsed_formula)

  if (tvtp) {
    out.lst$a <- as.matrix(data[, parsed_formula$a])
    colnames(out.lst$a) <- parsed_formula$a
    out.lst$b <- as.matrix(data[, parsed_formula$b])
    colnames(out.lst$b) <- parsed_formula$b
    out.lst$c <- as.matrix(data[, parsed_formula$c])
    colnames(out.lst$c) <- parsed_formula$c
  }

  out.lst$d <- as.matrix(data[, parsed_formula$d])
  colnames(out.lst$d) <- parsed_formula$d
  out.lst$e <- as.matrix(data[, parsed_formula$e])
  colnames(out.lst$e) <- parsed_formula$e

    out.lst$y <- data[ , parsed_formula$y]
  out.lst$n <- data[ , n]
  out.lst$t <- data[ , t]
  out.lst$N <- length(unique(out.lst$n))

  return(out.lst)
}




#' check_order
#'
#' @param data Included data frame
#' @param order_continuous Vector of parameter names declared to be continuous

check_data <- function(data, order_continuous, formula_continuous, formula_discrete, n_var, t_var,
                       state_varying_continuous = c(),
                       state_varying_discrete   = list()){
  var.con <- .split_formula(formula_continuous)
  var.dis <- .split_formula(formula_discrete)
  var <- c(var.con, var.dis)
  names_dta <- colnames(data)
  # check varying
  .check_varying(var.con, state_varying_continuous)
  .check_varying(var.dis, state_varying_discrete)
  # check inclusion
  .check_inclusion(names_dta, parsed_names = var)
  # check for n and t
  .check_nt(names_dta = names_dta, n_var = n_var, t_var = t_var)
  # check for order predictor
  .check_order(state_varying_continuous = state_varying_continuous, order_continuous = order_continuous)
}

.check_varying <- function(x, varying){
  if (is.list(x)){
    x <- unlist(a)
  }
  x[x == "1"] <- "Intercept"
  for (i in varying){
    if (!is.element(i, x)){
      stop(paste0("The varying predictor ", i, " is not among the predictors of the continuous model."))
    }
  }
}


.check_order <- function(state_varying_continuous, order_continuous){
  for (i in order_continuous){
    if (!is.element(i, state_varying_continuous)){
      stop(paste0("The predictor ", i, " which is supposed to be ordered across states is not indicated as varying by state."))
    }
  }
}

.check_nt <- function(names_dta, n_var, t_var){
  if (!(is.element(n_var, names_dta) & is.element(t_var, names_dta))){
    stop("Please supply indexes for n and t.")
  }
}

.check_inclusion <- function(names_dta, parsed_names){
  for (i in parsed_names){
    if (i != "Intercept"){
      if (!is.element(i, names_dta)) stop(paste0("The predictor ", i, " has not been found in the data."))
    }
  }
}


#' create_order_vector
#'
#' @param formula The formula list
#' @param order_continuous A character vector indicating which parameters are ordered.
#'
#' This function creates a 0/1 vector indicating whether a particular parameter that is varying across states is ordered.

create_order_vector <- function(data, order_continuous){
  intercept <- if (data$has_intercept[2] == 1) "Intercept" else NULL
  tmp <- c(intercept, colnames(data$x_e))
  out <- rep(0, length(tmp))
  out[match(order_continuous, tmp)] <- 1
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

  out.lst <- list(init_prob = NULL, tp = NULL, AR1 = NULL, alpha = NULL,
                  beta = NULL, gamma = NULL, delta = NULL, lambda = NULL)
  init_prob <- x[grepl("pi1\\[.+", names(x))]
  tp <- x[grepl("A\\[[0-9]+,[0-9]+,[0-9]+\\]", names(x))]

  # AR1
  phi <- x[grepl("phi\\[.+\\]", names(x))]# ar1
  names(phi) <- naming_state("AR1", K)

  # alpha
  alpha <- x[grepl("^alpha", names(x))]
  if (length(alpha) == 0){ alpha <- NULL } else { names(alpha) <- names_list[["x_d"]]}

  # beta
  beta <- x[grepl("beta", names(x))]
  if (length(beta) == 0){beta <- NULL} else {names(beta) <- naming_state(names_list[["x_e"]], K)}

  # out.lst para
  out.lst$init_prob <- init_prob
  out.lst$tp <- tp
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
  name_list <- c(names(para_names$init_prob), names(para_names$tp),
                 names(para_names$AR1), names(para_names$alpha), names(para_names$beta), "sigma" , res.para)
  return(name_list)
}


#' names_list_

naming_list <- function(x){
  out.lst <- list(alpha = NULL, beta = NULL, gamma = NULL, delta = NULL, lambda = NULL)
  if (x$has_intercept[1]) out.lst$alpha <- c("Intercept")
  if (x$has_intercept[2]) out.lst$beta <- c("Intercept")
  if (x$has_intercept[3]) out.lst$gamma <- c("Intercept")
  if (x$has_intercept[4]) out.lst$delta <- c("Intercept")
  if (x$has_intercept[5]) out.lst$lambda <- c("Intercept")

  if (length(x[["d"]]) != 0) out.lst$alpha <- c(out.lst$alpha, x[["d"]])
  if (length(x[["e"]]) != 0) out.lst$beta <- c(out.lst$beta, x[["e"]])
  if (length(x[["a"]]) != 0) out.lst$gamma <- c(out.lst$gamma, x[["a"]])
  if (length(x[["b"]]) != 0) out.lst$delta <- c(out.lst$delta, x[["b"]])
  if (length(x[["c"]]) != 0) out.lst$lambda <- c(out.lst$lambda, x[["c"]])
  return(out.lst)
}

#' default_stan_control
#'
#' Sets controls to default unless otherwise specified.


default_stan_control <- function (adapt_delta = NULL, max_treedepth = 15L) {
    adapt_delta = NULL
    if (is.null(adapt_delta)){
        adapt_delta <- 0.95
    }
    list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
}


#' pars_include
#'
#' Set parameters to include in output.

pars_include <- function(Mx_a = 0, Mx_b = 0, Mx_c = 0, Mx_d = 0, Mx_e = 0){
  pars <- c("beta_un", "beta_ord")
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






# MAD_sd
#
# @param x A numeric vector.
MAD_sd <- function(x){
  tilde_x <- median(x)
  x <- abs(x - tilde_x)
  out <- 1.482 * median(x)
  return(out)
}
