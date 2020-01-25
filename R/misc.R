#'
#'
#'
#' #' turn_type
#' #'
#'
#' form <- "y ~ 1 + x2 + as.factor(x1)"
#'
#' spl_form <- c("y", "1", "x2", "as.factor(x1)")
#'
#' x1 <- sample(c(1,2,3,4), 20, replace = TRUE)
#'
#' a <- "as.factor(x1)"
#'
#' gsub("\\)", "", gsub("as.factor\\(", "", a))
#'
#'
#' data <- as.data.frame(as.matrix(t(rmultinom(20, 1, rep(1/3, 3))) %*% c(1, 2, 3), ncol = 1, nrow = 20))
#' colnames(data) <- "x1"
#' for (i in a){
#'
#'   tmp <- a
#'
#' }
#'



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

########################################################
# Parsing the data
########################################################

#' data_parse
#'
#' @param formula_continuous
#' @param formula_discrete
#' @param state_varying_continuous
#' @param state_varying_discrte
#' @param data
#' @param n
#' @param t
#' @param j
#' @param q
#' @param K
#'
#' @examples
obs <- 200
N <- 20
T <- 10
J <- 4
Q <- 2
data <- data.frame(y = rnorm(obs), x1 = rnorm(obs), z1 = rnorm(obs), n = sort(rep(seq(1, N), T)), t = rep(seq(1, T), N))
jq <- data.frame(n = seq(1:N), j = sample(1:4, N, replace = TRUE))
jq <- merge(jq, data.frame(j = seq(1, J), q = sample(1:Q, J, replace = TRUE)), by = "j")
data <- merge(data, jq, by = "n")
formula_continuous <- "y ~ 1 + x1"
formula_discrete <- "z1"
state_varying_continuous <- "Intercept"
state_varying_discrete <- list(state = "z1", state_state = NULL)
K <- 3
data <- data_parse(formula_continuous = formula_continuous, formula_discrete = formula_discrete,
           state_varying_continuous = state_varying_continuous, state_varying_discrete = state_varying_discrete,
            data = data, n = "n", j = "j", q = "q", t = "t", K = K)

data_parse <- function(formula_continuous, formula_discrete,
                       state_varying_continuous, state_varying_discrete, data,
                       n = NULL, j = NULL, t = NULL, q = NULL, K = K){
  data_lst <- list(y = NULL, x_a = NULL, x_b = NULL, x_c = NULL, x_d = NULL, x_e = NULL, has_intercept = NULL, n = NULL, t = NULL, N = NULL)
  names_lst <- list(y = NULL, x_a = NULL, x_b = NULL, x_c = NULL, x_d = NULL, x_e = NULL)
  form_con <- if (!is.null(formula_continuous)) .split_formula(formula_continuous) else NULL
  form_dis <- if (!is.null(formula_discrete)) .split_formula(formula_discrete) else NULL

  # Intercept
  data_lst$has_intercept <- .has_intercept(form_con = form_con, form_dis = form_dis, state_varying_continuous = state_varying_continuous, state_varying_discrete = state_varying_discrete)

  # State process
  .j_define_out <- .j_define(data = data, j = j, n = n)
  data_lst$J <- .j_define_out$J
  data_lst$j <- .j_define_out$j
  data_lst$id_tp <- .state_probabilities(data = data, j = "j", q = "q")

  # NA
  na_df <- .sort_and_na_extend(data = data, n = "n", t = "t", j = "j")
  data <- na_df$data
  data_lst$id_miss <- na_df$id_miss

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
  data_lst$t <- data[, t]
  data_lst$n <- .n_adj(n = "n", data = data)
  data_lst$T <- max(unique(data_lst$t))
  data_lst$N <- length(unique(data_lst[["n"]]))
  data_lst$Q <- length(unique(data[, q]))
  data_lst$NT <- nrow(data)

  # K
  data_lst$K <- K

  return(list(data_lst = data_lst, names_lst = names_lst))
}

#' return_id_miss
#' @param df A dataframe which potentially contains missing values

.return_id_miss <- function(df){
  id_miss <- as.integer(apply(is.na(df), MARGIN = 1, function(x) TRUE %in% x))
  return(id_miss)
}

#' .na_extend
#' @param df A dataframe containing n and z
#' @param n The unit level indicator
#' @param t The timepoint indicator
#' @param j The group level indicator
#'
#' @examples
#' obs <- 200
#' N <- 5
#' J <- 3
#' data <- data.frame(x = rnorm(obs), n = sort(rep(seq(1, N), obs/N)), t = rep(seq(1, obs/N), N))
#' data <- merge(data, data.frame(n = seq(1, N), j = sample(seq(1, J), N, replace = T)))
#' data <- data[sample(seq(1,obs), round(0.8 * obs)) , ]
#' out <- .sort_and_na_extend(data = data, n = "n", j = "j", t = "t")

.sort_and_na_extend <- function(data, n = NULL, t = NULL, j = NULL){
  data <- data %>% rename(t = !!sym(t), n = !!sym(n), j = !!sym(j)) %>%
    arrange(j,n,t) %>%
    group_by(j) %>%
    complete(nesting(j), t = seq(min(t), max(t))) %>% ungroup()
  data <- as.data.frame(data)
  id_miss <- .return_id_miss(data)
  data[is.na(data)] <- 0
  return(list(data = data, id_miss = id_miss))
}


#' .j_define
#'
#' @param j Grouping variable
#' @param n Unit variable
#'
#' @details j Can be either individual specific (everyone has their own state process),
#' shared (everyone has the same state process) or based on indicators.
#' @examples
#' N <- 50
#' j <- "individual"
#' n <- sample(seq(1, 4), N, replace = TRUE)
#' .j_define(j = j, n = n)
#'
#' j <- "shared"
#' n <- sample(seq(1, 4), N, replace = TRUE)
#' .j_define(j = j, n = n)
#'
#' j <- sample(seq(1,2), N, replace = TRUE)
#' n <- sample(seq(1, 4), N, replace = TRUE)
#' .j_define(j = j, n = n)

.j_define <- function(data = data, j = j, n = n){
  j <- data[, j]
  n <- data[, n]
  if (j[1] == "individual"){
    out <- n # give everyone their own state process
  } else if (j[1] == "shared"){
    out <- rep(1, length(n)) # give everyone the same state process
  } else {
    out <- match(j, unique(j)) #
  }
  J <- length(unique(out)) # number of state processes
  return(list(j = out, J = J))
}

#' .state_probabilities
#'
#' @param j The state process indicator.
#' @param q The transition probability indicator.
#'
#' @details This function creates a vector of length J assigning the state
#' processes to specific transition probabilities.
#'
#' @examples
#' t <- data.frame(j = seq(1, 4), q = rep(seq(1, 2), 2))
#' N <- 20
#' data <- t[sample(seq(1, 4), N, replace = TRUE), ]
#' .state_probabilities(j = data[, 1], q = data[, 2])
#' .state_probabilities(j = data[, 1], q = NULL)
#' .state_probabilities(j = data[, 1], q = "individual")
#' .state_probabilities(j = data[, 1], q = "shared")

.state_probabilities <- function(data = data, j = j, q = q){
  q <- data[, q]
  j <- data[, j]
  if (is.null(q)){
    out <- unique(sort(j))
  } else {
    if (q[1] == "individual"){
      out <- seq(1, length(unique(j)))
    } else if (q[1] == "shared"){
      out <- rep(1, length(unique(j)))
    } else {
      data <- cbind(j, q)
      data <- data[order(data[, 1]),]
      out <- unique(data) # unique combinations
      if (length(out[, 1]) != length(unique(out[, 1]))){
        stop("Assigned transition probabilities are not unique.")
      }
      out <- out[, 2]
    }
  }
  return(out)
}

#' .add_intercept_name
#'
#' @param has_intercept Flag indicating whether the parameter set has an intercept.
#' @param names Names of the parameters
#'
#' @return Return the names of the parameters with `Intercept` appended at the front.

.add_intercept_name <- function(has_intercept, names){
  out <- names
  if (has_intercept == 1) out <- c("Intercept", out)
  return(out)
}

#' .n_adj
#'
#' @param n A vector of length N indicating group association.
#' @param data The data vector
#'
#' @return Returns a length n vector of 1, ..., J flags.

.n_adj <- function(n = NULL, data){
  if (is.null(n)){
    n <- rep(1, nrow(data))
  } else {
    n <- data[, n]
    unique.n <- unique(n)
    n <- match(n, unique.n)
    names(n) <- as.character(n)
  }
  return(n)
}

#' .factor_mat
#'
#' @param x
#' @param name
#'
#' @examples
#' N <- 20
#' x <- as.factor(sample(seq(1, 4), N, replace = TRUE))
#' .factor_mat(x, "NAME")

.factor_mat <- function(x, name){
  tmp.lvl <- levels(x)
  tmp.mat <- matrix(0, nrow = length(x), ncol = length(tmp.lvl))
  for (j in 1:length(tmp.lvl)) tmp.mat[, j][x == tmp.lvl[j]] <- 1
  colnames(tmp.mat) <- paste0(name, "_", tmp.lvl)
  if (all(rowSums(tmp.mat) == 1)) tmp.mat <- tmp.mat[, 2:length(tmp.lvl)]
  x <- tmp.mat
  if (is.null(dim(x))){
    x <- as.matrix(x)
    colnames(x) <- name
  }
  return(x)
}

#' .interaction
#'
#' @param x The variable name
#' @param data The data
#'
#' @examples
#' N <- 20
#' x <- "x1:x2"
#' data <- data.frame(x1 = as.factor(sample(seq(1, 4), N, replace = TRUE)), x2 = rnorm(N))
#' .interaction(x = x, data = data)

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

#' .create_predmat
#'
#' @param names Predictor names
#' @param data The dataframe
#'
#' @details Based on the parsed formula adjusts the data frame

.create_predmat <- function(names, data){
  out <- c()
  names <- names[!names == "Intercept"]
  for (i in names){
    tmp <- i
    if (grepl(":", tmp)){
      tmp.var <- .interaction(x = tmp, data = data)
    } else if (is.factor(data[, tmp])){ # flag for factor
      tmp.var <- data[, tmp]            # tmp for factor variable
      tmp.var <- .factor_mat(tmp.var, tmp) # split factor
    } else {
      tmp.var <- as.matrix(data[, tmp])
      colnames(tmp.var) <- i
    }
    out <- cbind(out, tmp.var) # append splitted factor
  }
  return(out)
}

#' .has_intercept

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
    } else if (is.element("Intercept", state_varying_discrete[["state_state"]])){
      has_intercept[5] <- 1
    } else {
      has_intercept[3] <- 1
    }
  }
  return(has_intercept)
}

#' .split_formula
#'
#' @param x The formula
#'
#' @details 1.
#'
#' @examples
#' x <- "y ~ 1 + x1#x2"
#' .split_formula(x)
#'
#' x <- "y ~ x1:x2"
#' .split_formula(x)

.split_formula <- function(x){
  if (!is.null(x)){
    x <- gsub("\\s", "", x) # Removes whitespace
    x <- unique(strsplit(x, split = "~|\\+" )[[1]]) # remove non-unique elements
    x[x == "1"] <- "Intercept" # add intercept
    out <- c() # container for variables
    for (i in x){
      tmp <- i # tmp
      if (grepl("#", i)) tmp <- c(strsplit(i, "#")[[1]], gsub("#", ":", i)) # split on interaction sign and sub interaction sign with singular sign
      out <- c(out, tmp) # append
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


########################################################
# Slicers
########################################################

#' slicer_time
#'
#' @param j The group indicator
#'
#' @details Returns a matrix that contains the start and stop for each state process
#' and further grabs elements from the matrix that contains the start and stop of individual
#' time points
#'
#' @examples
#' j <- sort(rep(seq(1, 10), 20))
#' slicer_time(j = j)

slicer_time <- function(j){
  J <- length(unique(j))
  slicer_T <- matrix(NA, ncol = 2, nrow = J)
  for (j in 1:J){
    range <- which(j == j)
    slicer_T[j,] <- c(min(range), max(range))
  }
  return(slicer_T)
}

#' start_stop_slicer
#'
#' @examples
#' N <- 360
#' j <- sample(seq(1, 3), N, replace = TRUE)
#' t <- sort(sample(rep(seq(1, 60), 6)))
#' start_stop_slicer(j = j, t = t)

start_stop_slicer <- function(j = j, t = t){
  J <- length(unique(j))
  T <- max(t) - min(t)
  start_stop <- matrix(NA, ncol = 2 * J, nrow = T)
  for (j in 1:J){
    for (t in 1:T){
      range <- which(t == t & j == j)
      if (length(range) != 0) start_stop[t, ((j * 2) - 1):(j * 2)] <- c(min(range), max(range))
    }
  }
  return(start_stop)
}



########################################################
# Stan data
########################################################

#' prepare_standata
#'
#'
#' order_continuous <- c("Intercept")

prepare_standata(data = data$data_lst, priors = list(priors = c(""), A_prior = c("")),
                 order_continuous = order_continuous, state_sigma = FALSE, tvtp = TRUE)

prepare_standata <- function(data = data, priors = priors, order_continuous = order_continuous,
                             state_sigma = state_sigma, tvtp = tvtp){

    # NT
    N <- data$N
    n <- data$n
    T <- data$T
    NT <- data$NT

    # model specifications
    K <- data$K
    has_intercept <- data$has_intercept
    id_miss <- data$id_miss

    # priors
    A_prior <- priors[["A_prior"]]
    priors <- priors[["priors"]]

    # Coerce to matrixes
    x_d <- data$x_d
    x_e <- data$x_e
    z <- cbind(data$x_a, data$x_b, data$x_c)

    # dimensions
    Mx_d <- if (!is.null(x_d)) ncol(x_d) else 0
    Mx_e <- if (!is.null(x_e)) ncol(x_e) else 0
    Mz   <- if (!is.null(z)) ncol(z) else 0

    # fix if NULL
    if (Mx_d == 0) x_d <- matrix(0, ncol = 0, nrow = NT)
    if (Mx_e == 0) x_e <- matrix(0, ncol = 0, nrow = NT)
    if (Mz == 0) z <- matrix(0, ncol = 0, nrow = NT)

    # time-varying transition probabilities
    pp1 <- if (!is.null(data$x_a)) ncol(data$x_a) else 0
    pp2 <- if (!is.null(data$x_b)) ncol(data$x_b) else 0
    pp3 <- if (!is.null(data$x_c)) ncol(data$x_c) else 0
    pp_lambda <- matrix(0, nrow = 3, ncol = Mz)
    if (pp1 != 0) pp_lambda[1, 1:pp1] <- 1
    if (pp2 != 0) pp_lambda[2, (pp1 + 1):(pp1 + pp2)] <- 1
    if (pp3 != 0) pp_lambda[3, (pp1 + pp2 + 1):(pp1 + pp2 + pp3)] <- 1
    pp_gamma <- pp_lambda[2, (pp1 + 1):ncol(pp_lambda)]
    pp_eta <- pp_lambda[3, (pp1 + 1):ncol(pp_lambda)]

    # order vector
    order_x_e <- .create_order_vector(data, order_continuous)

    # output
    y <- data$y

    # state process
    NS <- data$J
    NTP <- data$Q
    id_tp <- data$id_tp

    # slicer
    slicer_T <- slicer_time(j = data$j)
    start_stop <- start_stop_slicer(j = data$j, t = data$t)

    # standata
    standata <- list(
      N = N, # units
      T = T, # max of time points
      NT = NT, # N of observations
      NS = NS, # N of state processes
      NTP = NTP, # N of transition probabilities
      id_tp = id_tp,      # [NS] which state process belongs to which tp matrix
      slicer_T = slicer_T,   # min(N, T):max(N, T)
      startstop = start_stop,  # slicer for s units
      K = K,         # N of states
      has_intercept = has_intercept, # [5] 1: alpha, 2: beta; 3: gamma; 4: delta; 5: eta
      Mz = Mz, # N of tp predictors
      Mx_d = Mx_d, # N of fixed parameters of continuous process
      Mx_e = Mx_e, # N of continuous parameters of continuous process
      pp1 = pp1, pp2 = pp2, pp3 = pp3, # N of general, state, and state-state specific predictors
      pp_lambda = pp_lambda, # [3, Mz] which are varying at which level for tps
      pp_gamma = pp_gamma, # [pp2] 0/1 of varying at state
      pp_eta = pp_eta, # [pp2]   0/1 of varying at state-state
      z = z,      # [NTP * T, Mz] matrix of predictors of discrete process
      x_d = x_d,    # [NT, Mx_d] matrix of fixed predictors of continuous process
      x_e = x_e,    # [NT, Mx_e] matrix of varying predictors of continuous process
      y = y,      # [NT] vector of output
      state_sigma = as.integer(state_sigma), # 0: general,    1: state-specific
      tvtp = as.integer(tvtp),
      order_x_e = order_x_e, # [Mx_e + has_intercept[2]] 0: unordered, 1: ordered
      A_prior = A_prior, # [K] A_prior;
      priors = priors,  #[7,4];     // 1: Kind, 2: mean, 3: sd, 4: df; 1: normal, 2: cauchy, 3: student-t
      id_miss = id_miss  # 1: missing, 0: present / at least one observation
    )

  return(standata)
}

#' .create_order_vector
#'
#' @param formula The formula list
#' @param order_continuous A character vector indicating which parameters are ordered.
#'
#' This function creates a 0/1 vector indicating whether a particular parameter that is varying across states is ordered.

.create_order_vector <- function(data, order_continuous){
  intercept <- if (data$has_intercept[2] == 1) "Intercept" else NULL
  tmp <- c(intercept, colnames(data$x_e))
  out <- rep(0, length(tmp))
  out[match(order_continuous, tmp)] <- 1
  return(out)
}

#' initialization

# if (init.prior) {
#   initf2 <- function(chain_id = 1) {
#     list(
#       pi1 = matrix(rep(1/K, N * K), ncol = K, nrow = N),
#       A = array(rdirichlet(K * N_, rep(1, K)), dim = c(N_, K, K)),
#       phi = matrix(sort(rnorm(K, 0.5, 0.25), decreasing = FALSE), ncol = K, nrow = 1),
#       alpha = array(rnorm(Mx_d_, 0, 1)),
#       beta = matrix(rnorm(K * Mx_e_, 0, 1), ncol = K, nrow = Mx_e_),
#       sigma = 1)
#     }
#   init_ll <- lapply(1:n_chains, function(id) initf2(chain_id = id))
#   init = init_ll
# } else {
#   init = "random"
# }












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

########################################################
# Check data
########################################################

#' check_order
#'
#' @param data Included data frame
#' @param order_continuous Vector of parameter names declared to be continuous

check_data <- function(data, order_continuous, formula_continuous, formula_discrete, n, t,
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
  .check_nt(names_dta = names_dta, n = n, t = t)
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

.check_nt <- function(names_dta, n, t){
  if (!(is.element(n, names_dta) & is.element(t, names_dta))){
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

split_naming <- function(x, names_list, N, K, shared_TP = TRUE, state_SD = FALSE){

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

  # sigma
  sigma <- x[grepl("sigma", names(x))]
  if (state_SD) names(sigma) <- naming_state("sigma", K) else names(sigma) <- "sigma"

  # out.lst para
  out.lst$init_prob <- init_prob
  out.lst$tp <- tp
  out.lst$AR1 <- phi
  out.lst$alpha <- if (is.null(alpha)) NULL else alpha
  out.lst$beta <- if (is.null(beta)) NULL else beta
  out.lst$sigma <- sigma

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
                 names(para_names$AR1), names(para_names$alpha), names(para_names$sigma), names(para_names$beta), res.para)
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

pars_include <- function(pp1 = 0, pp2 = 0, pp3 = 0, Mx_d = 0, Mx_e = 0){
  pars <- c("beta_un", "beta_ord")
  if (pp1 == 0) pars <- c(pars, "gamma")
  if (pp2 == 0) pars <- c(pars, "delta")
  if (pp3 == 0) pars <- c(pars, "lambda")
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
