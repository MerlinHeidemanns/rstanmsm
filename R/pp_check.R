# pp_check

pp_check.stan_msm <- function(object, ..., type = ("mts")) {
 y <- object$data$y
 yrep <- object$fitted.values
 n <- object$data$n
 t <- object$data$t
  switch(match.arg(type),
         mts = ppc_mts(y, yrep, n, t)
         )
}

ppc_mts <- function(y, yrep, n, t, color = c("blue", "green")){
 #check_ignored_arguments
 group <- paste(fit$data$n, fit$data$t)
 y <- bayesplot:::validate_y(y)
 yrep <- bayesplot:::validate_yrep(yrep, y)
 group <- bayesplot:::validate_group(group, y)
 plot_data <- bayesplot:::ppc_group_data(y, yrep, group)
 plot_data <- plot_data %>% tidyr::separate(group, c("n", "t"), sep = "\\s") %>% dplyr::mutate(n = as.numeric(n), t = as.numeric(t))
  ggplot() +
    geom_line(data = plot_data %>% dplyr::filter(variable == "y"),
            aes(x = t, y = value, color = as.factor(n)), size = 1) +
    geom_line(data = plot_data %>% dplyr::filter(grepl("yrep", variable)),
            aes(x = t, y = value, group = interaction(as.factor(n),as.factor(variable)),color = as.factor(n)),
           alpha = 0.1, size = 0.1) +
    scale_colour_manual(values = color)
    guides(color = guide_legend(title = NULL)) +
    bayesplot_theme_get()
}

ppc_acf <- function(y, yrep, n, t){

}



## ----- Internal

#' validate_y
#' from bayesplot November 7th, 2019

validate_y <- function (y) {
    stopifnot(is.numeric(y))
    if (!(inherits(y, "ts") && is.null(dim(y)))) {
        if (!is_vector_or_1Darray(y)) {
            abort("'y' must be a vector or 1D array.")
        }
        y <- as.vector(y)
    }
    if (anyNA(y)) {
        abort("NAs not allowed in 'y'.")
    }
    unname(y)
}

validate_yrep <- function (yrep, y) {
    stopifnot(is.matrix(yrep), is.numeric(yrep))
    if (is.integer(yrep)) {
        if (nrow(yrep) == 1) {
            yrep[1, ] <- as.numeric(yrep[1, , drop = FALSE])
        }
        else {
            yrep <- apply(yrep, 2, as.numeric)
        }
    }
    if (anyNA(yrep)) {
        abort("NAs not allowed in 'yrep'.")
    }
    if (ncol(yrep) != length(y)) {
        abort("ncol(yrep) must be equal to length(y).")
    }
    unclass(unname(yrep))
}
