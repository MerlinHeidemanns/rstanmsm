# pp_check

#
# library(bayesplot)
#
# pp_check.stan_msm <- function(object, ..., type = c("multiple", "overlaid")) {
#   y <- object[["y"]]
#   yrep <- object[["fitted.values"]]
#   switch(match.arg(type),
#          multiple = ppc_hist(y, yrep[1:min(8, nrow(yrep)),, drop = FALSE]),
#          overlaid = ppc_dens_overlay(y, yrep))
# }
#
# ppc_mts <- function(y, yrep, group, stat = "mean", ..., facet_args = list())
# {
#   #check_ignored_arguments
#   y <- bayesplot:::validate_y(y)
#   yrep <- bayesplot:::validate_yrep(fit$fitted.values, fit$data$y)
#   group <- bayesplot:::validate_group(fit$data$n, fit$data$y)
#   plot_data <- bayesplot:::ppc_group_data(y, yrep, group)
#   plot_data$t <- rep(seq(1, 50), 1002)
#   plot_data$group2 <- paste(as.character(plot_data$group), as.character(plot_data$variable), sep = "_")
#   is_y <- plot_data$variable == "y"
#
# }
# group2 <- rep(seq(1, 50), 50)
#
# ggplot() +
#   geom_line(data = plot_data %>% dplyr::filter(variable == "y"), aes(x = t, y = value, color = group)) +
#   geom_line(data = plot_data %>% dplyr::filter(grepl("yrep",variable)), aes(x = t, y = value, group = group2, color = group),
#             alpha = 1/10, size = 1/5) +
#       guides(color = guide_legend(title = NULL),
#         fill = guide_legend(order = 1)) + bayesplot_theme_get() +
#         bayesplot:::dont_expand_y_axis() + bayesplot:::no_legend_spacing() + xaxis_title(FALSE) +
#         yaxis_text(FALSE) + yaxis_ticks(FALSE) + yaxis_title(FALSE)
#





bayesplot:::ppc_stat_grouped
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
