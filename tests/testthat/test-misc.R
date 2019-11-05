# pars_include
test_that("Correct parameters are included", {
  expect_equal(pars_include(Mx_a = 1, Mx_b = 1, Mx_c = 1, Mx_d = 1, Mx_e = 1) , c())
  expect_equal(pars_include(Mx_a = 10, Mx_b = 10, Mx_c = 0, Mx_d = 0, Mx_e = 30) , c("lambda", "alpha"))
  expect_equal(pars_include(Mx_a = 0, Mx_b = 1, Mx_c = 1, Mx_d = 1, Mx_e = 1) , c("gamma"))
  expect_equal(pars_include() , c("gamma", "delta", "lambda", "alpha", "beta"))
})

# naming_fun
test_that("Named output of correct length", {
  x <- c("pi1[1]", "pi1[2]", "A[1]", "alpha", "sigma", rep("logalpha[1]"), 10)
  para_names <- list(init_prob = c(0.1, 0.1), tp = c(0.2), alpha = c(8), sigma = c(0.2))
  names(para_names$init_prob) <- c("pi1[1]", "pi1[2]"); names(para_names$tp) <- c("A[1]"); names(para_names$alpha) <- c("shared_par"); names(para_names$sigma) <- c("sigma")
  expect_equal(length(naming_fun(x = x, para_names = para_names)), length(x))
  expect_equal(naming_fun(x = x, para_names = para_names)[6:16], x[6:16])
  expect_equal(naming_fun(x = x, para_names = para_names)[4], names(para_names$alpha))
})

# naming_states
test_that("Names are state-specific", {
  expect_equal(naming_state(c("x_e1", "x_e2"), 2), c("S1_x_e1", "S2_x_e1", "S1_x_e2", "S2_x_e2"))
})

# check n and t
test_that("N and T check", {
  names_dta <- c("a", "n", "t"); n_var <- "N"; t_var <- "t"
  expect_error(.check_nt(names_dta = names_dta, n_var = n_var, t_var = t_var), "Please supply indexes for n and t.")
})

# inclusion check
test_that("Inclusion check", {
  names_dta <- c("var1", "var2", "var3", "var4"); parsed_formula <- list(all.var = c("var5", "var4"))
  expect_error(.check_inclusion(names_dta = names_dta, parsed_formula = parsed_formula), "The predictor var5 has not been found in the data.")
})

# intercept check
test_that("Intercept detected", {
  names_dta <- "Intercept"; parsed_formula <- list(has_intercept = c(0, 1))
  expect_error(.check_intercept(names_dta = names_dta, parsed_formula = parsed_formula), "Please specify the intercept via 1 +.")
})

# create_order_vector
test_that("Creation of the order vector", {
  formula <- list(e = NULL, has_intercept = c(0, 0)); order_continuous = NULL
  expect_equal(length(create_order_vector(formula, order_continuous)), 0)
  formula <- list(e = rep("a", 20), has_intercept = c(0, 1)); order_continuous = NULL
  expect_equal(length(create_order_vector(formula, order_continuous)), 21)
  formula <- list(e = rep("Intercept", 2), has_intercept = c(0, 1)); order_continuous = c("Intercept")
  expect_equal(sum(create_order_vector(formula, order_continuous)), 1)
})



# two intercepts in continuous?


