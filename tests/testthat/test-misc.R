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
