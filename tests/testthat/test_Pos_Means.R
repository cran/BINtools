
################################################################
# TEST 1: Both groups with many forecasters.

test_that("Post means for mm case match", {
  set.seed(1)
  true_parameters <- list(
    mu_star = -0.8,
    mu_0 = -0.5,
    mu_1 = 0.2,
    gamma_0 = 0.1,
    gamma_1 = 0.3,
    rho_0 = 0.05,
    delta_0 = 0.1,
    rho_1 = 0.2,
    delta_1 = 0.3,
    rho_01 = 0.05
  )

  # Control has many, Treatment has many
  summary_parameters_mm = c("mu_star", "mu_0", "mu_1", "gamma_0", "gamma_1",
                            "rho_0", "delta_0", "rho_1","delta_1","rho_01")
  N = 50 # number of events
  N_0 = 100 # number of control group members
  N_1 = 100 # number of treatment group 1 members
  DATA_mm = simulate_data(true_parameters, N, N_0, N_1)
  # Fit the BIN model
  full_bayesian_fit = estimate_BIN(DATA_mm$Outcomes,DATA_mm$Control,DATA_mm$Treatment,
                                   warmup = 1000, iter = 2000, seed=1)
  # Computer posterior means
  post_means = unlist(lapply(rstan::extract(full_bayesian_fit), mean))
  # Obtained values
  obtained_values=unname(post_means[summary_parameters_mm])
  # Reference values
  reference_values=unlist(unname(true_parameters[summary_parameters_mm]))

  expect_true(all(abs(obtained_values-reference_values)<10))

})

################################################################
# TEST 2.1: One group with many forecasters. One group with one forecaster

test_that("Post means for m1 case match", {
  #skip_on_cran()
  set.seed(1)
  true_parameters <- list(
    mu_star = -0.8,
    mu_0 = -0.5,
    mu_1 = 0.2,
    gamma_0 = 0.1,
    gamma_1 = 0.3,
    rho_0 = 0.05,
    delta_0 = 0.1,
    rho_1 = 0.2,
    delta_1 = 0.3,
    rho_01 = 0.05
  )

  # Control has many, Treatment has one
  summary_parameters_m1 = c("mu_star", "mu_0", "mu_1", "gamma_0", "gamma_1",
                            "rho_0", "delta_0","delta_1","rho_01")
  N = 50 # number of events
  N_0 = 100 # number of control group members
  N_1 = 1 # number of treatment group 1 members
  DATA_m1 = simulate_data(true_parameters, N, N_0, N_1)
  # Fit the BIN model
  full_bayesian_fit = estimate_BIN(DATA_m1$Outcomes,DATA_m1$Control,DATA_m1$Treatment,
                                   warmup = 1000, iter = 2000, seed=1)
  # Computer posterior means
  post_means = unlist(lapply(rstan::extract(full_bayesian_fit), mean))
  # Obtained values
  obtained_values=unname(post_means[summary_parameters_m1])
  # Reference values
  reference_values=unlist(unname(true_parameters[summary_parameters_m1]))

  expect_true(all(abs(obtained_values-reference_values)<10))

})


################################################################
# TEST 2.2: One group with many forecasters.One group with one forecaster.

test_that("Post means for 1m case match", {
  #skip_on_cran()
  set.seed(1)
  true_parameters <- list(
    mu_star = -0.8,
    mu_0 = -0.5,
    mu_1 = 0.2,
    gamma_0 = 0.1,
    gamma_1 = 0.3,
    rho_0 = 0.05,
    delta_0 = 0.1,
    rho_1 = 0.2,
    delta_1 = 0.3,
    rho_01 = 0.05
  )

  # Control has one, Treatment has many
  summary_parameters_1m = c("mu_star", "mu_0", "mu_1", "gamma_0", "gamma_1",
                            "delta_0","delta_1","rho_01")
  N = 50 # number of events
  N_0 = 1 # number of control group members
  N_1 = 100 # number of treatment group 1 members
  DATA_1m = simulate_data(true_parameters, N, N_0, N_1)
  # Fit the BIN model
  full_bayesian_fit = estimate_BIN(DATA_1m$Outcomes,DATA_1m$Control,DATA_1m$Treatment,
                                   warmup = 1000, iter = 2000, seed=1)
  # Computer posterior means
  post_means = unlist(lapply(rstan::extract(full_bayesian_fit), mean))
  # Obtained values
  obtained_values=unname(post_means[summary_parameters_1m])
  # Reference values
  reference_values=unlist(unname(true_parameters[summary_parameters_1m]))

  expect_true(all(abs(obtained_values-reference_values)<10))

})

################################################################
# TEST 3: Both groups with one forecaster.

test_that("Post means for 11 case match", {
  #skip_on_cran()
  set.seed(1)
  true_parameters <- list(
    mu_star = -0.8,
    mu_0 = -0.5,
    mu_1 = 0.2,
    gamma_0 = 0.1,
    gamma_1 = 0.3,
    rho_0 = 0.05,
    delta_0 = 0.1,
    rho_1 = 0.2,
    delta_1 = 0.3,
    rho_01 = 0.05
  )

  summary_parameters_11 = c("mu_star", "mu_0", "mu_1", "gamma_0", "gamma_1",
                            "delta_0", "delta_1","rho_01")
  N = 50 # number of events
  N_0 = 1 # number of control group members
  N_1 = 1 # number of treatment group 1 members
  DATA_11 = simulate_data(true_parameters, N, N_0, N_1)
  # Fit the BIN model
  full_bayesian_fit = estimate_BIN(DATA_11$Outcomes,DATA_11$Control,DATA_11$Treatment,
                                   warmup = 1000, iter = 2000, seed=1)
  # Computer posterior means
  post_means = unlist(lapply(rstan::extract(full_bayesian_fit), mean))
  # Obtained values
  obtained_values=unname(post_means[summary_parameters_11])
  # Reference values
  reference_values=unlist(unname(true_parameters[summary_parameters_11]))

  expect_true(all(abs(obtained_values-reference_values)<10))

})



################################################################
# TEST 4: One group with many forecasters.

test_that("Post means for M0 case match", {
  #skip_on_cran()
  set.seed(1)
  true_parameters <- list(
    mu_star = -0.8,
    mu_0 = -0.5,
    mu_1 = 0.2,
    gamma_0 = 0.1,
    gamma_1 = 0.3,
    rho_0 = 0.05,
    delta_0 = 0.1,
    rho_1 = 0.2,
    delta_1 = 0.3,
    rho_01 = 0.05
  )

  summary_parameters_m = c("mu_star", "mu_0", "gamma_0", "rho_0", "delta_0")
  N = 50 # number of events
  N_0 = 100 # number of control group members
  N_1 = 0 # number of treatment group 1 members
  DATA_m = simulate_data(true_parameters, N, N_0, N_1)
  # Fit the BIN model
  full_bayesian_fit = estimate_BIN(DATA_m$Outcomes,DATA_m$Control, warmup = 1000, iter = 2000,seed=1)
  # Computer posterior means
  post_means = unlist(lapply(rstan::extract(full_bayesian_fit), mean))
  # Obtained values
  obtained_values=unname(post_means[summary_parameters_m])
  # Reference values
  reference_values=unlist(unname(true_parameters[summary_parameters_m]))

  expect_true(all(abs(obtained_values-reference_values)<10))

})

################################################################
# TEST 5: One group with one forecaster.

test_that("Post means for M0 case match", {
  #skip_on_cran()
  set.seed(1)
  true_parameters <- list(
    mu_star = -0.8,
    mu_0 = -0.5,
    mu_1 = 0.2,
    gamma_0 = 0.1,
    gamma_1 = 0.3,
    rho_0 = 0.05,
    delta_0 = 0.1,
    rho_1 = 0.2,
    delta_1 = 0.3,
    rho_01 = 0.05
  )

  summary_parameters_1 = c("mu_star", "mu_0", "gamma_0", "delta_0")
  N = 50 # number of events
  N_0 = 1 # number of control group members
  N_1 = 0 # number of treatment group 1 members
  DATA_1 = simulate_data(true_parameters, N, N_0, N_1)
  # Fit the BIN model
  full_bayesian_fit = estimate_BIN(DATA_1$Outcomes,DATA_1$Control, warmup = 1000, iter = 2000,seed=1)
  # Computer posterior means
  post_means = unlist(lapply(rstan::extract(full_bayesian_fit), mean))
  # Obtained values
  obtained_values=unname(post_means[summary_parameters_1])
  # Reference values
  reference_values=unlist(unname(true_parameters[summary_parameters_1]))

  expect_true(all(abs(obtained_values-reference_values)<10))

})



