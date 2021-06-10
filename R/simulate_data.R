
transform_parameters <- function(parameters)
{
  #transform the parameters: we observe only Z/sqrt(1-gamma), not Z
  #list2env(parameters, environment())

  gamma_1_p <- parameters$gamma_1 / sqrt(1 - parameters$gamma_1)
  gamma_0_p <- parameters$gamma_0 / sqrt(1 - parameters$gamma_0)

  v_1_p <- (parameters$gamma_1 + parameters$delta_1) / (1 - parameters$gamma_1)
  v_0_p <- (parameters$gamma_0 + parameters$delta_0) / (1 - parameters$gamma_0)

  rho_1_p <- parameters$rho_1 / (1 - parameters$gamma_1)
  rho_0_p <- parameters$rho_0 / (1 - parameters$gamma_0)
  rho_01_p <- parameters$rho_01 / sqrt((1 - parameters$gamma_0) * (1 - parameters$gamma_1))

  mu_1_p <- (parameters$mu_star + parameters$mu_1) / sqrt(1 - parameters$gamma_1)
  mu_0_p <- (parameters$mu_star + parameters$mu_0) / sqrt(1 - parameters$gamma_0)

  list(
    mu_star = parameters$mu_star, mu_0 = mu_0_p, mu_1 = mu_1_p,
    gamma_0 = gamma_0_p, gamma_1 = gamma_1_p,
    v_0 = v_0_p, v_1 = v_1_p,
    rho_0 = rho_0_p, rho_1 = rho_1_p,
    rho_01 = rho_01_p
  )
}

generate_covariance_matrix <- function(parameters, N_0, N_1){
  # generates the block covariance matrix
  #list2env(parameters, environment())

  O <- matrix(1, 1, 1)
  G_0 <- matrix(parameters$gamma_0, 1, N_0)
  G_1 <- matrix(parameters$gamma_1, 1, N_1)
  V_01 <- matrix(parameters$rho_01, N_0, N_1)
  V_0 <- diag(parameters$v_0 - parameters$rho_0, N_0, N_0) + matrix(parameters$rho_0, N_0, N_0)
  V_1 <- diag(parameters$v_1 - parameters$rho_1, N_1, N_1) + matrix(parameters$rho_1, N_1, N_1)
  rbind(cbind(O, G_0, G_1),
        cbind(t(G_0), V_0, V_01),
        cbind(t(G_1), t(V_01), V_1))
}

generate_mean_vector <- function(parameters, N_0, N_1){
  #list2env(parameters, environment())
  c(c(parameters$mu_star), rep(parameters$mu_0, N_0), rep(parameters$mu_1, N_1))
}

simulate_event <- function(z_0, z_sig, mu, sig, N_0, N_1)
{
  # could be made faster by pre-computing Cholesky
  # currently doing this way to avoid errors
  M = N_0 + N_1 + 1
  information_vec = sig[1, 2:M]

  mu_star = mu[1]
  mu_update = information_vec * (z_0 - mu_star)
  sig_update = -1*outer(information_vec, information_vec, '*') # conditioning on the outcome

  forecast_mu = mu[2:M] + mu_update
  forecast_sig = sig[2:M, 2:M] + sig_update

  Z = forecast_mu + chol(forecast_sig) %*% z_sig

  list(
    outcome = as.integer(z_0 > 0),
    ifp_id = stringi::stri_rand_strings(1, 6),
    control_probits = Z[1:N_0],
    treatment_probits = Z[(N_0+1):(N_0+N_1)]
  )
}

simulate_events <- function(parameters, N, N_0, N_1, rho_o = 0.0)
{
  p <- transform_parameters(parameters) # transform the parameters: we observe only Z/sqrt(1-gamma), not Z
  mean_vec <- generate_mean_vector(p, N_0, N_1)
  sigma_mat <- generate_covariance_matrix(p, N_0, N_1)

  Z_o = rnorm(N)
  Z_s = replicate(N, rnorm(N_0 + N_1), simplify=FALSE)

  if (!all(eigen(sigma_mat)$values > 0))
  {
    stop("Covariance matrix is not positive semidefinite.")
  }
  RR = length(rho_o)
  result = vector(mode = "list", length = RR)
  for (j in 1:RR)
  {
    mu_o = mean_vec[1]
    V_o <- diag(1 - rho_o[[j]], N, N) + matrix(rho_o[[j]], N, N)
    Z_os <- mu_o + chol(V_o) %*% Z_o

    events <- vector(mode = "list", length = N)
    for (i in 1:N)
    {
      events[[i]] = simulate_event(Z_os[[i]], Z_s[[i]], mean_vec, sigma_mat, N_0, N_1)
    }
    result[[j]] = events
    attr(result[[j]], "rho_o") = rho_o[[j]]
  }
  result
}

#' Simulate Data
#'
#' This function allows the user to generate synthetic data of two groups (control and treatment) of forecasters making probability predictions of binary events.
#' The function is mostly useful for testing and illustration purposes.
#'
#' @param parameters A list containing the true values of the parameters: mu_star,mu_0,mu_1,gamma_0,gamma_1,rho_0,delta_0,rho_1,delta_1 and rho_01
#' @param N Number of events
#' @param N_0 Number of forecasters in the control group
#' @param N_1 Number of forecasters in the treatment group
#' @param rho_o The level of dependence between event outcomes. (Default: the events are independent conditional on the model parameter values. This sets `rho_ = 0.0`)
#'
#' @details
#' See \code{\link{complete_summary}} for a description of the model parameters.
#' Not all combinations of parameters are possible.
#' In particular, the covariance parameters gamma and rho are dependent on each other and must result in a positive semi-definite covariance matrix for the outcomes and predictions.
#' To find a feasible set of parameters, we recommend users to experiment: begin with the desired levels of mu, gamma, and delta, and values of rho close to zero, and then increase rho until data can be generated without errors.
#'
#' @return List containing the simulated data.
#' The elements of the list are as follows.
#' \itemize{
#' \item Outcomes: Vector containing binary values that indicate the outcome of each event. The j-th entry is equal to 1 if the j-th event occurs and equal to 0 otherwise.
#' \item Control: List of vectors (one for each event) containing probability predictions made by the forecasters in the control group.
#' \item Treatment: List of vectors (one for each event) containing probability predictions made by the forecasters in the treatment group.
#' }
#'
#' @examples
#' \donttest{
#' simulate_data(list(mu_star = -0.8,mu_0 = -0.5,mu_1 = 0.2,gamma_0 = 0.1,gamma_1 = 0.3,
#' rho_0 = 0.05,delta_0 = 0.1,rho_1 = 0.2, delta_1 = 0.3,rho_01 = 0.05), 300,100,100)
#' }
#'
#'
#' @seealso \code{\link{estimate_BIN}}, \code{\link{complete_summary}}
#'
#' @export
simulate_data<-function(parameters,N, N_0, N_1, rho_o = 0.0){
  datasets <- simulate_events(parameters,N, N_0, N_1, rho_o)[[1]]
  control_probab=lapply(datasets, function(ifp) {
    control_prob = pnorm(ifp$control_probits)
    control_prob
  })
  treatment_probab=lapply(datasets, function(ifp) {
    treatment_prob = pnorm(ifp$treatment_probits)
    treatment_prob
  })
  outcomes=unlist(lapply(datasets, function(ifp) {
    outcome= ifp$outcome
    outcome
  }))

  data=list(outcomes,control_probab,treatment_probab)
  names(data)=c("Outcomes", "Control", "Treatment")
  return(data)
}


