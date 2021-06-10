compute_contributions <- function(mu_0, par_control, par_treatment)
{
  N_SIM = 10^6
  z_0 = rnorm(N_SIM, mu_0, 1)
  eps = rnorm(N_SIM)

  compute_mbs_analytically <- function(mu_0, parameters)
  {
    mu_1 = parameters$bias + mu_0
    gamma = parameters$information
    delta = parameters$information + parameters$noise
    var_o = pnorm(mu_0)
    msg = mu_1/sqrt(1-gamma)
    dg = delta/(1 - gamma)
    gg = gamma/sqrt(1 - gamma)

    Sigma_f = matrix(c(dg + 1, dg, dg, dg + 1), ncol = 2)
    var_f = mvtnorm::pmvnorm(upper=c(msg, msg), mean = c(0,0), sigma = Sigma_f)

    Sigma_of = matrix(c(1, gg, gg, dg + 1), ncol = 2)
    cov_of = mvtnorm::pmvnorm(upper=c(mu_0, msg), mean = c(0,0), sigma = Sigma_of)

    as.numeric(var_o + var_f - 2*cov_of)
  }


  compute_mbs_numerically <- function(mu_0, parameters)
  {
    mu_1 = parameters$bias + mu_0
    gamma = parameters$information
    delta = parameters$information + parameters$noise
    outcomes = as.integer(z_0 > 0)
    z_1 = mu_1 + gamma * (z_0 - mu_0) + sqrt(delta - gamma^2) * eps
    probs = pnorm(z_1/sqrt(1 - gamma))
    mean((probs - outcomes)^2)
  }

  compute_mbs = compute_mbs_analytically

  process_sequence <- function(par_sequence)
  {
    prev = par_control
    succ = par_control
    result = list()
    for (par_name in par_sequence)
    {
      succ[[par_name]] = par_treatment[[par_name]]
      contribution = compute_mbs(mu_0, prev) - compute_mbs(mu_0, succ)
      result[[par_name]] = contribution
      prev = succ
    }
    result
  }

  mbs_treatment = compute_mbs(mu_0, par_treatment)
  mbs_control = compute_mbs(mu_0, par_control)

  total = mbs_control - mbs_treatment
  contribution_perms = lapply(combinat::permn(c("bias", "information", "noise")), process_sequence) %>% bind_rows()
  contributions = contribution_perms %>%
    dplyr::summarise_all(mean)
  result = list(
    mean_brier_score_1 = mbs_treatment,
    mean_brier_score_0 = mbs_control,
    contribution_bias = contributions$bias,
    contribution_noise = contributions$noise,
    contribution_information = contributions$information
  )
}

#' Summary
#'
#' This function uses the return value of a call to the function \code{estimate_BIN} and produces a full BIN analysis based on that object.
#'
#' @param full_bayesian_fit The return value of a call to function \code{estimate_BIN}.
#'
#' @return List containing the parameter estimates of the model, the posterior inferences, and the analysis
#' of predictive performance.
#'
#' The elements of the list are as follows.
#' \itemize{
#' \item Parameter Estimates: Posterior means, standard deviations, and different quantiles of the model parameters and their differences.
#' The parameter values represent the following quantities.
#' \itemize{
#' \item mu_star: Base rate of the event outcome in the probit scale. E.g., if mu_star = 0, then, in expectation, Phi(0)*100% = 50% of the events happen, where Phi(.) is the CDF of a standard Gaussian random variable.
#' \item mu_0: The level of bias in the control group. This can be any number in the real-line. E.g., if mu_0 = 0.1, then the control group believes the base rate to be Phi(mu_star + 0.1).
#' \item mu_1: The level of bias in the treatment group; otherwise the interpretation is the same as above for mu_0.
#' \item gamma_0: The level of information in the control group. This is a number between 0 and 1, where 0 represents no information and 1 represents full information.
#' \item gamma_1: The level of information in the treatment group; otherwise the interpretation is the same as above for gamma_0.
#' \item delta_0: The level of noise in the control group. This is a positive value, with higher values indicating higher levels of noise. Noise is in the same scale as information. E.g., delta_0 = 1.0 says that the control group uses as many irrelevant signals as there are relevant signals in the universe. In this sense it represents a very high level of noise.
#' \item delta_1: The level of noise in the treatment group; otherwise the interpretation is the same as above for delta_0.
#' \item rho_0: The level of within-group dependence between forecasts of the control group. This is a positive value, with higher values indicating higher levels of dependence. The dependence can be interpreted to stem from shared irrelevant (noise) and/or relevant (information) signals.
#' \item rho_1: The level of within-group dependence between forecasts of the treatment group; otherwise the interpretation is the same as above for rho_0.
#' \item rho_01: The level of inter-group dependence between forecasts of the control and treatment groups; otherwise the interpretation is the same as above for rho_0.
#' }
#' \item Posterior Inferences: Posterior probabilities of events.
#' Compared to the control group, does the treatment group have: (i) less bias, (ii) more information, and (iii) less noise?
#' Intuitively, one can think of these probabilities as the Bayesian analogs of the p-values in classical hypothesis testing – the closer the probability is to 1, the stronger the evidence for the hypothesis.
#'
#' \item Control,Treatment: This compares the control group against the treatment group.
#' The value of the contribution gives
#' \itemize{
#' \item the mean Brier score of the control group;
#' \item the mean Brier score of the treatment group;
#' \item how the difference can be explained in terms of bias, noise, or information; and
#' \item in percentage terms, how does the change in bias, noise, or information (from control to treatment group) changes the Brier score.
#' }
#' \item Control,Perfect Accuracy: This compares the control group against a treatment group with perfect accuracy;
#' otherwise the interpretation is the same as above for 'Control,Treatment.'
#' }
#'
#'@examples
#' \donttest{
#' ## An example with one group
#' # a) Simulate synthetic data:
#' synthetic_data = simulate_data(list(mu_star = -0.8,mu_0 = -0.5,mu_1 = 0.2,gamma_0 = 0.1,
#' gamma_1 = 0.3,rho_0 = 0.05,delta_0 = 0.1,rho_1 = 0.2, delta_1 = 0.3,rho_01 = 0.05),300,100,0)
#' # b) Estimate the BIN-model on the synthetic data:
#' full_bayesian_fit = estimate_BIN(synthetic_data$Outcomes,synthetic_data$Control, warmup = 500,
#' iter = 1000)
#' # c) Analyze the results:
#' complete_summary(full_bayesian_fit)
#'}
#' \donttest{
#' ## An example with two groups
#' # a) Simulate synthetic data:
#' synthetic_data = simulate_data(list(mu_star = -0.8,mu_0 = -0.5,mu_1 = 0.2,gamma_0 = 0.1,
#' gamma_1 = 0.3, rho_0 = 0.05,delta_0 = 0.1, rho_1 = 0.2, delta_1 = 0.3,rho_01 = 0.05), 300,100,100)
#' # b) Estimate the BIN-model on the synthetic data:
#' full_bayesian_fit = estimate_BIN(synthetic_data$Outcomes,synthetic_data$Control,
#' synthetic_data$Treatment, warmup = 500, iter = 1000)
#' # c) Analyze the results:
#' complete_summary(full_bayesian_fit)
#'}
#'
#' @seealso \code{\link{simulate_data}}, \code{\link{estimate_BIN}}
#'
#' @export
complete_summary<-function(full_bayesian_fit){

  #pacman::p_load("dplyr")
  parameter_name<-NULL

  #List containing the summary of model k
  current_summary=list()

  if(dim(full_bayesian_fit)[3]<20){

    #Result Summary
    result_summary = rstan::summary(full_bayesian_fit)$summary %>%
      as.data.frame() %>% tibble::rownames_to_column("parameter_name") %>%
      dplyr::filter(parameter_name %in% c("mu_star",
                                   "mu_0", "gamma_0",
                                   "delta_0","rho_0")) %>%
      dplyr::mutate_if(is.numeric,function(x) round(x, 2)) %>%
      dplyr::select(-"n_eff",-"Rhat", -"se_mean")

    current_summary[["Parameter Estimates"]]<-result_summary
    posterior_samples <- rstan::extract(full_bayesian_fit)

    sample_values = lapply(posterior_samples, mean)
    par_control <- list(bias = sample_values$mu_0, information = sample_values$gamma_0,
                        noise = sample_values$delta_0)
    par_perfect <- list(bias = 0, information = 1 - 1e-04,
                        noise = 0)
    perfect_mbs_decomp <- compute_contributions(sample_values$mu_star,
                                                par_control, par_perfect)
    p_a_percentage_contributions <- list(perfect_accuracy_percentage_contribution_bias = perfect_mbs_decomp$contribution_bias *
                                           100/perfect_mbs_decomp$mean_brier_score_0, perfect_accuracy_percentage_contribution_noise = perfect_mbs_decomp$contribution_noise *
                                           100/perfect_mbs_decomp$mean_brier_score_0, perfect_accuracy_percentage_contribution_information = perfect_mbs_decomp$contribution_information *
                                           100/perfect_mbs_decomp$mean_brier_score_0)
    contribution_perfect_accuracy = list()
    contribution_perfect_accuracy[["Value of the contribution"]] <- perfect_mbs_decomp
    contribution_perfect_accuracy[["Percentage of control group Brier score"]] <- p_a_percentage_contributions
    current_summary[["Control, Perfect Accuracy"]] <- contribution_perfect_accuracy
  }
  else{

    #Result Summary
    result_summary = rstan::summary(full_bayesian_fit)$summary %>%
      as.data.frame() %>%
      tibble::rownames_to_column("parameter_name") %>%
      dplyr::filter(parameter_name %in% c("mu_star",
                                   "mu_0", "mu_1", "diff_bias",
                                   "gamma_0", "gamma_1", "diff_info",
                                   "delta_0", "delta_1", "diff_noise",
                                   "rho_0", "rho_1", "rho_01")) %>%
      dplyr::mutate_if(is.numeric, function(x) round(x, 2)) %>%
      dplyr::select(-"n_eff",-"Rhat", -"se_mean")

    current_summary[["Parameter Estimates"]]<-result_summary
    posterior_samples <- rstan::extract(full_bayesian_fit)

    ###  Posterior inferences:
    Posterior_inferences<-c("More information in treatment group","Less noise in treatment group","Less bias in treatment group")
    Posterior_Probability<-c(mean(posterior_samples$gamma_0 < posterior_samples$gamma_1),mean(posterior_samples$delta_1 < posterior_samples$delta_0),mean(abs(posterior_samples$bias_1) < abs(posterior_samples$bias_0)))
    Posterior_Inferences<-data.frame(Posterior_inferences,Posterior_Probability)
    current_summary[["Posterior Inferences"]]<-Posterior_Inferences

    ### Predictive Performance & Contributions

    sample_values = lapply(posterior_samples, mean)
    par_control <- list(
      bias = sample_values$mu_0,
      information = sample_values$gamma_0,
      noise = sample_values$delta_0
    )

    par_treatment <- list(
      bias = sample_values$mu_1,
      information = sample_values$gamma_1,
      noise = sample_values$delta_1
    )

    par_perfect <- list(
      bias = 0,
      information = 1 - 1e-4,
      noise = 0
    )

    # Decomposition that compares two groups
    exp_mbs_decomp <- compute_contributions(sample_values$mu_star, par_control, par_treatment)
    # contribution_bias: contribution of bias to expected Brier score
    # contribution_information, contribution_noise: same for information and noise
    # mean_brier_score_0, mean_brier_score_1: expected Brier score for each of the two groups, given corresponding mu, gamma, delta
    #Percentages
    t_percentage_contributions<-list(
      treatment_percentage_contribution_bias=exp_mbs_decomp$contribution_bias*100/exp_mbs_decomp$mean_brier_score_0,
      treatment_percentage_contribution_noise=exp_mbs_decomp$contribution_noise*100/exp_mbs_decomp$mean_brier_score_0,
      treatment_percentage_contribution_information=exp_mbs_decomp$contribution_information*100/exp_mbs_decomp$mean_brier_score_0
    )

    contribution_treatment=list()
    contribution_treatment[["Value of the contribution"]]<-exp_mbs_decomp
    contribution_treatment[["Percentage of control group Brier score"]]<-t_percentage_contributions
    current_summary[["Control,Treatment"]]<-contribution_treatment

    # Decomposition of the maximum achievable improvement
    perfect_mbs_decomp <- compute_contributions(sample_values$mu_star, par_control, par_perfect)
    #Percentages
    p_a_percentage_contributions<-list(
      perfect_accuracy_percentage_contribution_bias=perfect_mbs_decomp$contribution_bias*100/perfect_mbs_decomp$mean_brier_score_0,
      perfect_accuracy_percentage_contribution_noise=perfect_mbs_decomp$contribution_noise*100/perfect_mbs_decomp$mean_brier_score_0,
      perfect_accuracy_percentage_contribution_information=perfect_mbs_decomp$contribution_information*100/perfect_mbs_decomp$mean_brier_score_0
    )

    contribution_perfect_accuracy=list()
    contribution_perfect_accuracy[["Value of the contribution"]]<-perfect_mbs_decomp
    contribution_perfect_accuracy[["Percentage of control group Brier score"]]<-p_a_percentage_contributions
    current_summary[["Control, Perfect Accuracy"]]<-contribution_perfect_accuracy
  }

  return(current_summary)
}
