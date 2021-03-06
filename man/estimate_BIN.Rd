% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_BIN.R
\name{estimate_BIN}
\alias{estimate_BIN}
\title{Estimate a BIN (Bias, Information, Noise) model}
\usage{
estimate_BIN(
  Outcomes,
  Control,
  Treatment = NULL,
  initial = list(mu_star = 0, mu_0 = 0, mu_1 = 0, gamma_0 = 0.4, gamma_1 = 0.4, delta_0
    = 0.5, rho_0 = 0.27, delta_1 = 0.5, rho_1 = 0.27, rho_01 = 0.1),
  warmup = 2000,
  iter = 4000,
  seed = 1
)
}
\arguments{
\item{Outcomes}{Vector of binary values indicating the outcome of each event. The j-th entry is equal to 1 if the j-th event occurs and equal to 0 otherwise.}

\item{Control}{List of vectors containing the predictions made for each event by forecasters in the control group. The j-th vector contains predictions for the j-th event.}

\item{Treatment}{(Default:\code{NULL}) List of vectors containing the predictions made for each event by forecasters in the treatment group. The j-th vector contains predictions for the j-th event.}

\item{initial}{A list containing the initial values for the parameters mu_star,mu_0,mu_1,gamma_0,gamma_1,delta_0,rho_0,delta_1,rho_1,and rho_01.
(Default: \code{list(mu_star = 0,mu_0 = 0,mu_1 = 0,gamma_0 = 0.4,gamma_1 = 0.4,
delta_0 = 0.5,rho_0 = 0.27, delta_1 = 0.5,rho_1 = 0.27,rho_01 = 0.1)})}

\item{warmup}{The number of initial iterations used for ``burnin.''
These values are not included in the analysis of the model. (Default:\code{2000})}

\item{iter}{Total number of iterations.
Must be larger than warmup. (Default:\code{4000})}

\item{seed}{(Default: \code{1})}
}
\value{
Model estimation is performed with the statistical programming language called \href{https://mc-stan.org/}{Stan}.
The return object is a Stan model.
This way the user can apply available diagnostics tools in other packages, such as \href{https://mc-stan.org/rstan/}{rstan}, to analyze the final results.
}
\description{
This function allows the user to compare two groups (treatment and control) of forecasters in terms of their bias, information, and noise levels.
Model estimation is performed with a Markov Chain Monte Carlo (MCMC) approach called Hamiltonian Monte Carlo.
}
\examples{
\donttest{
## An example with one group
# a) Simulate synthetic data:
synthetic_data = simulate_data(list(mu_star = -0.8,mu_0 = -0.5,mu_1 = 0.2,gamma_0 = 0.1,
gamma_1 = 0.3,rho_0 = 0.05,delta_0 = 0.1,rho_1 = 0.2, delta_1 = 0.3,rho_01 = 0.05),300,100,0)
# b) Estimate the BIN-model on the synthetic data:
full_bayesian_fit = estimate_BIN(synthetic_data$Outcomes,synthetic_data$Control, warmup = 500,
iter = 1000)
# c) Analyze the results:
complete_summary(full_bayesian_fit)
}
\donttest{
## An example with two groups
# a) Simulate synthetic data:
synthetic_data = simulate_data(list(mu_star = -0.8,mu_0 = -0.5,mu_1 = 0.2,gamma_0 = 0.1,
gamma_1 = 0.3, rho_0 = 0.05,delta_0 = 0.1, rho_1 = 0.2, delta_1 = 0.3,rho_01 = 0.05), 300,100,100)
# b) Estimate the BIN-model on the synthetic data:
full_bayesian_fit = estimate_BIN(synthetic_data$Outcomes,synthetic_data$Control,
synthetic_data$Treatment, warmup = 500, iter = 1000)
# c) Analyze the results:
complete_summary(full_bayesian_fit)
}

}
\seealso{
\code{\link{simulate_data}}, \code{\link{complete_summary}}
}
