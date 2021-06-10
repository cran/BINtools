  // N_0 = 1, N_1 > 1
  functions{
    real information_diversity_probit_likelihood(int N_0, int N_1, int outcome, 
    real M_0, real M_1, real V_0, real V_1, real C_0, real C_1, real C_01, 
    real mu_star, real mu_0, real mu_1, real gamma_0, real gamma_1, real rho_0, real v_0, real v_1, real rho_01){
      
      real D_0 = v_0 + (N_0 - 1) * rho_0;
      real D_1 = v_1;
      real W = D_0 * D_1 - N_0 * N_1 * rho_01 ^ 2;
      real W_0 = W * (rho_0 - v_0);

      // invert the covariance matrix (the inverse will have similar structure)
      real rho_0_i_0 = rho_0 * v_1 - N_1 * rho_01 ^ 2;
      real rho_0_i = rho_0_i_0 / W_0;
      real delta_0_i = (rho_0_i_0 - W) / W_0;
      real delta_1_i = (v_0 + (N_0 - 1)*rho_0) / W;
      real rho_01_i = -1 * rho_01 / W;
      
      // compute the determinant of the covariance matrix
      real logW = log(W);
      real logdet = logW + (N_0 - 1) * log(v_0 - rho_0);
      
      // compute the quadratic form as a function of sufficient statistics
      real cov_0 = delta_0_i * (V_0 - 2 * mu_0 * M_0 + N_0 * mu_0^2);
      real cov_00 = rho_0_i * (C_0 - 2 * mu_0 * (N_0 - 1) * M_0 + N_0 * (N_0 - 1) * mu_0 ^ 2);
      real cov_1 = delta_1_i * (V_1 - 2 * mu_1 * M_1 + N_1 * mu_1^2);
      real cov_01 = 2 * rho_01_i * (C_01 - N_0 * mu_0 * M_1 - N_1 * mu_1 * M_0 + N_0 * N_1 * mu_0 * mu_1);
      real logcov = cov_0 + cov_00 + cov_1 + cov_01;
      
      real D_0_i = delta_0_i + (N_0 - 1) * rho_0_i;
      real D_1_i = delta_1_i;
      
      // compute the probit term: conditional mean and variance
      real mu_star_cond = mu_star + gamma_0 * D_0_i * (M_0 - N_0 * mu_0) + gamma_1 * D_1_i * (M_1 - N_1 * mu_1) + gamma_0 * N_0 * rho_01_i * (M_1 - N_1 * mu_1) + gamma_1 * N_1 * rho_01_i * (M_0 - N_0 * mu_0);
      real ssq_star_cond = 1 - N_0 * gamma_0 ^ 2 * D_0_i - N_1 * gamma_1 ^ 2 * D_1_i - 2 * N_0 * N_1 * gamma_0 * gamma_1 * rho_01_i;
      real sigma_star_cond = sqrt(ssq_star_cond);
      real logprobit = bernoulli_lpmf(outcome | 1 - Phi(-mu_star_cond / sigma_star_cond));
      
      // finally, return the full log-likelihood
      // Jacobian correction for probits is not needed--the Jacobian only depends on data
      return -0.5 * (logdet + logcov + (N_0 + N_1) * log(2 * pi())) + logprobit;
    }
  }
  data {
    int<lower = 1> N; // number of events
    int<lower = 1, upper = 1> N_0[N]; // number of forecasts from control group
    int<lower = 1 > N_1[N]; // number of forecasts from treatment group
    int<lower = 0, upper = 1> outcome[N]; // binary outcome for each event
    real M_0[N]; // sum of probits for control group
    real M_1[N]; // sum of probits for treatment group
    real V_0[N]; // sum of squared probits for control group
    real V_1[N]; // sum of squared probits for treatment group
    real C_0[N]; // sum of out-of-diagonal probit cross-products for control group
    real C_1[N]; // sum of out-of-diagonal probit cross-products for treatment group
    real C_01[N]; // sum of probit cross-products between control and treatment
  }
  parameters{
    real mu_star;
    real mu_0;
    real mu_1;
    real<lower=0, upper = 1> gamma_0;
    real<lower=0, upper = 1> gamma_1;
    real<lower=0> rho_1;
    real<lower=0> delta_0;
    real<lower=0> delta_1;
    real<lower=0> rho_01;
  }
  transformed parameters{
      
      real v_0 = gamma_0 + delta_0;
      real v_1 = gamma_1 + delta_1;
      
      real bias_0 = fabs(mu_0);
      real bias_1 = fabs(mu_1);
  
      real diff_bias = bias_0 - bias_1;
      real diff_info = gamma_0 - gamma_1;
      real diff_noise = delta_0 - delta_1;
        
      // transform the parameters: we observe only Z/sqrt(1-gamma), not Z
      real sg1 = sqrt(1 - gamma_1);
      real sg0 = sqrt(1 - gamma_0);
      real g1 = 1 - gamma_1;
      real g0 = 1 - gamma_0;
      
      real gamma_1_ = gamma_1 / sg1;
      real gamma_0_ = gamma_0 / sg0;
      
      real v_1_ = v_1 / g1;
      real v_0_ = v_0 / g0;
      
      real rho_1_ = rho_1 / g1;
      
      real rho_01_ = rho_01 / (sg1 * sg0);
      real mu_1_ = (mu_star + mu_1) / sg1;
      real mu_0_ = (mu_star + mu_0) / sg0;
  }
  model{
  
  
  for (i in 1:N){
   target += information_diversity_probit_likelihood(N_1[i], N_0[i], outcome[i], M_1[i], M_0[i], V_1[i], V_0[i], C_1[i], C_0[i], C_01[i], 
   mu_star, mu_1_, mu_0_, gamma_1_, gamma_0_, rho_1_, v_1_, v_0_, rho_01_);
  }
  
  }
