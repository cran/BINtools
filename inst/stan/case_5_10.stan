// N_0 = 1, N_1 = 0
  functions{
    real information_diversity_probit_likelihood(int N_0, int N_1, int outcome, 
    real M_0, real M_1, real V_0, real V_1, real C_0, real C_1, real C_01, 
    real mu_star, real mu_0, real gamma_0, real v_0){
      
      real logdet = log(v_0);
      real logcov = (V_0 - 2 * mu_0 * M_0 + mu_0^2) / v_0;

      // compute the probit term: conditional mean and variance
      real mu_star_cond = mu_star +  gamma_0 * (M_0 - mu_0) / v_0;
      real ssq_star_cond = 1 - gamma_0 ^ 2/ v_0;
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
    int<lower = 0, upper = 0> N_1[N]; // number of forecasts from treatment group
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
    real<lower=0, upper = 1> gamma_0;
    real<lower=0> delta_0;
  }
  transformed parameters{
      
      real v_0 = gamma_0 + delta_0;
      
      real bias_0 = fabs(mu_0);
  
      real diff_bias = 0;
      real diff_info = 0;
      real diff_noise = 0;
        
      // transform the parameters: we observe only Z/sqrt(1-gamma), not Z
      real sg0 = sqrt(1 - gamma_0);
      real g0 = 1 - gamma_0;
      real gamma_0_ = gamma_0 / sg0;
      real v_0_ = v_0 / g0;
      real mu_0_ = (mu_star + mu_0) / sg0;
  }
  model{
  
  
  for (i in 1:N){
   target += information_diversity_probit_likelihood(N_0[i], N_1[i], outcome[i], M_0[i], M_1[i], V_0[i], V_1[i], C_0[i], C_1[i], C_01[i], 
   mu_star, mu_0_, gamma_0_, v_0_);
  }
  
  }
