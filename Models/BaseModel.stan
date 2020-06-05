functions {
  real soft_uniform_lpdf(real x, real lb, real ub) {
    return(log(inv_logit(x - lb) - inv_logit(x - ub)) - log(ub - lb));
  }
}

data {
  int<lower = 0> N; // Total number of observations (missing and non-missing)
  int<lower = 0> N_obs; // Number of non-missing observations
  int<lower = 0> N_pt; // Number of patients
  int<lower = 0> t_max[N_pt]; // Vector of time-series length (number of days) for each patient 
  int<lower = 1, upper = N> idx_obs[N_obs]; // Index of non-missing observations
  real<lower = 0, upper = 10> S_obs[N_obs]; // Observed severity score
  real<lower = 0, upper = 1> Treat[N]; // Daily treatment usage
  int<lower = 0> horizon; // Time horizon (in days) for prediction
}

transformed data {
  int N_pred = N + N_pt * horizon; // Number of observations for posterior predictive check (fit + prediction)
  int start[N_pt]; // Index of first observation for patient each patient
  int end[N_pt]; // Index of last observation for patient each patient
  int N_mis = N - N_obs; // Number of missing observations
  int idx_mis[N_mis]; // Index of missing observations
  
  if (N != sum(t_max)) {
    reject("N should be equal to sum(t_max)")
  }
  
  // Start and end of each time-series
  for (k in 1:N_pt) {
    if (k == 1) {
      start[k] = 1;
    } else {
      start[k] = end[k - 1] + 1;
    }
    end[k] = start[k] - 1 + t_max[k];
  }
  
  // Index of missing observations
  {
    int id = 1;
    int obs[N] = rep_array(0, N);
    obs[idx_obs] = rep_array(1, N_obs);
    for (i in 1:N) {
      if (obs[i] == 0) {
        idx_mis[id] = i;
        id += 1;
      }
    }
  }
  
}

parameters {
  real<lower = 0, upper = 10> S_mis[N_mis]; // Missing S
  real<lower = -0.5, upper = 0.5> err[N_obs]; // Rounding error
  real<lower = 0> sigma_S; // Standard deviation of the Gaussian
  real b_S; // Intercept
  
  real mu_wS; // Population autocorrelation logit mean
  real<lower = 0> sigma_wS; // Population autocorrelation logit standard deviation
  real eta_wS[N_pt]; // Non-centered parametrisation for autocorrelation
  
  real mu_T; // Population mean for responsiveness to treatment
  real<lower = 0> sigma_T; // Population standard deviation for responsiveness to treatment
  real eta_T[N_pt]; // Non-centered parametrisation for responsiveness to treatment
  
  real<lower = 0> sigma_P; // Population standard deviation for P
  real<lower = 0> eta_P[N_pt]; // Non-centered parametrisation for log P
  real<lower = 0> eta_R[N]; // Non-centered parametrisation for R
}

transformed parameters {
  real S[N]; // Latent severity (before rounding)
  real wS[N_pt]; // Patient autocorrelation
  real wT[N_pt]; // Patient responsiveness to treatment
  real P[N_pt]; // Scale of distribution for R (pathogens load)
  real R[N]; // Flare intensity
  
  // Define S: mix observe and missing values, rounding process
  for (i in 1:N_obs) {
    if (S_obs[i] == 0) {
      S[idx_obs[i]] = S_obs[i] + (0.25 + 0.5 * err[i]); //cf. bounds at 0
    } else if (S_obs[i] == 10) {
      S[idx_obs[i]] = S_obs[i] - (0.25 + 0.5 * err[i]); // cf. bounds at 10
    } else {
      S[idx_obs[i]] = S_obs[i] + err[i];
    }
  }
  S[idx_mis] = S_mis;
  
  for (k in 1:N_pt) {
    wS[k] = inv_logit(mu_wS + sigma_wS * eta_wS[k]);
    wT[k] = mu_T + sigma_T * eta_T[k];
    P[k] = sigma_P * eta_P[k];
    for (t in start[k]:end[k]) {
      R[t] = P[k] * eta_R[t];
    }
  }

}

model {
  eta_wS ~ std_normal();
  eta_T ~ std_normal();
  eta_P ~ std_normal();
  eta_R ~ exponential(1);
  
  b_S ~ normal(0, 2);
  sigma_S ~ normal(0, 1.5);
  mu_wS ~ normal(0, 1);
  sigma_wS ~ normal(0, 1.5);
  mu_T ~ normal(0, 1);
  sigma_T ~ normal(0, 0.5);
  sigma_P ~ normal(0, 1);
  
  {
    real mu;
    for (k in 1:N_pt) {
      // Loop over patients
      (b_S + wT[k]) ~ normal(0, 2); // prior on "constant term" (can't be too big)
      for (t in (start[k] + 1):end[k]) {
        // Loop over time
        mu = wS[k] * S[t - 1] + wT[k] * Treat[t - 1] + b_S + R[t];
        mu ~ soft_uniform(-1, 11); // Regularising prior ensuring mu mostly in 0-10
        S[t] ~ normal(mu, sigma_S) T[0, 10]; 
      }
    }
  }
  
}

generated quantities {
  vector[N_pred] S_pred;
  
  {
    int i = 1; // Indexing S[t]
    int i_pred = 1; // Indexing S_pred[t]
    real S_prev; // S[t - 1]
    real T_prev; // Treat[t - 1]
    
    for (k in 1:N_pt) {
      S_pred[i_pred] = S[i]; // Initialisation
      for (t in 2:(t_max[k] + horizon)) {
        i_pred += 1;
        if (t <= t_max[k]) {
          i += 1;
          S_prev = S[i - 1]; // Fit
          T_prev = Treat[i - 1];
        } else if (t == (t_max[k] + 1)) {
          S_prev = S[i]; // First prediction
          T_prev = Treat[i];
        } else {
          S_prev = S_pred[i_pred - 1]; // Remaining predictions
          T_prev = 0; // Assume no treatment
        }
        S_pred[i_pred] = exp_mod_normal_rng(wS[k] * S_prev + wT[k] * T_prev + b_S, sigma_S, 1 / P[k]);
      }
      i_pred += 1;
      i += 1;
    }
  }
  
}
