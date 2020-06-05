data {
  int<lower = 0> N; // Total number of observations (missing and non-missing)
  int<lower = 0> N_obs; // Number of non-missing observations
  int<lower = 0> N_pt; // Number of patients
  int<lower = 0> t_max[N_pt]; // Vector of time-series length (number of days) for each patient 
  int<lower = 1, upper = N> idx_obs[N_obs]; // Index of non-missing observations
  real<lower = 0, upper = 10> S_obs[N_obs]; // Observed severity score
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
  real<lower = 0> sigma_S; // Standard deviation of the Gaussian
}

transformed parameters {
  real S[N];
  S[idx_obs] = S_obs;
  S[idx_mis] = S_mis;
}

model {
  sigma_S ~ normal(0, 1.5);
  
  for (k in 1:N_pt) {
    // Loop over patients
    for (t in (start[k] + 1):end[k]) {
      // Loop over time
      S[t] ~ normal(S[t - 1], sigma_S) T[0, 10]; 
    }
  }
  
}

generated quantities {
  vector[N_pred] S_pred;
  
  {
    int i = 1; // Indexing S[t]
    int i_pred = 1; // Indexing S_pred[t]
    real S_prev; // S[t - 1]
    
    for (k in 1:N_pt) {
      S_pred[i_pred] = S[i]; // Initialisation
      for (t in 2:(t_max[k] + horizon)) {
        i_pred += 1;
        if (t <= t_max[k]) {
          i += 1;
          S_prev = S[i - 1]; // Fit
        } else if (t == (t_max[k] + 1)) {
          S_prev = S[i]; // First prediction
        } else {
          S_prev = S_pred[i_pred - 1]; // Remaining predictions
        }
        S_pred[i_pred] = normal_rng(S_prev, sigma_S);
      }
      i_pred += 1;
      i += 1;
    }
  }
  
}
