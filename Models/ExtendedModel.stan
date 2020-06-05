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
  int<lower = 0> horizon; // Time horizon (in days) for prediction
  
  real<lower = 0, upper = 1> FLG[N_pt]; // Presence of filaggrin mutation
  real<lower = 0, upper = 1> Sex[N_pt]; // Sex
  real<lower = 0> Age[N_pt]; // Age
  real<lower = 0, upper = 1> White[N_pt]; // Ethnicity (white skin)
  real<lower = 0, upper = 1> Home[N]; // Sleeping at home
  real<lower = 1, upper = 4> Conf[N_pt]; // Confidence in the reported treatment quantity
  real<lower = 0, upper = 1> SU[N]; // Daily step-up (SU) usage
  real<lower = 0, upper = 1> CS[N]; // Daily corticosteroids (CS) usage
  real<lower = 0, upper = 1> CI[N]; // Daily calcineurin inhibitors (CI) usage
  real<lower = 0> N_CS[N_pt]; // Number of days CS were used during the first 12 weeks
  real<lower = 0> Q_CS_Mild[N_pt]; // Total quantity of mild CS during the first 12 weeks
  real<lower = 0> Q_CS_Mod[N_pt]; // Total quantity of moderate CS during the first 12 weeks
  real<lower = 0> Q_CS_Pot[N_pt]; // Total quantity of potent CS during the first 12 weeks
  real<lower = 0> Q_CS_VPot[N_pt]; // Total quantity of very potent CS during the first 12 weeks
  real<lower = 0> N_CI[N_pt]; // Number of days CI were used during the first 12 weeks
  real<lower = 0> Q_CI_Mild[N_pt]; // Total quantity of mild CI during the first 12 weeks
  real<lower = 0> Q_CI_Mod[N_pt]; // Total quantity of moderate CI
}

transformed data {
  int N_pred = N + N_pt * horizon; // Number of observations for posterior predictive check (fit + prediction)
  int start[N_pt]; // Index of first observation for patient each patient
  int end[N_pt]; // Index of last observation for patient each patient
  int N_mis = N - N_obs; // Number of missing observations
  int idx_mis[N_mis]; // Index of missing observations
  int bin[16, 4]; // Interactions between binary variables (used to define priors)
  real sigma_Q = 0.25; // Standard deviation for treatment quantities estimates (set values since not identifiable as parameter)
  
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
  
  // Dummy variables (decimal to binary)
  {
    int ct = 0;
    for (i in 1:16) {
      ct = i - 1;
      for (j in 1:4) {
        bin[i, j] = ct % 2;
        ct = ct / 2;
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
  
  real<lower = 0> sigma_P; // Population standard deviation for P
  real<lower = 0> eta_P[N_pt]; // Non-centered parametrisation for log P
  real<lower = 0> eta_R[N]; // Non-centered parametrisation for R
  
  // Population risk factors
  real w_FLG; // Filaggrin
  real w_Sex; // Sex
  real w_Age; // Age
  real w_White; // Ethnicity (white skin)
  real w_Home; // Sleeping at home
  
  // Treatments
  // Step-Up
  real mu_SU; // Population mean responsiveness to SU
  real<lower = 0> sigma_SU; // Population standard deviation responsiveness to SU
  real eta_SU[N_pt]; // cf. non-centered parametrisation for SU responsiveness
  // Corticosteroids
  real mu_CS; // Population mean for the intrinsic responsiveness to CS
  real<lower = 0> sigma_CS; // Population standard deviation for the intrinsic responsiveness to CS
  real eta_CS[N_pt]; // cf. non-centered parametrisation for CS intrinsic responsiveness
  real eta_CS_Mild[N_pt]; // cf. non-centered parametrisation for the quantity of mild CS
  real eta_CS_Mod[N_pt]; // cf. non-centered parametrisation for the quantity of moderate CS
  real eta_CS_Pot[N_pt]; // cf. non-centered parametrisation for the quantity of potent CS
  real eta_CS_VPot[N_pt]; // cf. non-centered parametrisation for the quantity of very potent CS
  real w_CS_Mild; // coefficient for mild CS on patient responsiveness
  real w_CS_Mod; // coefficient for moderate CS on patient responsiveness
  real w_CS_Pot; // coefficient for potent CS on patient responsiveness
  real w_CS_VPot; // coefficient for very potent CS on patient responsiveness
  // Calcineurin inhibitors
  real mu_CI; // Population mean for the intrinsic responsiveness to CI
  real<lower = 0> sigma_CI; // Population standard deviation for the intrinsic responsiveness to CI
  real eta_CI[N_pt]; // cf. non-centered parametrisation for CI intrinsic responsiveness
  real eta_CI_Mild[N_pt]; // cf. non-centered parametrisation for the quantity of mild CI
  real eta_CI_Mod[N_pt]; // cf. non-centered parametrisation for the quantity of moderate CI
  real w_CI_Mild; // coefficient for mild CI on patient responsiveness
  real w_CI_Mod; // coefficient for moderate CI on patient responsiveness
  
}

transformed parameters {
  real S[N]; // Latent severity (before rounding)
  real wS[N_pt]; // Patient autocorrelation
  real P[N_pt]; // Scale of distribution for R (pathogens load)
  real R[N]; // Flare intensity
  real risk[N_pt]; // Sum of patient-dependent risk factors
  
  // Step-up
  real w_SU[N_pt]; // Patient responsiveness to SU
  // Corticosteroids
  real w_CS[N_pt]; // Patient responsiveness to CS
  real b_CS[N_pt]; // Patient intrinsic responsiveness to CS
  real q_CS_Mild[N_pt]; // Estimated average daily quantity of mild CS
  real q_CS_Mod[N_pt]; // Estimated average daily quantity of moderate CS
  real q_CS_Pot[N_pt]; // Estimated average daily quantity of potent CS
  real q_CS_VPot[N_pt]; // Estimated average daily quantity of very potent CS
  // Calcineurin inhibitors
  real w_CI[N_pt]; // Patient responsiveness to CI
  real b_CI[N_pt]; // Patient intrinsic responsiveness to CI
  real q_CI_Mild[N_pt]; // Estimated average daily quantity of mild CI
  real q_CI_Mod[N_pt]; // Estimated average daily quantity of moderate CI
  
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
    // Estimated daily quantities of CS
    if (Q_CS_Mild[k] == 0 || N_CS[k] == 0) {
      q_CS_Mild[k] = 0;
    } else {
      q_CS_Mild[k] = Q_CS_Mild[k] * exp(eta_CS_Mild[k] * sigma_Q / sqrt(Conf[k])) / N_CS[k];
    }
    if (Q_CS_Mod[k] == 0 || N_CS[k] == 0) {
      q_CS_Mod[k] = 0;
    } else {
      q_CS_Mod[k] = Q_CS_Mod[k] * exp(eta_CS_Mod[k] * sigma_Q / sqrt(Conf[k])) / N_CS[k];
    }
    if (Q_CS_Pot[k] == 0 || N_CS[k] == 0) {
      q_CS_Pot[k] = 0;
    } else {
      q_CS_Pot[k] = Q_CS_Pot[k] * exp(eta_CS_Pot[k] * sigma_Q / sqrt(Conf[k])) / N_CS[k];
    }
    if (Q_CS_VPot[k] == 0 || N_CS[k] == 0) {
      q_CS_VPot[k] = 0;
    } else {
      q_CS_VPot[k] = Q_CS_VPot[k] * exp(eta_CS_VPot[k] * sigma_Q / sqrt(Conf[k])) / N_CS[k];
    }
    // Estimated daily quantities of CI
    if (Q_CI_Mild[k] == 0 || N_CI[k] == 0) {
      q_CI_Mild[k] = 0;
    } else {
      q_CI_Mild[k] = Q_CI_Mild[k] * exp(eta_CI_Mild[k] * sigma_Q / sqrt(Conf[k])) / N_CI[k];
    }
    if (Q_CI_Mod[k] == 0 || N_CI[k] == 0) {
      q_CI_Mod[k] = 0;
    } else {
      q_CI_Mod[k] = Q_CI_Mod[k] * exp(eta_CI_Mod[k] * sigma_Q / sqrt(Conf[k])) / N_CI[k];
    }
    // Responsiveness
    w_SU[k] = mu_SU + eta_SU[k] * sigma_SU;
    b_CS[k] = mu_CS + eta_CS[k] * sigma_CS;
    w_CS[k] = b_CS[k] +
              w_CS_Mild * q_CS_Mild[k] +
              w_CS_Mod * q_CS_Mod[k] +
              w_CS_Pot * q_CS_Pot[k] +
              w_CS_VPot * q_CS_VPot[k];
    b_CI[k] = mu_CI + eta_CI[k] * sigma_CI;
    w_CI[k] = b_CI[k] +
              w_CI_Mild * q_CI_Mild[k] +
              w_CI_Mod * q_CI_Mod[k];
    
    risk[k] = w_FLG * FLG[k] + w_Sex * Sex[k] + w_Age * Age[k] + w_White * White[k];
    wS[k] = inv_logit(mu_wS + sigma_wS * eta_wS[k]);
    P[k] = sigma_P * eta_P[k];
    for (t in start[k]:end[k]) {
      R[t] = P[k] * eta_R[t];
    }
  }

}

model {
  b_S ~ normal(0, 2);
  sigma_S ~ normal(0, 1.5);
  mu_wS ~ normal(0, 1);
  sigma_wS ~ normal(0, 1.5);
  sigma_P ~ normal(0, 1);
  
  // Prior for the risk factors
  to_array_1d({w_FLG, w_Sex, w_White, w_Home}) ~ normal(0, 0.5);
  w_Age ~ normal(0, 0.1);
  // Prior for treatment parameters
  to_array_1d({mu_SU, mu_CS, mu_CI}) ~ normal(0, 1);
  to_array_1d({sigma_SU, sigma_CS, sigma_CI}) ~ normal(0, 0.5); // gamma(1.8, 4.5); // cf. 10% of the mass below 0.1 and 5% above 1 (because identifiability issue when sigma_Q >> sigma_CS)
  to_array_1d({w_CS_Mild, w_CS_Mod, w_CS_Pot, w_CS_VPot, w_CI_Mild, w_CI_Mod}) ~ normal(0, 0.5);
  
  // Non-centered parametrisation
  eta_wS ~ std_normal();
  eta_P ~ std_normal();
  eta_R ~ exponential(1);
  eta_SU ~ std_normal();
  eta_CS ~ std_normal();
  eta_CI ~ std_normal();
  eta_CS_Mild ~ std_normal();
  eta_CS_Mod ~ std_normal();
  eta_CS_Pot ~ std_normal();
  eta_CS_VPot ~ std_normal();
  eta_CI_Mild ~ std_normal();
  eta_CI_Mod ~ std_normal();
  
  {
    real mu;
    for (k in 1:N_pt) {
      // Loop over patients
      for (j in 1:16) {
        // Prior on constant term (same as prior for intercept)
        b_S + risk[k] + w_SU[k] * bin[j, 1] + w_CS[k] * bin[j, 2] + w_CI[k] * bin[j, 3] + w_Home * bin[j, 4] ~ normal(0, 2);
      }
      for (t in (start[k] + 1):end[k]) {
        // Loop over time
        mu =  wS[k] * S[t - 1] + b_S + R[t] +
              w_Home * Home[t] + w_SU[k] * SU[t - 1] + w_CS[k] * CS[t - 1] + w_CI[k] * CI[t - 1] + risk[k];
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
    real Home_prev; // Home[t - 1]
    real SU_prev; // SU[t - 1]
    real CS_prev; // CS[t - 1]
    real CI_prev; // CI[t - 1]
    
    for (k in 1:N_pt) {
      S_pred[i_pred] = S[i]; // Initialisation
      for (t in 2:(t_max[k] + horizon)) {
        i_pred += 1;
        if (t <= t_max[k]) {
          i += 1;
          S_prev = S[i - 1]; // Fit
          Home_prev = Home[i - 1];
          SU_prev = SU[i - 1];
          CS_prev = CS[i - 1];
          CI_prev = CI[i - 1];
        } else if (t == (t_max[k] + 1)) {
          S_prev = S[i]; // First prediction
          Home_prev = Home[i];
          SU_prev = SU[i];
          CS_prev = CS[i];
          CI_prev = CI[i];
        } else {
          S_prev = S_pred[i_pred - 1]; // Remaining predictions
          Home_prev = 1; // Assume staying at home
          SU_prev = 0; // Assume no step-up
          CS_prev = 0; // Same
          CI_prev = 0; // Same
        }
        S_pred[i_pred] = exp_mod_normal_rng(wS[k] * S_prev + b_S + risk[k] + w_Home * Home_prev +
                                            w_SU[k] * SU_prev + w_CS[k] * CS_prev + w_CI[k] * CI_prev,
                                            sigma_S, 1 / P[k]);
      }
      i_pred += 1;
      i += 1;
    }
  }
  
}
