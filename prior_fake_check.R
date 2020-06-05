# Notes -------------------------------------------------------------------

# Prior predictive checks and fake data check

# Fake data are generated as missing values to be constrained between 0 and 10, but still need to be rounded
# To check the prior distribution, you should specify a short time series (e.g. n_dur = 5) and 1 patient (n_pt = 1) is enough
# This is because for longer time series, since the score is bounded by 0 and 10, some areas of the priors cannot be accessed easily due to numerical constraints
# This limit the predictive distribution available for fake data check (but it's better than nothing)
# Before running the fake data check, check that the sampling of the prior predictive distribution is OK
# When testing on different samples from the prior predictive distribution we can recover the true parameters and patient parameters have a good coverage probability

# Initialisation ----------------------------------------------------------

rm(list = ls()) # Clear workspace (but bette to restart R entirely)

seed <- 1744834965 # seed also used for Stan
set.seed(seed) # Reproducibility

library(HuraultMisc) # Functions shared across projects
library(rstan)
rstan_options(auto_write = TRUE) # Save compiled model
options(mc.cores = parallel::detectCores()) # Parallel computing
library(ggplot2)
library(cowplot)
source("functions.R") # Additional functions

#### OPTIONS
mdl_name <- "ExtendedModel"
n_pt <-  1 # 1, 20 # Number of fake patients
n_dur <- 5 # 5, 100 # Time series length
prop_missing <- 0.3 # Proportion of missing values
run_prior <- TRUE
run_fake <- FALSE
n_chains <- 4
n_it <- 2000
####

mdl_name <- match.arg(mdl_name, c("RandomWalk", "Autoregression", "BaseModel", "ExtendedModel"))
stan_code <- file.path("Models", paste0(mdl_name, ".stan"))

if (mdl_name == "RandomWalk") {
  param_pop <- c("sigma_S")
  param_ind <- c()
  param_other <- c("S")
} else if (mdl_name == "Autoregression") {
  param_pop <- c("b_S", "sigma_S", "mu_wS", "sigma_wS", "mu_T", "sigma_T")
  param_ind <- c("wS", "wT")
  param_other <- c("S")
} else if (mdl_name == "BaseModel") {
  param_pop <- c("b_S", "sigma_S", "mu_wS", "sigma_wS", "mu_T", "sigma_T", "sigma_P")
  param_ind <- c("wS", "wT", "P")
  param_other <- c("S")
} else if (mdl_name %in% c("ExtendedModel")) {
  main_param <- c("b_S", "sigma_S", "mu_wS", "sigma_wS", "sigma_P")
  param_demo <- c("w_FLG", "w_Sex", "w_Age", "w_White", "w_Home")
  param_SU <- c("mu_SU", "sigma_SU")
  param_CS <- c("mu_CS", "sigma_CS", "w_CS_Mild", "w_CS_Mod", "w_CS_Pot", "w_CS_VPot")
  param_CI <- c("mu_CI", "sigma_CI", "w_CI_Mild", "w_CI_Mod")
  param_pop <- c(main_param, param_demo, param_SU, param_CS, param_CI)
  param_ind <- c("wS", "P", "risk", "w_SU", "w_CS", "b_CS", "w_CI", "b_CI",
                 "q_CS_Mild", "q_CS_Mod", "q_CS_Pot", "q_CS_VPot", "q_CI_Mild", "q_CI_Mod")
  param_other <- c("S")
}
param <- c(param_pop, param_ind, param_other)

prior_file <- file.path("Results", paste0("prior_", mdl_name, ".rds"))
par0_file <- file.path("Results", paste0("par0_", mdl_name, ".rds"))
fake_file <- file.path("Results", paste0("fake_", mdl_name, ".rds"))

if (any(c(run_prior, run_fake))) {
  compiled_model <- stan_model(stan_code)
}

# Processing --------------------------------------------------------------

# Data
pt <- 1:n_pt
df <- expand.grid(Patient = pt, Day = 1:n_dur)
df <- df[order(df$Patient, df$Day), ]
df[["S"]] <- NA
df[df$Day == 1, "S"] <- sample(0:10, n_pt, replace = TRUE)
if (mdl_name != "ExtendedModel") {
  df[["Treat"]] <- do.call(c,
                           lapply(1:n_pt,
                                  function(x) {
                                    generate_treatment(c(rbeta(1, 2, 3), rbeta(1, 3, 2)), n_dur)
                                  }))
  
  format_stan_data <- function(df, lbl) {
    list(
      N = length(df[[lbl]]),
      N_obs = sum(!is.na(df[[lbl]])),
      N_pt = length(unique(df[["Patient"]])),
      t_max = array(aggregate(Day ~ Patient, df, length)$Day),
      idx_obs = array(which(!is.na(df[[lbl]]))),
      S_obs = array(na.omit(df[[lbl]])),
      Treat = array(df[["Treat"]]),
      horizon = 0
    )
  }
  
}

if (mdl_name == "ExtendedModel") {
  for (x in c("CS", "CI", "SU", "Home")) {
    df[[x]] <- do.call(c,
                       lapply(1:n_pt,
                              function(x) {
                                generate_treatment(c(rbeta(1, 2, 3), rbeta(1, 3, 2)), n_dur)
                              }))
  }
  
  dp <- data.frame(
    Patient = 1:n_pt,
    Age = round(abs(rnorm(n_pt, 0, 5))),
    Sex = rbinom(n_pt, 1, .5), 
    White = rbinom(n_pt, 1, .5),
    FLG = rbinom(n_pt, 1, .5),
    N_CS = array(sample(1:n_dur, n_pt, replace = TRUE)),
    N_CI = array(sample(1:n_dur, n_pt, replace = TRUE)),
    Confidence = sample(1:4, n_pt, replace = TRUE)
  )
  dp <- cbind(dp,
              data.frame(
                # number of applications * using treatment (0-1) * daily quantity
                CS_Mild = dp[["N_CS"]] * rbinom(n_pt, 1, .25) * round(abs(rnorm(n_pt, 0, 2)), 1), # total appli
                CS_Mod = dp[["N_CS"]] * rbinom(n_pt, 1, .25) * round(abs(rnorm(n_pt, 0, 2)), 1),
                CS_Pot = dp[["N_CS"]] * rbinom(n_pt, 1, .25) * round(abs(rnorm(n_pt, 0, 2)), 1),
                CS_VPot = dp[["N_CS"]] * rbinom(n_pt, 1, .25) * round(abs(rnorm(n_pt, 0, 2)), 1),
                CI_Mild = dp[["N_CI"]] * rbinom(n_pt, 1, .25) * round(abs(rnorm(n_pt, 0, 2)), 1),
                CI_Mod = dp[["N_CI"]] * rbinom(n_pt, 1, .25) * round(abs(rnorm(n_pt, 0, 2)), 1)
              )
  )
  
  format_stan_data <- function(dt, dp, lbl) {
    list(
      N = length(dt[[lbl]]),
      N_obs = sum(!is.na(dt[[lbl]])),
      N_pt = length(unique(dt[["Patient"]])),
      t_max = array(aggregate(Day ~ Patient, dt, length)[["Day"]]),
      idx_obs = array(which(!is.na(dt[[lbl]]))),
      S_obs = array(na.omit(dt[[lbl]])),
      horizon = 0,
      
      FLG = array(dp[["FLG"]]),
      Sex = array(dp[["Sex"]]),
      Age = array(dp[["Age"]]),
      White = array(dp[["White"]]),
      Home = array(dt[["Home"]]),
      Conf = array(dp[["Confidence"]]),
      SU = dt[["SU"]],
      CS = dt[["CS"]],
      CI = dt[["CI"]],
      N_CS = array(dp[["N_CS"]]),
      Q_CS_Mild = array(dp[["CS_Mild"]]),
      Q_CS_Mod = array(dp[["CS_Mod"]]),
      Q_CS_Pot = array(dp[["CS_Pot"]]),
      Q_CS_VPot = array(dp[["CS_VPot"]]),
      N_CI = array(dp[["N_CI"]]),
      Q_CI_Mild = array(dp[["CI_Mild"]]),
      Q_CI_Mod = array(dp[["CI_Mod"]])
    )
  }
  
}

# Prior predictive checks -----------------------------------------------------

if (mdl_name != "ExtendedModel") {
  data_prior <- format_stan_data(df, "S")
} else {
  data_prior <- format_stan_data(df, dp, "S")
}

if (run_prior) {
  fit_prior <- sampling(compiled_model,
                        data = data_prior,
                        pars = param,
                        iter = n_it,
                        chains = n_chains,
                        seed = seed,
                        control = list(adapt_delta = .9))
  saveRDS(fit_prior, file = prior_file)
  par0 <- extract_parameters(fit_prior,
                             param = param,
                             param_ind = param_ind,
                             param_obs = c("S"),
                             param_pred = c(),
                             pt = pt,
                             data_stan = data_prior)
  saveRDS(par0, file = par0_file)
} else {
  fit_prior <- readRDS(prior_file)
  par0 <- readRDS(par0_file)
}

# Analyse results
if (FALSE) {
  check_hmc_diagnostics(fit_prior)
  pairs(fit_prior, pars = param_pop)
  
  # Distribution of parameters
  plot(fit_prior, pars = param_pop)
  plot(fit_prior, pars = c(param_pop, paste0(param_ind, "[1]")), plotfun = "hist")
  
  # Posterior predictive distribution
  lapply(pt[1:min(length(pt), 5)],
         function(i) {
           ggplot(data = subset(par0, Patient == i & Variable == "S"),
                  aes(x = Day, y = Mean, ymin = `5%`, ymax = `95%`)) +
             geom_line() +
             geom_ribbon(alpha = .5) +
             scale_y_continuous(breaks = 0:10, limits = c(0, 10)) +
             theme_bw(base_size = 20) +
             theme(panel.grid.minor.y = element_blank())
         })
}

# Fake data check -------------------------------------------------------

s_meas <- extract(fit_prior, pars = "S")[[1]]

draw <- sample(1:nrow(s_meas), 1) # Take one draw from predictive distribution

# Fake trajectory
df[["S_fake"]] <- s_meas[draw, ]
df[["S_fake"]] <- round(df[["S_fake"]]) # round
df[as.logical(rbinom(nrow(df), 1, prop_missing)), "S_fake"] <- NA # Missing values

# Extract true parameters values
true_param <- HuraultMisc::extract_parameters_from_draw(fit_prior, param, draw)
true_param[["Patient"]] <- NA
id <- (true_param[["Parameter"]] %in% param_ind)
true_param[id, "Patient"] <- pt[true_param[id, "Index"]]

# Look at the data
lapply(pt[1:5],
       function(patientID) {
         ggplot(data = subset(df, Patient == patientID),
                aes(x = Day, y = S_fake)) +
           geom_path() +
           scale_y_continuous(limits = c(0, 10), breaks = 0:10) +
           labs(y = "Severity (fake)") +
           theme_bw(base_size = 15) +
           theme(panel.grid.minor.y = element_blank())
       })

# Fit model with fake data
if (mdl_name != "ExtendedModel") {
  data_fake <- format_stan_data(df, "S_fake")
} else {
  data_fake <- format_stan_data(df, dp, "S_fake")
}
param <- c(param, "S_pred")
if (run_fake) {
  fit_fake <- sampling(compiled_model,
                       data = data_fake,
                       pars = param,
                       iter = n_it,
                       chains = n_chains,
                       seed = seed,
                       control = list(adapt_delta = 0.9))
  saveRDS(fit_fake, file = fake_file)
} else {
  fit_fake <- readRDS(fake_file)
}

# Analyse results
if (FALSE) {
  
  check_hmc_diagnostics(fit_fake)
  pairs(fit_fake, pars = param_pop)
  
  par_fake <- extract_parameters(fit_fake,
                                 param = param,
                                 param_ind = param_ind,
                                 param_obs = c(),
                                 param_pred = c("S_pred"),
                                 pt = pt,
                                 data_stan = data_fake)
  
  # Can we recover population parameters
  tmp <- merge(subset(par_fake, Variable %in% c(param_pop, param_ind)),
               change_colnames(true_param, c("Parameter", "Value"), c("Variable", "True")),
               by = c("Variable", "Patient"))
  tmp$Patient <- factor(tmp$Patient, levels = pt)
  ggplot(data = subset(tmp, Variable %in% param_pop),
         aes(x = Variable)) +
    geom_pointrange(aes(y = Mean, ymin = `5%`, ymax = `95%`)) +
    geom_point(aes(y = True), colour = "#E69F00", size = 2) +
    coord_flip() +
    labs(x = "", y = "Estimate") +
    theme_bw(base_size = 15)
  
  # Can we recover patient parameters
  lapply(intersect(c("wS", "wT", "b_CS", "b_CI", "w_SU", "P"), param_ind),
         function(var_name) {
           # Coefficient plot
           tmp <- subset(tmp, Variable == var_name)
           tmp$Patient <- factor(tmp$Patient, levels = tmp[order(tmp$Mean), "Patient"])
           p1 <- ggplot(data = tmp,
                        aes(x = Patient)) +
             geom_pointrange(aes(y = Mean, ymin = `5%`, ymax = `95%`)) +
             geom_point(aes(y = True), colour = "#E69F00") +
             coord_flip() +
             labs(x = "", y = "Estimate") +
             theme_bw(base_size = 15)
           
           # Coverage plot
           p2 <- HuraultMisc::plot_coverage(extract(fit_fake, pars = var_name)[[1]],
                                            true_param[true_param[["Parameter"]] == var_name, "Value"])
           
           plot_grid(p1, p2, ncol = 2)
         })
  
  # Posterior predictive checks
  ppc <- prepare_ppc(fit_fake, change_colnames(df, "S_fake", "Severity"), par_fake, predictions_dictionary(pt, data_fake))
  lapply(sample(pt, 5),
         function(pid) {
           plot_ppc(ppc, patientID = pid)
         })
  
}
