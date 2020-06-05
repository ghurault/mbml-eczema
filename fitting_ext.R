# Notes -------------------------------------------------------------------

# Fitting the extended model (with SWET)

# Initialisation ----------------------------------------------------------

rm(list = ls()) # Clear Workspace

seed <- 1744834965 # seed also used for stan
set.seed(seed) # Reproducibility

library(TanakaData) # Contains data and data processing functions
library(HuraultMisc) # Functions shared across projects
library(rstan)
rstan_options(auto_write = TRUE) # Save compiled model
options(mc.cores = parallel::detectCores()) # Parallel computing
library(ggplot2)
library(cowplot)
source("functions.R") # Additional functions
source("functions_data.R") # Data processing functions

#### OPTIONS
run <- FALSE
n_chains <- 6 # max cores=48
n_it <- 3000
####

stan_code <- "Models/ExtendedModel.stan"
res_file <- "Results/fit_ExtendedModel_SWET.rds"
par_file <- "Results/par_ExtendedModel_SWET.rds"
par0_file <- "Results/par0_ExtendedModel.rds"

main_param <- c("b_S", "sigma_S", "mu_wS", "sigma_wS", "sigma_P")
param_demo <- c("w_FLG", "w_Sex", "w_Age", "w_White", "w_Home")
param_SU <- c("mu_SU", "sigma_SU")
param_CS <- c("mu_CS", "sigma_CS", "w_CS_Mild", "w_CS_Mod", "w_CS_Pot", "w_CS_VPot")
param_CI <- c("mu_CI", "sigma_CI", "w_CI_Mild", "w_CI_Mod")
param_pop <- c(main_param, param_demo, param_SU, param_CS, param_CI)
param_ind <- c("wS", "P", "risk", "w_SU", "w_CS", "b_CS", "w_CI", "b_CI",
               "q_CS_Mild", "q_CS_Mod", "q_CS_Pot", "q_CS_VPot", "q_CI_Mild", "q_CI_Mod")
param_other <- c("S", "S_pred")
param <- c(param_pop, param_ind, param_other)

# Functions ---------------------------------------------------------------

plot_coef <- function(fit, parNames, CI = c(.05, .95), limits = NULL, size = 1) {
  # Plot (population) coefficients
  #
  # Args:
  # fit: stanfit object
  # parNames: vector of names of the parameters to plot
  # CI: (optional) vector of length two indicating the credible interval lower and upper bounds
  # limits: (optional) vector of length 2 for x axis limits
  # size: (optional) size of estimates
  #
  # Returns:
  # Ggplot of coefficient estimates
  
  library(ggplot2)
  
  tmp <- rstan::extract(fit, pars = parNames)
  d <- do.call(rbind,
               lapply(1:length(tmp),
                      function(i) {
                        data.frame(Mean = mean(tmp[[i]]),
                                   Lower = quantile(tmp[[i]], probs = min(CI)),
                                   Upper = quantile(tmp[[i]], probs = max(CI)),
                                   Variable = names(tmp)[i])
                      }))
  
  d$Variable <- sapply(as.character(d$Variable),
                       function(str) {
                         if (substr(str, 1, 2) == "w_") {
                           substr(str, 3, nchar(str)) # remove w prefix
                         }
                       })
  d$Variable <- factor(d$Variable, levels = sort(d$Variable, decreasing = TRUE))
  
  ggplot(data = d, aes(x = Variable, y = Mean, ymin = Lower, ymax = Upper)) +
    geom_pointrange(size = size) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    coord_flip(ylim = limits) +
    labs(y = "Estimate", x = "", title = "") +
    theme_bw(base_size = 20) +
    theme(panel.grid.minor.x = element_blank())
}

plot_responsiveness <- function(fit, patientID) {
  # Plot population and patient responsiveness for each treatment
  #
  # Args:
  # fit: stanfit object
  # patientID: ID of the patient to consider
  #
  # Returns:
  # Ggplot of responsiveness
  
  library(ggplot2)
  
  parNames <- c("b_CS", "b_CI", "w_SU")
  parLabels <- c("Corticosteroids", "Calcineurin Inhibitors", "Step-Up")
  
  sub_fit <- extract(fit, pars = paste0(parNames,"[", which(pt == patientID), "]"))
  pat_distri <- do.call(rbind,
                        lapply(1:length(parNames),
                               function(i) {
                                 d0 <- density(sub_fit[[i]], kernel = "gaussian", n = 128)
                                 data.frame(Value = d0$x,
                                            Density = d0$y,
                                            Parameter = parNames[i],
                                            Distribution = "Patient-level",
                                            Patient = patientID)
                               }))
  
  pop_distri <- do.call("rbind",
                        lapply(1:length(parNames),
                               function(i) {
                                 x <- seq(-2, 2, .001)
                                 tr <- strsplit(parNames[i], "_")[[1]][2]
                                 y <- dnorm(x, mean = par[par$Variable == paste("mu", tr, sep = "_"), "Mean"],
                                            sd = par[par$Variable == paste("sigma", tr, sep = "_"), "Mean"])
                                 data.frame(Value = x,
                                            Density = y,
                                            Parameter = parNames[i],
                                            Distribution = "Population-level",
                                            Patient = NA)
                               }))
  
  tmp <- rbind(pop_distri, pat_distri)
  tmp$Parameter <- factor(tmp$Parameter, levels = parNames, labels = parLabels)
  
  ggplot(data = tmp, aes(x = Value, y = Density, fill = Distribution)) +
    geom_area(alpha = 0.3) +
    geom_line() +
    geom_vline(xintercept = 0, linetype = "dotted", size = 1) +
    facet_grid(rows = vars(Parameter)) +
    scale_fill_manual(values = c("#000000", "#E69F00")) +
    scale_y_continuous(limits = c(0, NA), expand = expand_scale(mult = c(0, .05))) +
    coord_cartesian(xlim = c(-1, 1)) +
    labs(x = "Estimate", fill = "") +
    theme_classic(base_size = 20) +
    theme(legend.position = "top")
}

format_stan_data <- function(dt, dp) {
  list(
    N = length(dt[["Severity"]]),
    N_obs = sum(!is.na(dt[["Severity"]])),
    N_pt = length(unique(dt[["Patient"]])),
    t_max = aggregate(Day ~ Patient, dt, length)[["Day"]],
    idx_obs = which(!is.na(dt[["Severity"]])),
    S_obs = na.omit(dt[["Severity"]]),
    horizon = 0,
    
    FLG = dp[["FLG"]],
    Sex = dp[["Sex"]],
    Age = dp[["Age"]],
    White = dp[["White"]],
    Home = dt[["Home"]],
    Conf = dp[["Confidence"]],
    SU = dt[["SU"]],
    CS = dt[["CS"]],
    CI = dt[["CI"]],
    N_CS = dp[["N_CS"]],
    Q_CS_Mild = dp[["CS_Mild"]],
    Q_CS_Mod = dp[["CS_Mod"]],
    Q_CS_Pot = dp[["CS_Pot"]],
    Q_CS_VPot = dp[["CS_VPot"]],
    N_CI = dp[["N_CI"]],
    Q_CI_Mild = dp[["CI_Mild"]],
    Q_CI_Mod = dp[["CI_Mod"]]
  )
}

# Processing -------------------------------------------------------------------

l <- process2_SWET(SWET)
dt <- l[[1]]
dp <- l[[2]]

pt <- unique(dp[["Patient"]])

if (file.exists(par0_file)) {
  par0 <- readRDS(par0_file)
}

data_stan <- format_stan_data(dt, dp)

# Fit Stan model ----------------------------------------------------------

if (run) {
  fit <- stan(file = stan_code,
              data = data_stan,
              iter = n_it,
              chains = n_chains,
              pars = param,
              seed = seed,
              init = 0,
              control = list(adapt_delta = 0.9))
  saveRDS(fit, file = res_file)
  par <- extract_parameters(fit,
                            param = param,
                            param_ind = param_ind,
                            param_obs = c("S"),
                            param_pred = c("S_pred"),
                            pt = pt,
                            data_stan = data_stan)
  saveRDS(par, file = par_file)
} else {
  fit  <- readRDS(res_file)
  par <- readRDS(par_file)
}

# Results ----------------------------------------------------------------

if (FALSE) {
  
  # Parameters -----------------------------------------------------------------
  
  # shinystan::launch_shinystan(fit) # Shinystan
  check_hmc_diagnostics(fit)
  
  pairs(fit, pars = main_param)
  pairs(fit, pars = param_demo)
  pairs(fit, pars = param_SU)
  pairs(fit, pars = param_CS)
  pairs(fit, pars = param_CI)
  
  # print(fit, pars = param_pop)
  plot(fit, pars = param_pop)
  # HuraultMisc::plot_prior_posterior(par0, par, param_pop)
  # HuraultMisc::check_model_sensitivity(par0, par, c(param_pop, param_ind[!grepl("q_", param_ind)]))
  
  # for (i in param_ind) {print(fit, pars = i)}
  lapply(param_ind, function(x) {plot(fit, pars = x)})
  
  # Distribution density of patient-dependent parameters for single draws
  lapply(param_ind, function(i) {PPC_group_distribution(fit, i, 50)})
  
  # Number of patients for who Prob(theta < 0 | x) > 0.95 ("significance") 
  p <- do.call(rbind,
               lapply(c("w_SU", "w_CS", "b_CS", "w_CI", "b_CI"),
                      function(x) {
                        tmp <- extract(fit, pars = x)[[1]]
                        
                        p_neg <- apply(tmp, 2, function(x) {mean(x < 0)})
                        
                        data.frame(Patient = pt, Probability = p_neg, Parameter = x)
                        # pt[which(apply(tmp, 2, function(x) {mean(x < 0)}) > .95)]
                      }))
  aggregate(Probability ~ Parameter, p, function(x) {sum(x > .95)})
  # Number of treatment applications
  N_treat <- Reduce(function(...) {merge(..., by = "Patient")},
                    lapply(c("SU", "CS", "CI"),
                           function(x) {
                             aggregate(formula(paste0(x, " ~ Patient")), dt, function(y) {sum(y == 1)})
                           }))
  # sum(N_treat$CS > 0)

  # "Significance" of risk factors
  lapply(extract(fit, pars = param_demo), function(x) {mean(x < 0)})
  
  # Plot risk factors and patient responsiveness
  plot_grid(plot_grid(NULL,
                      plot_coef(fit, param_demo, CI = c(.025, .975), limits = c(-.25, .25), size = 1.5),
                      nrow = 2, rel_heights = c(.05, .95)),
            plot_responsiveness(fit, 6006), # can change patient
            nrow = 1, labels = "auto")
  # ggsave("Plots/FigCovariatesSWEText.jpg", width = 20, height = 10, units = "cm", dpi = 300, scale = 2.2)
  
  # Posterior predictive checks -------------------------------------------------------------------
  
  ppc <- prepare_ppc(fit, dt, par, predictions_dictionary(pt, data_stan))
  lapply(sample(pt, 5),
         function(pid) {
           plot_ppc(ppc, patientID = pid)
         })
  
  if (FALSE) {
    # SWET PPC figure
    pl <- lapply(c("1061", "2033", "2098", "6013"),
                 function(pid) {
                   plot_ppc(ppc, patientID = pid)
                 })
    plot_grid(get_legend(pl[[1]] + theme(legend.position = "top")),
              plot_grid(plotlist = lapply(pl, function(x) {x + theme(legend.position = "none")}),
                        ncol = 2,
                        labels = "auto"),
              ncol = 1,
              rel_heights = c(.1, .9))
    # ggsave("Plots/FigPPCSWEText.jpg", width = 30, height = 15, units = "cm", dpi = 300, scale = 1.3)
  }
  
}
