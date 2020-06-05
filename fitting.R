# Notes -------------------------------------------------------------------

# Fitting our model (BaseModel) and different baselines (RandomWalk, Autoregression) with Flares and SWET

# Initialisation ----------------------------------------------------------

rm(list = ls()) # Clear Workspace

seed <- 1744834965 # seed also used for stan
set.seed(seed) # Reproducibility

library(TanakaData) # Contains data and data processing functions
library(HuraultMisc) # Functions shared across projects
library(rstan)
rstan_options(auto_write = TRUE) # Save compiled model
options(mc.cores = parallel::detectCores()) # Parallel computing
source("functions.R") # Additional functions
source("functions_data.R") # Data processing functions

#### OPTIONS
mdl_name <- "BaseModel"
dataset <- "SWET"
score <- "Bother"
run <- FALSE
n_chains <- 6
n_it <- 3000
####

mdl_name <- match.arg(mdl_name, c("RandomWalk", "Autoregression", "BaseModel"))
stan_code <- file.path("Models", paste0(mdl_name, ".stan"))

dataset <- match.arg(dataset, c("Flares", "SWET"))
score <- match.arg(score, c("Bother", "Scratch"))
stopifnot(score == "Bother" | dataset == "Flares")

if (mdl_name %in% c("RandomWalk")) {
  param_pop <- c("sigma_S")
  param_ind <- c()
  param_other <- c("S", "S_pred")
} else if (mdl_name == "Autoregression") {
  param_pop <- c("b_S", "sigma_S", "mu_wS", "sigma_wS", "mu_T", "sigma_T")
  param_ind <- c("wS", "wT")
  param_other <- c("S", "S_pred")
} else if (mdl_name == "BaseModel") {
  param_pop <- c("b_S", "sigma_S", "mu_wS", "sigma_wS", "mu_T", "sigma_T", "sigma_P")
  param_ind <- c("wS", "wT", "P")
  param_other <- c("S", "S_pred")
}
param <- c(param_pop, param_ind, param_other)

suff <- paste0(mdl_name, "_", dataset, "_", score, ".rds")
res_file <- file.path("Results", paste0("fit_", suff))
par_file <- file.path("Results", paste0("par_", suff))
par0_file <- file.path("Results", paste0("par0_", mdl_name, ".rds"))

# Functions ---------------------------------------------------------------

plot_patient_coef <- function(fit, parNames, pt, CI = c(.05, .95)) {
  # Plot patient coefficient estimates from stan model (custom function).
  # Patient are ordered according to the first parameter in parNames.
  #
  # Args:
  # fit: stanfit object
  # parNames: vector of names of the patient-dependent parameter to plot
  # pt: vector of patient ID (same order as the patients in the model)
  # CI: (optional) vector of length two indicating the credible interval lower and upper bounds
  #
  # Returns:
  # List of ggplot of patient coefficient estimates
  
  library(ggplot2)
  
  # Extract and summarise posterior
  tmp <- rstan::extract(fit, pars = parNames)
  d <- do.call(rbind,
               lapply(1:length(tmp),
                      function(i) {
                        data.frame(Patient = factor(pt, levels = rev(pt)),
                                   Mean = apply(tmp[[i]], 2, mean),
                                   Lower = apply(tmp[[i]], 2, function(x) {quantile(x, probs = min(CI))}),
                                   Upper = apply(tmp[[i]], 2, function(x) {quantile(x, probs = max(CI))}),
                                   Variable = names(tmp)[i])
                      }))
  
  # Order patients by the mean value of parNames[1]
  par1 <- subset(d, Variable == parNames[1])
  d$Patient <- factor(d$Patient, levels = par1$Patient[order(par1$Mean, decreasing = TRUE)])
  
  # Plot
  lapply(1:length(parNames),
         function(i) {
           ggplot(data = subset(d, Variable == parNames[i]),
                  aes(x = Patient, y = Mean, ymin = Lower, ymax = Upper)) +
             geom_pointrange() +
             coord_flip() +
             theme_bw(base_size = 20) + 
             theme(panel.grid.minor.x = element_blank(),
                   axis.text.y = element_blank())
         })
}

format_stan_data <- function(df) {
  with(df,
       list(N = length(Severity),
            N_obs = sum(!is.na(Severity)),
            N_pt = length(unique(Patient)),
            t_max = aggregate(Day ~ Patient, FUN = length)$Day,
            idx_obs = which(!is.na(Severity)),
            S_obs = na.omit(Severity),
            Treat = Treatment, # not used in RandomWalk
            horizon = 0))
}

# Processing -------------------------------------------------------------------

if (dataset == "Flares") {
  df <- process_Flares(load_Flares(), score)
} else if (dataset == "SWET") {
  df <- process1_SWET(SWET)
}

if (file.exists(par0_file)) {
  par0 <- readRDS(par0_file)
}

pt <- unique(df[["Patient"]])

data_stan <- format_stan_data(df)

# Fit Stan model ----------------------------------------------------------

if (run) {
  fit <- stan(file = stan_code,
              data = data_stan,
              iter = n_it,
              chains = n_chains,
              pars = param,
              seed = seed,
              # init = 0,
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
  fit <- readRDS(res_file)
  par <- readRDS(par_file)
}

# Results ----------------------------------------------------------------

if (FALSE) {
  
  ## Diagnostics
  
  # shinystan::launch_shinystan(fit) # Shinystan
  check_hmc_diagnostics(fit)
  pairs(fit, pars = param_pop)
  plot(fit, pars = param_pop, plotfun = "trace")
  
  ## Population parameters
  
  # print(fit, pars = param_pop) # might take some time, results already in par
  # plot(fit, pars = param_pop)
  HuraultMisc::plot_prior_posterior(par0, par, param_pop)
  HuraultMisc::check_model_sensitivity(par0, par, c(param_pop, param_ind))
  
  ## Patient-dependent parameters
  if (FALSE) {
    # Figure patient-dependent parameters
    pl <- plot_patient_coef(fit, c("wS", "wT", "P"), pt)
    cowplot::plot_grid(pl[[1]] + labs(y = expression(paste("Persistence (", w[S]^(k), ")", sep = ""))) + coord_flip(ylim = c(0, 1)),
                       pl[[2]] + labs(y = expression(paste("Responsiveness to treatment (", w[T]^(k), ")", sep = ""))) + coord_flip(ylim = c(-1.5, 1.5)),
                       pl[[3]] + labs(y = expression(paste("Flare triggers (", P^(k), ")", sep = ""))) + coord_flip(ylim = c(0, 5.5)),
                       nrow = 1, labels = "auto")
    # ggsave("Plots/FigPatientCoefFlaresBother.jpg", width = 20, height = 15, units = "cm", dpi = 300, scale = 1.8)
    # ggsave("Plots/FigPatientCoefSWET.jpg", width = 20, height = 15, units = "cm", dpi = 300, scale = 2.2)

    # Distribution density of patient-dependent parameters for single draws
    lapply(param_ind, function(i) {PPC_group_distribution(fit, i, 50)})
  }
  
  ## Posterior predictive checks
  ppc <- prepare_ppc(fit, df, par, predictions_dictionary(pt, data_stan))
  lapply(sample(pt, 5),
         function(pid) {
           plot_ppc(ppc, patientID = pid)
         })
  
  if (FALSE) {
    # Flares PPC figure
    pl <- lapply(c(15, 13, 35, 52),
                 function(pid) {
                   plot_ppc(ppc, patientID = pid)
                 })
    plot_grid(get_legend(pl[[1]] + theme(legend.position = "top")),
              plot_grid(plotlist = lapply(pl, function(x) {x + theme(legend.position = "none")}),
                        ncol = 2,
                        labels = "auto"),
              ncol = 1,
              rel_heights = c(.1, .9))
    # ggsave("Plots/FigPPCFlaresBother.jpg", width = 30, height = 15, units = "cm", dpi = 300, scale = 1.3)
  }

}
