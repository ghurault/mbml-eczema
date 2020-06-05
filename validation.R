# Notes -------------------------------------------------------------------

# Run the forward chaining

# Initialisation ----------------------------------------------------------

rm(list = ls()) # Clear Workspace

seed <- 1744834965 # seed also used for stan
set.seed(seed) # Reproducibility

library(TanakaData) # Contains data and data processing functions
library(HuraultMisc) # Functions shared across projects
library(rstan)
rstan_options(auto_write = TRUE) # Save compiled model
options(mc.cores = parallel::detectCores()) # Parallel computing
library(foreach)
library(doParallel)
source("functions.R") # Additional functions
source("functions_data.R")

#### OPTIONS
mdl_name <- "BaseModel"
dataset <- "Flares"
score <- "Bother"
t_horizon <- 7
run <- FALSE
n_chains <- 6 # max cores=48
n_it <- 3000
n_cluster <- 6 # floor(parallel::detectCores() / n_chains)
####

mdl_name <- match.arg(mdl_name, c("Uniform", "Historical", "RandomWalk", "Autoregression", "BaseModel"))
stan_code <- file.path("Models", paste0(mdl_name, ".stan"))

dataset <- match.arg(dataset, c("Flares", "SWET"))
score <- match.arg(score, c("Bother", "Scratch"))
stopifnot(score == "Bother" | dataset == "Flares")

suff <- paste(mdl_name, dataset, score, sep = "_")
dir_name <- file.path("Results", paste0("val_", suff)) # temporary directory
res_file <- file.path("Results", paste0("val_", suff, ".rds"))
cal_file <- file.path("Results", paste0("cal_", suff, ".rds"))

if (mdl_name == "RandomWalk") {
  param_pop <- c("sigma_S")
  param_ind <- c()
  param_other <- c("S_mis", "S_pred")
} else if (mdl_name == "Autoregression") {
  param_pop <- c("b_S", "sigma_S", "mu_wS", "sigma_wS", "mu_T", "sigma_T")
  param_ind <- c("wS", "wT")
  param_other <- c("S_mis", "S_pred")
} else if (mdl_name == "BaseModel") {
  param_pop <- c("b_S", "sigma_S", "mu_wS", "sigma_wS", "mu_T", "sigma_T", "sigma_P")
  param_ind <- c("wS", "wT", "P")
  param_other <- c("S_mis", "S_pred")
}

is_stan_model <- !(mdl_name %in% c("Uniform", "Historical"))

# Functions ---------------------------------------------------------------

process_predictions <- function(fit, df_test) {
  # Process predictions:
  # - Identify test replications
  # - Compute probability table from samples
  # - Compute lpd and RPS
  #
  # Args:
  # fit: Stanfit object
  # df_test: Test dataset
  #
  # Returns:
  # Prediction dataframe
  
  do.call(rbind,
          lapply(unique(df_test[["Patient"]]),
                 # Deal patient by patient as need to remove wrong trajectories
                 function(pid) {
                   tmp <- subset(df_test, Patient == pid)
                   lbl <- paste0("S_pred[", tmp[["Index"]], "]")
                   ps <- rstan::extract(fit, pars = lbl)
                   # Put in matrix form
                   if (length(ps) > 1) {
                     ps <- do.call(cbind, ps)
                   } else {
                     ps <- matrix(ps[[1]], ncol = 1)
                   }
                   # Compute probability table
                   prob <- compute_pmf(ps, pred = TRUE)
                   # Join probability table
                   colnames(prob) <- paste0("P(S=", colnames(prob), ")")
                   tmp <- cbind(tmp, prob)
                   # Compute lpd, RPS and expected value
                   tmp <- summarise_predictions(tmp)
                   return(tmp)
                 }))
  
}

summarise_predictions <- function(res) {
  # Add lpd, RPS and E(S) column to prediction dataframe
  #
  # Args:
  # res: Prediction dataframe
  #
  # Returns:
  # Prediction dataframe
  
  prob <- prediction_matrices(res)$Forecast
  for (i in 1:nrow(res)) {
    res$lpd[i] <- log(prob[i, res$Severity[i] + 1])
    res$RPS[i] <- compute_RPS(as.numeric(prob[i, ]), res$Severity[i] + 1)
  }
  res[["E(S)"]] <- prob %*% (0:10) # Expected value
  return(res)
}

# Processing -------------------------------------------------------------------

if (is_stan_model) {
  param <- c(param_pop, param_ind, param_other)
  
  format_stan_data <- function(df) {
    with(df,
         list(N = length(Severity),
              N_obs = sum(!is.na(Severity)),
              N_pt = length(unique(Patient)),
              t_max = aggregate(Day ~ Patient, FUN = length)$Day,
              idx_obs = which(!is.na(Severity)),
              S_obs = na.omit(Severity),
              Treat = Treatment, # not used in RandomWalk
              horizon = t_horizon))
  }
  
  if (run) {
    compiled_model <- rstan::stan_model(stan_code)
  }
}

if (dataset == "Flares") {
  df <- process_Flares(load_Flares(), score)
} else if (dataset == "SWET") {
  df <- process1_SWET(SWET)
}

# df <- subset(df, Patient %in% unique(df$Patient)[1:20])

max_day <- aggregate(Day ~ Patient, df, max)
max_it <- floor((max(max_day$Day) - 1) / t_horizon) # -1 so that there is at least one prediction in the last iteration

# Forward chaining --------------------------------------------------------

if (run) {
  
  duration <- Sys.time()
  cl <- makeCluster(n_cluster)
  registerDoParallel(cl)
  
  writeLines(c(""), "log.txt")
  dir.create(dir_name)
  
  out <- foreach(it = max_it:0) %dopar% {
    # Need to reload functions and libraries
    library(rstan)
    rstan_options(auto_write = TRUE) # Save compiled model
    options(mc.cores = parallel::detectCores()) # Parallel computing
    source("functions.R") # need to reload functions
    
    sink("log.txt", append = TRUE)
    cat(paste("Starting model", it, "\n"))
    
    ####
    
    # Hold out data
    train_days <- 1:(it * t_horizon + 1)
    test_days <- (it * t_horizon + 1):((it + 1) * t_horizon) + 1
    df_train <- df[df$Day %in% train_days, ]
    pt <- unique(df_train[["Patient"]])
    df_test <- na.omit(df[df$Day %in% test_days, c("Patient", "Day", "Severity")])
    
    # Compute prediction horizon
    last_obs <- aggregate(Day ~ Patient, df_train[!is.na(df_train$Severity), ], max) # Last observed value
    colnames(last_obs)[colnames(last_obs) == "Day"] <- "Last_obs"
    df_test <- merge(df_test, last_obs, by = "Patient", all.x = TRUE, all.y = FALSE)
    df_test[["Horizon"]] <- df_test[["Day"]] - df_test[["Last_obs"]]
    df_test[["Last_obs"]] <- NULL
    
    if (is_stan_model) {
      # Fit Stan model
      data_stan <- format_stan_data(df_train)
      fit <- sampling(compiled_model,
                      data = data_stan,
                      iter = n_it,
                      chains = n_chains,
                      # init = 0,
                      pars = param)
      
      # Get index of test predictions
      df_test <- merge(df_test, predictions_dictionary(pt, data_stan), all.x = TRUE, all.y = FALSE)
      
      # Extract and process predictions
      tmp_pred <- process_predictions(fit, df_test)
      tmp_pred$Iteration <- it
      
      # Extract parameters
      tmp_par <- extract_parameters(fit,
                                    param = param,
                                    param_ind = param_ind,
                                    param_obs = c("S"),
                                    param_pred = c("S_pred"),
                                    pt = pt,
                                    data_stan = data_stan)
      tmp_par$Iteration <- it
      
    } else if (mdl_name == "Uniform") {
      prob <- matrix(1 / 11, ncol = 11, nrow = nrow(df_test))
      colnames(prob) <- paste0("P(S=", 0:10, ")")
      tmp_pred <- cbind(df_test, prob)
      tmp_pred <- summarise_predictions(tmp_pred)
      tmp_pred$Iteration <- it

      tmp_par <- NULL
    } else if (mdl_name == "Historical") {
      p <- table(c(0:10, df_train$Severity))
      p <- p / sum(p)
      prob <- matrix(rep(p, nrow(df_test)), ncol = length(p), byrow = TRUE)
      colnames(prob) <- paste0("P(S=", 0:10, ")")
      tmp_pred <- cbind(df_test, prob)
      tmp_pred <- summarise_predictions(tmp_pred)
      tmp_pred$Iteration <- it
      
      tmp_par <- NULL
    }
    
    # Save results (better to save in the loop in case something breaks)
    saveRDS(list(Prediction = tmp_pred, Parameters = tmp_par),
            file = file.path(dir_name, paste0("val_", it, ".rds")))
    
    ####
    
    cat(paste("Ending model", it, "\n"))
    NULL # return
  }
  stopCluster(cl)
  (duration = Sys.time() - duration)
  
  # Recombine results
  files <- list.files(dir_name)
  if (length(files) < max_it + 1) {
    warning("Number of files (", length(files), ") less than the number of iterations (", max_it + 1, "). Some runs may have failed.")
  }
  res_parallel <- lapply(files,
                         function(f) {
                           readRDS(file.path(dir_name, f))
                         })
  res <- do.call("rbind", lapply(res_parallel, function(x) {x$Prediction}))
  par <- do.call("rbind", lapply(res_parallel, function(x) {x$Parameters}))
  saveRDS(list(Prediction = res, Parameters = par), file = res_file)
  
}
