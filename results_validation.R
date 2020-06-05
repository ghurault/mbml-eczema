# Notes -------------------------------------------------------------------

# Analyse validation results

# Initialisation ----------------------------------------------------------

rm(list = ls()) # Clear Workspace

seed <- 1744834965 # seed also used for stan
set.seed(seed) # Reproducibility

library(TanakaData) # Contains data and data processing functions
library(HuraultMisc) # Functions shared across projects
library(ggplot2)
library(cowplot)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
source("functions.R") # Additional functions
source("functions_data.R") # Data processing functions

#### OPTIONS
dataset <- "Flares"
score <- "Bother"
mdl_names <- c("Uniform", "Historical", "RandomWalk", "Autoregression", "BaseModel")
new_names <- c("Uniform", "Historical", "Random Walk", "Autoregression", "Our Model")
####

dataset <- match.arg(dataset, c("Flares", "SWET"))
score <- match.arg(score, c("Bother", "Scratch"))
stopifnot(score == "Bother" | dataset == "Flares")

res_files <- file.path("Results", paste0("val_", mdl_names, "_", dataset, "_", score, ".rds"))

# Processing --------------------------------------------------------------------

if (dataset == "Flares") {
  df <- process_Flares(load_Flares(), score)
} else if (dataset == "SWET") {
  df <- process1_SWET(SWET)
}

pt <- unique(df[["Patient"]])

# Learning curve ----------------------------------------------------------

RPS <- do.call(rbind,
               lapply(1:length(mdl_names),
                      function(i) {
                        res <- readRDS(res_files[i])$Prediction
                        perf <- estimate_performance(res, metric = "RPS")$Fit
                        perf[["Model"]] <- mdl_names[i]
                        return(perf)
                      }))
# RPS$Model <- factor(RPS$Model, levels = rev(mdl_names))
RPS$Model <- factor(RPS$Model, levels = rev(mdl_names), labels = rev(new_names))

if (FALSE) {
  # RPS comparison
  tmp <- reshape2::dcast(subset(RPS, Horizon == 1),
                         Iteration ~ Model,
                         value.var = "Mean")
  with(tmp, (BaseModel - Uniform) / Uniform)
}

p1 <- ggplot(data = subset(RPS, Horizon == 1),
             aes(x = Iteration, y = Mean, ymin = Mean - SE, ymax = Mean + SE, colour = Model, fill = Model)) +
  geom_line(size = 1.5) +
  geom_ribbon(alpha = .5) +
  scale_colour_manual(values = cbbPalette) +
  scale_fill_manual(values = cbbPalette) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Training week", y = "RPS", colour = "", fill = "") +
  theme_bw(base_size = 15) +
  theme(legend.position = "top")
p1
# saveRDS(p1, file = file.path("Plots", paste0("RPS_", dataset, "_", score, ".rds")))

lpd <- do.call(rbind,
               lapply(1:length(mdl_names),
                      function(i) {
                        res <- readRDS(res_files[i])$Prediction
                        perf <- estimate_performance(res, metric = "lpd")$Fit
                        perf[["Model"]] <- mdl_names[i]
                        return(perf)
                      }))
lpd$Model <- factor(lpd$Model, levels = rev(mdl_names), labels = rev(new_names))

brk <- c(.1, .2, .3)
p2 <- ggplot(data = subset(lpd, Horizon == 1),
             aes(x = Iteration, y = Mean, ymin = Mean - SE, ymax = Mean + SE, colour = Model, fill = Model)) +
  geom_line(size = 1.5) +
  geom_ribbon(alpha = .5) +
  scale_colour_manual(values = cbbPalette) +
  scale_fill_manual(values = cbbPalette) +
  scale_y_continuous(limits = c(NA, log(max(brk))), breaks = log(brk), labels = paste0("log(", brk, ")")) +
  labs(x = "Training week", y = "Log predictive probability", colour = "", fill = "") +
  theme_bw(base_size = 15) +
  theme(legend.position = "top")
p2
# saveRDS(p2, file = file.path("Plots", paste0("lpd_", dataset, "_", score, ".rds")))

if (FALSE) {
  
  AUROC <- do.call(rbind,
                   lapply(1:length(mdl_names),
                          function(i) {
                            res <- readRDS(res_files[i])$Prediction
                            if (mdl_names[i] %in% c("Uniform", "Historical")) {
                              perf <- data.frame(Iteration = unique(res$Iteration),
                                                 Horizon = 1, # doesn't matter for now
                                                 Fit = 0.5,
                                                 SE = 0)
                            } else {
                              perf <- estimate_discrimination(res)$Fit
                            }
                            perf[["Model"]] <- mdl_names[i]
                            return(perf)
                          }))
  AUROC$Model <- factor(AUROC$Model, levels = rev(mdl_names), labels = rev(new_names))
  
  p3 <- ggplot(data = subset(AUROC, Horizon == 1),
         aes(x = Iteration, y = Fit, ymin = Fit - SE, ymax = Fit + SE, colour = Model, fill = Model)) +
    geom_line(size = 1.5) +
    geom_ribbon(alpha = .5) +
    scale_colour_manual(values = cbbPalette) +
    scale_fill_manual(values = cbbPalette) +
    coord_cartesian(ylim = c(0.5, 1)) +
    labs(x = "Training week", y = "AUROC", colour = "", fill = "") +
    theme_bw(base_size = 20) +
    theme(legend.position = "top")
  p3
  # saveRDS(p3, file = file.path("Plots", paste0("AUROC_", dataset, "_", score, ".rds")))
  
  MAE <- do.call(rbind,
                 lapply(1:length(mdl_names),
                        function(i) {
                          res <- readRDS(res_files[i])$Prediction
                          res$error <- abs(res$`E(S)` - res$Severity) # Absolute error for MAE or squared error for MSE
                          perf <- estimate_performance(res, metric = "error")$Fit
                          perf[["Model"]] <- mdl_names[i]
                          return(perf)
                        }))
  MAE$Model <- factor(MAE$Model, levels = rev(mdl_names), labels = rev(new_names))
  
  ggplot(data = subset(MAE, Horizon == 1),
         aes(x = Iteration, y = Mean, ymin = Mean - SE, ymax = Mean + SE, colour = Model, fill = Model)) +
    geom_line(size = 1.5) +
    geom_ribbon(alpha = .5) +
    scale_colour_manual(values = cbbPalette) +
    scale_fill_manual(values = cbbPalette) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(x =  "Training week", y = "MAE", colour = "", fill = "") +
    theme_bw(base_size = 20) +
    theme(legend.position = "top")
}

# Calibration curves -------------------------------------------------------------
# Plot without CI for more readable plot
# NB: calibration curves are obtained by averaging across multiple iterations and might be biased (like uncorrected performance)

res <- readRDS(res_files[mdl_names == "BaseModel"])$Prediction

it <- 0 # median(res$Iteration) # Iteration after which to plot calibration
plot_calibration(subset(res, Iteration > it), cumulative = FALSE, pool = FALSE, CI = NULL, score = score)

# lpd_diff ----------------------------------------------------------------

res1 <- readRDS(res_files[mdl_names == "BaseModel"])$Prediction
res0 <- readRDS(res_files[mdl_names == "Autoregression"])$Prediction

var_id <- c("Patient", "Day", "Iteration", "Horizon", "Severity")

tmp <- merge(res1[, c(var_id, "lpd")],
             res0[, c(var_id, "lpd")],
             by = var_id)
tmp[["lpd_diff"]] <- tmp[["lpd.x"]] - tmp[["lpd.y"]]

lpd_diff <- estimate_performance(tmp, "lpd_diff")$Fit

ggplot(data = subset(lpd_diff, Horizon == 1),
       aes(x = Iteration, y = Mean, ymin = Mean - SE, ymax = Mean + SE)) +
  geom_line(size = 1.5) +
  geom_ribbon(alpha = .5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = expression(Delta*"lpd")) +
  theme_bw(base_size = 20)

# Parameters and Trajectories ---------------------------------------------

par <- readRDS(res_files[mdl_names == "BaseModel"])$Parameters

# Plot parameter estimate (w_S here) as function of training iteration
lapply(pt[1:5],
       function(pid) {
         ggplot(data = subset(par, Variable == "wS" & Patient == pid),
                aes(x = Iteration, y = Mean, ymin = `5%`, ymax = `95%`)) +
           geom_line() +
           geom_ribbon(alpha = .5) +
           coord_cartesian(ylim = c(0, 1)) +
           labs(y =  expression(w["S"])) +
           theme_bw(base_size = 20)
       })
