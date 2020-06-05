# Figures in the manuscript:
# - Model: not in R
# - PPC Flares Bother: in fitting.R
# - Performance: here
# - PPC SWET ext: in fitting_ext.R
# - Covariates SWET: in fitting_ext.R

# Figures in the supplementary:
# - Example data trajectories: here
# - Forward chaining: here
# - Missing Flares: here
# - Missing SWET: here
# - Factor Graph: not in R
# - Patient-dependent parameters Flares: in fitting.R
# - Patient-dependent parameters SWET: in fitting.R
# - Calibration curves: here
# - Results Flares Scratch: here
# - Simpson's paradox: here

# Initialisation ----------------------------------------------------------

rm(list = ls()) # Clear Workspace

library(TanakaData) # Contains data and data processing functions
library(HuraultMisc) # Functions shared across projects
library(ggplot2)
library(cowplot)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
source("functions.R") # Additional functions

# Figure performance ------------------------------------------------------

RPS_Flares <- readRDS("Plots/RPS_Flares_Bother.rds")
lpd_Flares <- readRDS("Plots/lpd_Flares_Bother.rds")
RPS_SWET <- readRDS("Plots/RPS_SWET_Bother.rds")
lpd_SWET <- readRDS("Plots/lpd_SWET_Bother.rds")

plot_grid(get_legend(RPS_Flares + theme(legend.position = "top")),
          plot_grid(plotlist = lapply(list(RPS_Flares + labs(title = "Flares"),
                                           RPS_SWET + labs(title = "SWET"),
                                           lpd_Flares + labs(title = "Flares"),
                                           lpd_SWET + labs(title = "SWET")),
                                      function(x) {x + theme(legend.position = "none")}),
                    ncol = 2, labels = "auto"),
          ncol = 1, rel_heights = c(.1, .9)) # NULL for second dataset
# ggsave("Plots/FigPerformance.jpg", width = 20, height = 15, units = "cm", dpi = 300, scale = 1.3)

if (FALSE) {
  AUROC_Flares <- readRDS("Plots/AUROC_Flares_Bother.rds")
  AUROC_SWET <- readRDS("Plots/AUROC_SWET_Bother.rds")
  
  plot_grid(get_legend(AUROC_Flares + theme(legend.position = "top")),
            plot_grid(plotlist = lapply(list(AUROC_Flares + labs(title = "Flares"),
                                             AUROC_SWET + labs(title = "SWET")),
                                        function(x) {
                                          x + theme(legend.position = "none")
                                        }),
                      ncol = 2, labels = "auto"),
            ncol = 1, rel_heights = c(.1, .9))
  # ggsave("Plots/AUROC_perf.jpg", width = 20, height = 12, units = "cm", dpi = 300, scale = 1.3)
}

# Calibration curves ------------------------------------------------------

datasets <- c("Flares", "SWET")
res_files <- file.path("Results", paste0("val_BaseModel_", datasets, "_Bother.rds"))

pl <- lapply(1:length(datasets),
             function(i) {
               plot_calibration(subset(readRDS(res_files[i])$Prediction, Iteration > 0),
                                cumulative = FALSE,
                                pool = FALSE,
                                CI = NULL,
                                score = "Bother") +
                 labs(title = datasets[i])
             })
plot_grid(plot_grid(plotlist = lapply(pl, function(x) {x + theme(legend.position = "none")}),
                    nrow = 1, labels = "auto"),
          get_legend(pl[[1]] + theme(legend.position = "right")),
          nrow = 1, rel_widths = c(.85, .15))
# ggsave("Plots/SuppCalibration.jpg", width = 20, height = 10, units = "cm", dpi = 300, scale = 1.5)

# Results Scratch ----------------------------------------------------------

RPS_Flares <- readRDS("Plots/RPS_Flares_Scratch.rds")
lpd_Flares <- readRDS("Plots/lpd_Flares_Scratch.rds")
res <- readRDS("Results/val_BaseModel_Flares_Scratch.rds")$Prediction

cal_plot <- plot_calibration(subset(res, Iteration > 0),
                             cumulative = FALSE,
                             pool = FALSE,
                             CI = NULL,
                             score = "Scratch")

perf_plot <- plot_grid(get_legend(RPS_Flares + theme(legend.position = "top")),
                       plot_grid(plotlist = lapply(list(RPS_Flares, lpd_Flares),
                                                   function(x) {
                                                     x + theme(legend.position = "none")
                                                   }),
                                 nrow = 1, labels = c("a", "b")),
                       ncol = 1, rel_heights = c(.1, .9))

plot_grid(perf_plot,
          plot_grid(NULL, cal_plot, NULL,
                    nrow = 1, rel_widths = c(.1, .8, .1), labels = c("", "c", "")),
          ncol = 1)
# ggsave("Plots/SuppScratchPerformance.jpg", width = 20, height = 20, units = "cm", dpi = 300, scale = 1.2)

# Forward chaining ------------------------------------------------------------------

HuraultMisc::illustrate_forward_chaining()
# ggsave("Plots/SuppForwardChaining.tiff", width = 25, height = 15, units = "cm", dpi = 300)

# Missing score ----------------------------------------------------

df <- load_Flares()
score <- "Bother"
df$Observed <- factor(!is.na(df[[score]]), levels = c(TRUE, FALSE))

ggplot(data = df, aes(y = Patient, x = Date, fill = Observed)) +
  geom_tile() +
  scale_fill_manual(values = cbbPalette, labels = c("Observed", "Missing")) +
  labs(fill = score) +
  theme_classic(base_size = 30) +
  theme(axis.text.y = element_blank())
# ggsave("Plots/SuppMissingFlares.tiff", width = 40, height = 40, units = "cm", dpi = 300)

TanakaData::plot_missing_SWET(SWET, date = TRUE) +
  theme_classic(base_size = 40) +
  theme(axis.text.y = element_blank())
# ggsave("Plots/SuppMissingSWET.jpg", width = 20, height = 40, units = "cm", dpi = 300, limitsize = FALSE, scale = 2.5)

# Simpson Paradox -------------------------------------------

res <- readRDS("Results/val_BaseModel_Flares_Bother.rds")$Prediction

it <- 29 # Iteration cut-off

dp <- aggregate(Iteration ~ Patient, res, max)
dp <- change_colnames(dp, "Iteration", "LastIteration")

d1 <- aggregate(RPS ~ Iteration, res, mean)
d1$Label <- "All"

d2 <- aggregate(RPS ~ Iteration, subset(res, Patient %in% dp$Patient[dp$LastIteration < it]), mean)
d2$Label <- paste0("Leave_Before_", it)

d3 <- aggregate(RPS ~ Iteration, subset(res, Patient %in% dp$Patient[dp$LastIteration >= it]), mean)
d3$Label <- paste0("Leave_After_", it)

ggplot(data = do.call(rbind, list(d1, d2, d3)),
       aes(x = Iteration, y = RPS, colour = Label)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +
  geom_vline(xintercept = it - 1, linetype = "dashed") +
  # scale_y_continuous(limits = c(0, 0.2)) +
  scale_colour_manual(values = c("#000000", "#E69F00", "#56B4E9")) +
  labs(colour = "") +
  theme_bw(base_size = 15) +
  theme(legend.position = "top")
# ggsave("Plots/SuppSimpsonParadox.jpg", width = 13, height = 8, units = "cm", dpi = 300, scale = 1.5)

# Plot trajectory data ----------------------------------------------------

# patient 6006 (same as Fig patient responsiveness) to illustrate the complexity of the data
# patient 2064 to illustrate missing values

pl <- lapply(c("6006", "2064"),
             function(patientID) {
               df <- subset(SWET, Patient == patientID)
               df <- factor_to_numeric(df, c("CS", "CI", "SU", "Home"))
               
               var_names <- c("Bother", "SU", "CS", "CI", "Home")
               var_lbl <- c("Bother", "Step-up", "Corticosteroids", "Calcineurin \nInhibitors", "Sleeping \nat home")
               
               pl <- lapply(1:length(var_names),
                            function(i) {
                              p <- ggplot(data = df, aes_string(x = "Day", y = var_names[i])) +
                                geom_path(size = 1.5) +
                                geom_point(size = 0.8) + # show observed values surrounded by missings
                                labs(y = var_lbl[i]) +
                                theme_bw(base_size = 15) +
                                theme(panel.grid.minor.y = element_blank(),
                                      axis.title.y = element_text(angle = 0, vjust = 0.5))
                              if (var_names[i] == "Bother") {
                                p <- p + scale_y_continuous(limits = c(0, 10), breaks = 0:10)
                              } else {
                                p <- p +
                                  scale_y_continuous(limits = c(0, 1), breaks = 0:1, labels = c("No", "Yes")) +
                                  theme(axis.title.x = element_blank()) # Remove "Day" to save space
                              }
                              return(p)
                            })
               
               plot_grid(plotlist = pl,
                         rel_heights = c(2, 1, 1, 1, 1),
                         ncol = 1, align = "v")
             })

plot_grid(plotlist = pl, ncol = 2, labels = "auto")
# ggsave("Plots/SuppDataTrajectories.jpg", width = 40, height = 20, units = "cm", dpi = 300)
