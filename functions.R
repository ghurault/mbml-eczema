# Simulate data -----------------------------------------------------------

generate_treatment <- function(p, tmax) {
  # Generate treatment from a markov chain
  #
  # Args:
  # p: vector of length 2 containing the probabilities of using treatment when treatment was not used (p01) or used (p11) the day before respectively
  # To simulate "sticky" (realistic behaviour), p[1] << p[2]
  # tmax: length of time series
  #
  # Returns:
  # Time series of treatment usage
  
  stopifnot(is.numeric(p),
            length(p) == 2,
            min(p) >= 0,
            max(p) <= 1,
            is.numeric(tmax),
            length(tmax) == 1,
            tmax > 0,
            tmax == round(tmax))
  
  Treat <- rep(0, tmax)
  for (i in 2:tmax) {
    if (Treat[i - 1] == 0){
      Treat[i] <- rbinom(1, 1, p[1])
    } else {
      Treat[i] <- rbinom(1, 1, p[2])
    }
  }
  return(Treat)
}

# Parameters --------------------------------------------------------------

get_index <- function(pt, t_max) {
  # Translate unique indices (used in Stan model) into (patient, day) pair
  #
  # Args:
  # pt: Vector of patient ID (same order as the patient parameters in stanfit)
  # t_max: Vector of time-series length (including missing values) for each patient
  #
  # Returns:
  # Dataframe with columns: Patient, Day, Index
  
  stopifnot(length(pt) == length(t_max),
            is.numeric(t_max),
            all(t_max == round(t_max)))
  
  out <- data.frame(Patient = rep(pt, t_max),
                    Day = do.call(c, lapply(t_max, function(x) {1:x})))
  out[["Patient"]] <- as.character(out[["Patient"]])
  out[["Index"]] <- 1:nrow(out)
  return(out)
}

observations_dictionary <- function(pt, data_stan) {
  # Alias of get_index to convert observation parameters' indices into (patient, day) pair
  #
  # Args:
  # pt: vector of patients ID (same order as the patient parameters in stanfit)
  # data_stan: data input to the stan function
  #
  # Returns:
  # Same as get_index
  
  stopifnot(is.list(data_stan),
            "t_max" %in% names(data_stan))
  
  get_index(pt, data_stan$t_max)
}

predictions_dictionary <- function(pt, data_stan) {
  # Alias of get_index to convert prediction parameters' indices into (patient, day) pair 
  #
  # Args:
  # pt: vector of patients ID (same order as the patient parameters in stanfit)
  # data_stan: data input to the stan function
  #
  # Returns:
  # Same as get_index
  
  stopifnot(is.list(data_stan),
            all(c("t_max", "horizon") %in% names(data_stan)))
  
  get_index(pt, data_stan$t_max + data_stan$horizon)
}

extract_parameters <- function(fit, param, param_ind, param_obs, param_pred, pt, data_stan) {
  # Extract parameters' summary
  #
  # Args:
  # fit: stanfit object
  # param: parameters to extract
  # param_ind: individual parameters in param
  # param_obs
  # pt: vector of patients ID (same order as the patient parameters in stanfit)
  # data_stan: data input to the stan function
  #
  # Returns: dataframe containing posterior summary statistics of the parameters 
  
  stopifnot(is.character(c(param, param_ind, param_obs, param_pred)),
            is.list(data_stan),
            length(pt) == data_stan$N_pt)
  # fit checked in HuraultMisc::summary_statistics
  
  par <- HuraultMisc::summary_statistics(fit, param)
  par$Patient <- NA
  par$Day <- NA
  
  pt <- as.character(pt)
  
  ## Patient-dependent parameter
  for (i in intersect(param_ind, param)) {
    idx <- which(par$Variable == i)
    par$Patient[idx] <- pt[par$Index[idx]]
  }
  
  ## Patient and time-dependent parameter (observation parameters)
  dict1 <- observations_dictionary(pt, data_stan)
  for (i in intersect(param_obs, param)) {
    idx <- sort(which(par$Variable == i))
    par[idx, c("Patient", "Day")] <- dict1[, c("Patient", "Day")]
  }
  
  ## Missing score
  par_mis <- intersect("S_mis", param)
  if (length(par_mis) > 0) {
    stopifnot("idx_obs" %in% names(data_stan))
    
    dict1$missing <- TRUE
    dict1$missing[data_stan$idx_obs] <- FALSE
    idx_mis <- which(dict1$missing)
    for (i in par_mis) {
      idx <- which(par$Variable == i)
      id_mis <- idx_mis[par$Index[idx]]
      par[idx, c("Patient", "Day")] <- dict1[id_mis, c("Patient", "Day")]
    }
  }

  ## Prediction parameters
  dict2 <- predictions_dictionary(pt, data_stan)
  for (i in intersect(param_pred, param)) {
    idx <- sort(which(par$Variable == i))
    par[idx, c("Patient", "Day")] <- dict2[, c("Patient", "Day")]
  }
  
  return(par)
}

# Trajectories, predictions and PPC ------------------------------------------------------------

compute_pmf <- function(ps, pred = FALSE) {
  # Compute probability mass function from matrix of posterior samples:
  # - Truncate: discard value outside 0-10
  # - (Truncate trajectories: discard future when past is outside 0-10)
  # - Discretise
  # - Correct rounding at the bounds
  # - Add artificial sample (to avoid probability to be exactly 0)
  #
  # Args:
  # ps: matrix of posterior samples
  # pred: logical indicating whether we are in "prediction mode"
  #
  # Returns:
  # Matrix: observations in rows, pmf in columns
  
  stopifnot(is.matrix(ps),
            is.logical(pred))
  
  # Truncate values outside 0-10
  ps[ps < 0 | ps > 10] <- NA # Truncate value outside 0-10
  
  # If prediction mode, only accept samples when past samples are in 0-10 as well
  if (pred & ncol(ps) > 1) {
    for (i in 1:(ncol(ps) - 1)) {
      ps[is.na(ps[, i]), (i + 1):(ncol(ps))] <- NA
    }
  }
  
  # Discretise/Round
  ps <- round(ps)
  
  # Compute probability table
  pmf <- do.call(rbind,
                 lapply(1:ncol(ps),
                        function(i) {
                          x <- na.omit(ps[, i])
                          n_samp <- length(x)
                          x <- c(x, x[x %in% c(0, 10)]) #  Double the number of samples at the bounds (correct for half support)
                          x <- c(x, 0:10) # Add artificial samples to avoid probability of 0
                          tbl <- table(x) / length(x) # Probability table
                          # If original x doesn't contain enough samples return NA
                          if (n_samp <= 100) {
                            tbl <- tbl * NA
                          }
                          return(tbl)
                        }))
  
  return(pmf)
}

prepare_ppc <- function(fit, df, par, idx_pred) {
  # Prepare a dataset for posterior predictive checks plot
  #
  # Args:
  # fit: Stanfit object
  # df: Dataframe of the data
  # par: Dataframe of parameters' summary
  # idx_pred: Dataframe indexing predictions (cf. predictions_dictionary)
  #
  # Returns:
  # Dataframe with columns: Patient, Day, Severity (data), Fit (Posterior fit), 0 to 10 (Predictions probability)
  
  stopifnot(class(fit) == "stanfit",
            is.data.frame(df),
            all(c("Patient", "Day", "Severity") %in% colnames(df)),
            is.data.frame(par),
            all(c("Patient", "Day", "Mean", "Variable") %in% colnames(par)),
            nrow(subset(par, Variable == "S")) > 0,
            is.data.frame(idx_pred),
            all(c("Patient", "Day", "Index") %in% colnames(idx_pred)))
  
  stopifnot(nrow(par) > 0)
  
  # Compute predictions
  yrep <- extract(fit, pars = "S_pred")[[1]]
  prob <- compute_pmf(yrep, pred = FALSE)
  pred <- cbind(idx_pred, prob)
  
  
  
  # Merge df and fit (in par)
  out <- merge(df[, c("Patient", "Day", "Severity")],
               par[par[["Variable"]] == "S", c("Patient", "Day", "Mean")],
               by = c("Patient", "Day"), all = TRUE)
  # Merge results with predictions
  out <- merge(out, pred, by = c("Patient", "Day"), all.x = TRUE, all.y = FALSE)
  # Minor processing
  out <- out[order(out[["Patient"]], out[["Day"]]), ]
  out <- change_colnames(out, "Mean", "Fit")
  
  return(out)
}

plot_ppc <- function(ppc, patientID) {
  # Plot Posterior predictive distribution (density/mass plot)
  #
  # Args:
  # ppc: Dataframe, output from prepare_ppc
  # patientID: Patient ID
  #
  # Returns:
  # Ggplot
  
  library(ggplot2)
  palette <- c("#FFFFFF", RColorBrewer::brewer.pal(n = 6, "Blues")) # white-blue
  gamma <- 1.5
  col_correction <- scales::rescale(seq(0, 1, .1)^gamma)
  
  stopifnot(is.data.frame(ppc),
            all(c("Patient", "Day", "Severity", "Fit", 0:10) %in% colnames(ppc)),
            patientID %in% unique(ppc[["Patient"]]))
  
  tmp <- subset(ppc, Patient == as.character(patientID))
  
  ## Process trajectories (observed and fit)
  traj <- tmp[, c("Patient", "Day", "Severity", "Fit")]
  
  # If we have a pattern of points Observed Missing Observed (O-M-O), the O-M line will be marked as Observed and the M-O line will be marked as Missing
  # instead, we want both line to be marked as Missing
  traj[["Observed"]] <- !is.na(traj[["Severity"]])
  traj$ColourLine <- as.logical(traj$Observed)
  traj$ColourLine[which(!as.logical(traj$ColourLine)) - 1] <- FALSE
  traj$ColourLine <- factor(traj$ColourLine, levels = c(TRUE, FALSE))
  
  # Only show fit when missing
  traj$Fit[traj$Observed] <- round(traj$Fit[traj$Observed])
  # Offset by 1 to align with factors
  traj$Fit <- traj$Fit + 1
  
  ## Process predictions
  pred <- tmp[, c("Patient", "Day", 0:10)]
  pred<- reshape2::melt(pred,
                        id.vars = c("Patient", "Day"),
                        variable.name = "S",
                        value.name = "Probability")
  
  ## Plot
  # Heatmap
  p <- ggplot() +
    geom_tile(data = pred, aes(x = Day, y = S, fill = Probability)) +
    scale_fill_gradientn(colours = palette, limits = c(0, 1), breaks = c(0, .5, 1), values = col_correction) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA))
  # Overlay trajectory
  p <- p +
    geom_line(data = traj, aes(x = Day, y = Fit, colour = ColourLine, group = 1), lwd = 1.5) + # Trajectory (observed and missing)
    geom_point(data = traj[filter(as.numeric(as.logical(traj$Observed)), c(-1, 1, -1), method = "convolution", sides = 2) == 1, ],
               aes(x = Day, y = Fit), colour = "black") + # show the observed points in the case of a pattern M-O-M
    scale_color_manual("", labels = c("Observed", "Missing"), values = c("black","grey"))
  # Formatting
  p <-  p +
    labs(y = "Severity score (S)") +
    theme_classic(base_size = 20) +
    theme(legend.position = "top")
  
  return(p)
}

# Validation process ------------------------------------------------------------------

prediction_matrices <- function(res) {
  # Forecast and outcome matrices
  #
  # Args:
  # res: Prediction dataframe
  #
  # Returns:
  # List containing Forecast and outcome matrices
  
  lbl <- paste0("P(S=", 0:10, ")")
  
  stopifnot(is.data.frame(res),
            all(c(lbl, "Severity") %in% colnames(res)))
  
  f <- as.matrix(res[, lbl])
  o <- f * 0
  for (i in 1:nrow(o)) {
    o[i, res$Severity[i] + 1] <- 1
  }
  list(Forecast = f, Outcome = o)
}

plot_calibration <- function(res, cumulative = FALSE, pool = FALSE, CI = NULL, score = "Severity", ...) {
  # Calibration plot obtained by LOWESS smoothing
  #
  # Args:
  # res: dataframe of predictions/actual (typically the online learning output)
  # cumulative: whether to compute the calibration of the cumulative forecasts
  # pool: whether to show the calibration of the scores pooled together or stratified (0, 1, 2, ..., 10)
  # CI: confidence level (e.g. 0.95). If NULL, confidence bounds are not shown
  # score: name of the score
  # ...: arguments of loess
  #
  # Returns:
  # Ggplot of calibration curves
  
  library(ggplot2)
  
  stopifnot(is.logical(cumulative),
            is.logical(pool),
            is.null(CI) || (is.numeric(CI) && CI > 0 && CI < 1),
            is.character(score))

  l <- prediction_matrices(na.omit(res))
  f <- l$Forecast
  o <- l$Outcome
  
  if (cumulative) {
    f <- t(apply(f, 1, cumsum))
    o <- t(apply(o, 1, cumsum))
  }
  
  if (pool) {
    tmp <- HuraultMisc::compute_calibration(reshape2::melt(f, id.vars = c())$value,
                                            reshape2::melt(o, id.vars = c())$value, CI = CI, ...)
    tmp$Label <- "Pooled"
  } else {
    if (!cumulative) {x <- 1:ncol(f)} else {x <- 1:(ncol(f) - 1)}
    tmp <- do.call(rbind,
                   lapply(x,
                          function(i) {
                            cal <- HuraultMisc::compute_calibration(f[, i],
                                                                    o[, i], CI = CI, ...)
                            cal$Label <- i - 1
                            return(cal)
                          }))
    if (!cumulative) {
      tmp$Label <- factor(tmp$Label, levels = x - 1, labels = paste(score, "=", x - 1))
    } else {
      tmp$Label <- factor(tmp$Label, levels = x - 1, labels = paste(score, "<=", x - 1))
    }
  }
  
  p <- ggplot(data = tmp,
              aes(x = Forecast, y = Frequency, colour = Label, fill = Label)) +
    geom_line(size  = 1.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") + # Reference
    labs(x = "Predicted Probability", y = "Observed Frequency", colour = "", fill = "") +
    scale_x_continuous(breaks = seq(0, 1, .1), limits = c(0, 1)) +
    scale_y_continuous(breaks = seq(0, 1, .1), limits = c(0, 1)) +
    annotate("text", x = .2, y = .8, label = "Underconfident", size = 6) +
    annotate("text", x = .8, y = .2, label = "Overconfident", size = 6) +
    theme_bw(base_size = 20) + theme(panel.grid.minor = element_blank())
  
  if (!pool) {
    palette <- c('#A50026', '#D73027', '#F46D43', '#FDAE61', '#FEE090',
                 '#FFFFBF', '#E0F3F8', '#ABD9E9', '#74ADD1', '#4575B4', '#313695')
    p <- p +
      scale_colour_manual(values = palette) +
      scale_fill_manual(values = palette)
  }
  
  if (!is.null(CI)) {
    p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.5)
  }
  
  return(p)
}

compute_RPS <- function(forecast, outcome) {
  # Compute RPS (for a single forecast)
  #
  # Args:
  # forecast: Vector of length N
  # outcome: Index of true outcome (between 1 and N)
  #
  # Returns:
  # RPS
  
  stopifnot(is.vector(forecast),
            length(forecast) > 1)
  
  if (any(is.na(c(forecast, outcome)))) {
    return(NA)
  } else {
    stopifnot(all(forecast >=0 & forecast <= 1),
              round(sum(forecast), 2) == 1,
              outcome %in% 1:length(forecast))
    dummy_outcome <- 0 * forecast
    dummy_outcome[outcome] <- 1
    RPS <- sum((cumsum(forecast) - cumsum(dummy_outcome))^2) / (length(forecast) - 1)
    return(RPS)
  }
}

estimate_performance <- function(res, metric) {
  # Estimate the evolution of a performance metric while controlling for patient-dependence, missing observations and prediction horizon...
  # ... using a generative additive model.
  #
  # Args:
  # res: dataframe of predictions/actual, output of the forward chaining
  # metric: metric to estimate the performance for (RPS or lpd)
  #
  # Returns:
  # List containing the model and the dataframe of the fit
  
  library(gamm4)
  
  stopifnot(is.data.frame(res),
            all(c(metric, "Horizon", "Iteration", "Patient") %in% colnames(res)))
  
  # Prepare dataset
  tmp <- res
  tmp$HorizonMinus1 <- tmp$Horizon - 1 # Prediction horizon (0 corresponds to the first step ahead so that the intercept is the prediction one step ahead)
  tmp$Iteration0 <- as.numeric(tmp$Iteration == 0)
  tmp$IterationGreaterThan0 <- 1 - tmp$Iteration0
  tmp$Iteration1 <- tmp$Iteration * tmp$IterationGreaterThan0 # Interaction
  
  # Prediction dataframe
  pred <- unique(tmp[, c( "Iteration", "Horizon", "HorizonMinus1", "Iteration0", "IterationGreaterThan0", "Iteration1")])
  pred <- na.omit(pred) # possible to have missing in tmp (Horizon NA at Iteration 0 because Severity at Day 1 missing but not Treatment)
  pred$Patient <- 0 # new patient
  
  # GAM
  # NB: average of smoothing spline 0 so intercept represents an "average performance"
  f <- formula(paste(metric, " ~ Iteration0 + IterationGreaterThan0:HorizonMinus1 + s(Iteration1, bs = 'cr')"))
  mdl <- gamm4(f, random = ~ (1 | Patient), data = tmp)
  # print(mdl); gam.check(mdl$gam)
  
  mdl_fit <- predict(mdl$gam, newdata = pred, se.fit = TRUE)
  pred$Mean <- mdl_fit$fit
  pred$SE <- mdl_fit$se.fit
  
  pred[, c("HorizonMinus1", "Iteration0", "Iteration1", "IterationGreaterThan0", "Patient")] <- NULL
  pred <- pred[order(pred$Iteration, pred$Horizon), ]
  
  return(list(Model = mdl, Fit = pred))
}
