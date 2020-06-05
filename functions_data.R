# Flares ------------------------------------------------------------------

process_Flares <- function(df, score = c("Bother", "Scratch")) {
  # Process Flares dataset for the BaseModel
  #
  # Args:
  # df: Dataframe of the Flares dataset
  # score: Score to consider (bother or scratch)
  #
  # Returns:
  # Dataframe
  
  score <- match.arg(score)
  stopifnot(is.data.frame(df),
            all(c(score, "Treatment") %in% colnames(df)))
  
  df <- change_colnames(df, score, "Severity")
  df[, c("Date", setdiff(c("Bother", "Scratch"), score))]
  df$Treatment[is.na(df$Treatment)] <- 0 # Imputing missing values by 0
  df <- factor_to_numeric(df, "Treatment")
  
  return(df)
}

# SWET --------------------------------------------------------------------

process1_SWET <- function(df) {
  # Process SWET dataset for the BaseModel
  #
  # Args:
  # df: SWET main dataframe
  #
  # Returns:
  # Dataframe
  
  stopifnot(is.data.frame(df),
            all(c("Bother", "CS") %in% colnames(df)))
  
  df$Severity <- df[, "Bother"]
  df <- change_colnames(df, "CS", "Treatment")
  df[, c("CI", "SU", "Home")] <- NULL
  df$Treatment[is.na(df$Treatment)] <- 0 # Impute missing treatment by 0
  df <- factor_to_numeric(df, "Treatment")
  return(df)
}

process2_SWET <- function(dt) {
  # Process SWET dataset for the extended model
  #
  # Args:
  # dt: SWET main dataframe
  #
  # Returns:
  # List of: 1) dataframe of time series data (trajectories) 2) dataframe of patients' data
  
  stopifnot(is.data.frame(dt),
            all(c("Bother", "CS", "Home", "CI", "SU") %in% colnames(dt)))
  
  ## Time data
  dt <- change_colnames(dt, "Bother", "Severity")
  dt <- factor_to_numeric(dt, c("CS", "Home", "CI", "SU"))
  # Impute missings
  dt$CS[is.na(dt$CS)] <- 0
  dt$CI[is.na(dt$CI)] <- 0
  dt$SU[is.na(dt$SU)] <- 0
  dt$Home[is.na(dt$Home)] <- 1
  
  ## Patient data
  # Demographics data
  demo <- demographics_SWET()
  demo[, c("ethnic", "othdet")] <- NULL
  levels(demo$Sex) <- c(0, 1) # 0 female, 1 male
  # Treatment data
  treat <- treatment_SWET(format = "wide")
  treat$Week <- NULL
  treat$Confidence <- factor(treat$Confidence,
                             levels = c("not at all sure", "not sure", "sure", "very sure"),
                             labels = 1:4)
  treat$Confidence[is.na(treat$Confidence)] <- 1 # only 2 missings (impute with not at all sure)
  # Merge demo and treat
  dp <- merge(demo, treat, by = "Patient", all = FALSE) # Keep data in both demo and treat only
  dp <- factor_to_numeric(dp, c("Sex", "White", "FLG", "Confidence"))
  dp[is.na(dp)] <- 0 # impute missing with 0 (quantity/reference category)
  
  # Number of CS and CI applications in the first 12 weeks
  dp <- merge(dp, aggregate(CS ~ Patient, subset(dt, Day <= 12 * 7), sum), by = "Patient", all = FALSE)
  dp <- merge(dp, aggregate(CI ~ Patient, subset(dt, Day <= 12 * 7), sum), by = "Patient", all = FALSE)
  dp <- change_colnames(dp, c("CS", "CI"), c("N_CS", "N_CI"))

  # Check consistency between dp and dt (patients)
  pt <- intersect(unique(dt$Patient), unique(dp$Patient))
  dt <- dt[dt$Patient %in% pt, ]
  dp <- dp[dp$Patient %in% pt, ]
  
  # Reorder
  dt <- dt[order(dt$Patient), ]
  dp <- dp[order(dp$Patient), ]
  
  return(list(TimeData = dt, PatientData = dp))
}
