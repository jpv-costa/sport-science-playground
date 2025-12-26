# scripts/replicate_peerj_velocity.R
# Replication of Paulsen et al. (2025) PeerJ Article
# "Exercise type, training load, velocity loss threshold, and sets affect
#  the relationship between lifting velocity and perceived repetitions in
#  reserve in strength-trained individuals"
#
# Original: https://peerj.com/articles/19797/
# Data: https://doi.org/10.7717/peerj.19797/supp-1

# ==============================================================================
# Setup
# ==============================================================================

library(readxl)
library(tidyverse)
library(lme4)
library(data.table)

cat("=== Paulsen et al. (2025) Velocity-RIR Replication ===\n\n")

# ==============================================================================
# Data Loading
# ==============================================================================

data_path <- "data/external/peerj_2025_velocity/peerj_velocity_rir_data.xlsx"
cat("Loading data from:", data_path, "\n")

# Read all sheets
sheets <- excel_sheets(data_path)
cat("Available sheets:", paste(sheets, collapse = ", "), "\n\n")

# Load main data (usually first sheet)
df <- read_xlsx(data_path, sheet = 1)
cat("Loaded", nrow(df), "rows,", ncol(df), "columns\n")
cat("Columns:", paste(names(df), collapse = ", "), "\n\n")

# ==============================================================================
# Data Exploration
# ==============================================================================

cat("=== Data Summary ===\n")

# Check for key variables - look for mean_velocity specifically
vel_col <- names(df)[grepl("mean.velocity|mean_velocity", names(df), ignore.case = TRUE)][1]
if (!is.na(vel_col)) {
  cat("Velocity column:", vel_col, "\n")
  cat("  Range:", round(range(df[[vel_col]], na.rm = TRUE), 3), "m/s\n")
  cat("  Mean:", round(mean(df[[vel_col]], na.rm = TRUE), 3), "m/s\n")
}

if ("pRIR" %in% names(df) || "RIR" %in% names(df) || "rir" %in% tolower(names(df))) {
  rir_col <- names(df)[grepl("rir|pRIR", names(df), ignore.case = TRUE)][1]
  cat("\nRIR column:", rir_col, "\n")
  cat("  Range:", range(df[[rir_col]], na.rm = TRUE), "\n")
  cat("  Mean:", round(mean(df[[rir_col]], na.rm = TRUE), 2), "\n")
}

# Check for exercise type
if (any(grepl("exercise|type", names(df), ignore.case = TRUE))) {
  ex_col <- names(df)[grepl("exercise|type", names(df), ignore.case = TRUE)][1]
  cat("\nExercise column:", ex_col, "\n")
  cat("  Types:", paste(unique(df[[ex_col]]), collapse = ", "), "\n")
}

# Check for participant ID
if (any(grepl("subject|participant|id", names(df), ignore.case = TRUE))) {
  id_col <- names(df)[grepl("subject|participant|id", names(df), ignore.case = TRUE)][1]
  cat("\nParticipant column:", id_col, "\n")
  cat("  N participants:", length(unique(df[[id_col]])), "\n")
}

# Check for load
if (any(grepl("load|%|1RM", names(df), ignore.case = TRUE))) {
  load_col <- names(df)[grepl("load|%|1RM", names(df), ignore.case = TRUE)][1]
  cat("\nLoad column:", load_col, "\n")
  cat("  Values:", paste(unique(df[[load_col]]), collapse = ", "), "\n")
}

# ==============================================================================
# Data Processing
# ==============================================================================

cat("\n\n=== Data Processing ===\n")

# Standardize column names
df_clean <- df %>%
  rename_with(tolower) %>%
  rename_with(~gsub(" ", "_", .x)) %>%
  rename_with(~gsub("\\.", "_", .x))

cat("Cleaned column names:", paste(names(df_clean), collapse = ", "), "\n\n")

# Print first few rows
cat("First 5 rows:\n")
print(head(df_clean, 5))

# ==============================================================================
# Correlation Analysis (per paper methodology)
# ==============================================================================

cat("\n\n=== Velocity-RIR Correlation Analysis ===\n")

# Identify velocity and RIR columns after cleaning
vel_cols <- names(df_clean)[grepl("mean_velocity|^mv$", names(df_clean), ignore.case = TRUE)]
rir_cols <- names(df_clean)[grepl("perceived.*rir|prir", names(df_clean), ignore.case = TRUE)]

cat("Velocity columns found:", paste(vel_cols, collapse = ", "), "\n")
cat("RIR columns found:", paste(rir_cols, collapse = ", "), "\n")

if (length(vel_cols) > 0 && length(rir_cols) > 0) {
  vel_var <- vel_cols[1]
  rir_var <- rir_cols[1]

  # Overall correlation
  overall_cor <- cor(df_clean[[vel_var]], df_clean[[rir_var]],
                     use = "complete.obs", method = "pearson")
  cat("\nOverall Pearson correlation (velocity vs RIR):", round(overall_cor, 4), "\n")
  cat("R-squared:", round(overall_cor^2, 4), "\n")

  # Per-participant correlations if participant ID exists
  id_cols <- names(df_clean)[grepl("subject|participant|id", names(df_clean), ignore.case = TRUE)]

  if (length(id_cols) > 0) {
    id_var <- id_cols[1]
    cat("\nPer-participant correlations:\n")

    participant_cors <- df_clean %>%
      group_by(across(all_of(id_var))) %>%
      summarise(
        n = n(),
        r = cor(get(vel_var), get(rir_var), use = "complete.obs"),
        r_sq = r^2,
        .groups = "drop"
      )

    cat("  Mean r:", round(mean(participant_cors$r, na.rm = TRUE), 4), "\n")
    cat("  SD r:", round(sd(participant_cors$r, na.rm = TRUE), 4), "\n")
    cat("  Mean r²:", round(mean(participant_cors$r_sq, na.rm = TRUE), 4), "\n")
    cat("  Range r:", round(range(participant_cors$r, na.rm = TRUE), 4), "\n")
  }
}

# ==============================================================================
# Linear Mixed Effects Model (matching paper methodology)
# ==============================================================================

cat("\n\n=== Linear Mixed Effects Model ===\n")

# Check for required columns
required_found <- TRUE
model_vars <- list()

# Try to identify key variables
for (pattern in c("velocity|mv", "rir", "subject|participant|id", "exercise|type", "load|%|1rm", "set")) {
  cols <- names(df_clean)[grepl(pattern, names(df_clean), ignore.case = TRUE)]
  if (length(cols) > 0) {
    model_vars[[pattern]] <- cols[1]
  }
}

cat("Variables identified for modeling:\n")
print(unlist(model_vars))

# Fit basic model if we have the key variables
if (length(model_vars) >= 3) {
  tryCatch({
    vel_var <- model_vars[[1]]
    rir_var <- model_vars[[2]]
    id_var <- model_vars[[3]]

    formula_str <- paste0(rir_var, " ~ ", vel_var, " + (1|", id_var, ")")
    cat("\nFitting model:", formula_str, "\n\n")

    model <- lmer(as.formula(formula_str), data = df_clean)
    cat("Model Summary:\n")
    print(summary(model))

    # Extract fixed effects
    fe <- fixef(model)
    cat("\n\nFixed Effects Interpretation:\n")
    cat("  For every 0.1 m/s decrease in velocity:\n")
    cat("  Change in perceived RIR:", round(-0.1 * fe[2], 2), "\n")

  }, error = function(e) {
    cat("Could not fit mixed model:", e$message, "\n")
  })
} else {
  cat("Insufficient variables identified for mixed model fitting\n")
  cat("Please inspect the data structure manually\n")
}

# ==============================================================================
# Summary Statistics by Exercise Type (if available)
# ==============================================================================

ex_cols <- names(df_clean)[grepl("exercise|type", names(df_clean), ignore.case = TRUE)]
if (length(ex_cols) > 0 && length(vel_cols) > 0 && length(rir_cols) > 0) {
  cat("\n\n=== Summary by Exercise Type ===\n")

  ex_var <- ex_cols[1]
  vel_var <- vel_cols[1]
  rir_var <- rir_cols[1]

  summary_by_ex <- df_clean %>%
    group_by(across(all_of(ex_var))) %>%
    summarise(
      n = n(),
      mean_velocity = mean(get(vel_var), na.rm = TRUE),
      sd_velocity = sd(get(vel_var), na.rm = TRUE),
      mean_rir = mean(get(rir_var), na.rm = TRUE),
      sd_rir = sd(get(rir_var), na.rm = TRUE),
      correlation = cor(get(vel_var), get(rir_var), use = "complete.obs"),
      .groups = "drop"
    )

  print(summary_by_ex)
}

# ==============================================================================
# Save Results
# ==============================================================================

cat("\n\n=== Saving Results ===\n")

results <- list(
  data = df_clean,
  n_observations = nrow(df_clean),
  columns = names(df_clean)
)

saveRDS(results, "data/processed/peerj_velocity_replication_results.rds")
cat("Results saved to data/processed/peerj_velocity_replication_results.rds\n")

cat("\n=== REPLICATION COMPLETE ===\n")
