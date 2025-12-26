# R/loaders/pelland_data_loader.R
# Loader: Pelland et al. Meta-Analysis Data
#
# SOLID Principles Applied:
# - SRP: Single responsibility - load and preprocess Pelland data (actor: data engineer)
# - OCP: Open for extension with different preprocessing strategies

box::use(
  R6[R6Class],
  ../calculators/correlation_imputer[CorrelationImputer, CorrelationImputationResult]
)

#' Pelland Data Loader
#'
#' Loads and preprocesses data from Pelland et al. (2024) meta-regression.
#' Handles all data transformations including correlation imputation and
#' effect size calculation.
#'
#' @description
#' Service class for loading the Pelland dataset from Excel and performing
#' all necessary preprocessing steps to prepare for meta-regression.
#'
#' @export
PellandDataLoader <- R6Class(

  classname = "PellandDataLoader",

  public = list(

    # =========================================================================
    # Constructor
    # =========================================================================

    #' @description Create loader with data path
    #' @param data_path Path to the Pelland Excel data file
    initialize = function(data_path) {
      private$.data_path <- data_path
      private$.correlation_imputer <- CorrelationImputer$new()
    },

    # =========================================================================
    # Public Methods
    # =========================================================================

    #' @description Load and preprocess all data
    #' @return Named list with processed data frames
    load = function() {
      if (!requireNamespace("readxl", quietly = TRUE)) {
        stop("Package 'readxl' is required to load Excel data")
      }
      if (!requireNamespace("metafor", quietly = TRUE)) {
        stop("Package 'metafor' is required for effect size calculation")
      }

      # Load raw data
      raw_data <- private$.load_raw_data()

      # Preprocess data
      preprocessed <- private$.preprocess_data(raw_data)

      # Calculate correlations and impute missing
      with_correlations <- private$.calculate_and_impute_correlations(preprocessed)

      # Calculate effect sizes
      with_effect_sizes <- private$.calculate_effect_sizes(with_correlations$data)

      # Subset by outcome
      strength_data <- private$.subset_by_outcome(with_effect_sizes, "Strength")
      hypertrophy_data <- private$.subset_by_outcome(with_effect_sizes, "Hypertrophy")

      list(
        full_data = with_effect_sizes,
        strength = strength_data,
        hypertrophy = hypertrophy_data,
        correlation_strength = with_correlations$imputation_strength,
        correlation_hypertrophy = with_correlations$imputation_hypertrophy,
        n_observations = nrow(with_effect_sizes)
      )
    },

    #' @description Get summary statistics for loaded data
    #' @param data Loaded data list from load() method
    #' @return Named list with summary statistics
    summarize = function(data) {
      list(
        total_effects = data$n_observations,
        strength_effects = nrow(data$strength),
        hypertrophy_effects = nrow(data$hypertrophy),
        strength_studies = length(unique(data$strength$study)),
        hypertrophy_studies = length(unique(data$hypertrophy$study)),
        imputed_correlation_strength = data$correlation_strength$imputed_correlation,
        imputed_correlation_hypertrophy = data$correlation_hypertrophy$imputed_correlation
      )
    }
  ),

  private = list(

    .data_path = NULL,
    .correlation_imputer = NULL,

    .load_raw_data = function() {
      # Convert tibble to data.frame for consistent subsetting behavior
      as.data.frame(readxl::read_xlsx(private$.data_path))
    },

    .preprocess_data = function(data) {
      # Calculate pre-post SDs from SEs where needed
      data$pre.sd <- ifelse(is.na(data$pre.se), data$pre.sd,
                            data$pre.se * sqrt(data$n))
      data$post.sd <- ifelse(is.na(data$post.se), data$post.sd,
                             data$post.se * sqrt(data$n))

      # Convert p to t (change scores)
      # replmiss replaces NA in first arg with corresponding value from second
      t_from_p <- stats::qt(data$p.value / 2, df = data$n - 1, lower.tail = FALSE)
      data$t.value <- ifelse(is.na(data$t.value), t_from_p, data$t.value)

      # Convert t to SE (change scores)
      se_from_t <- ifelse(
        is.na(data$delta.mean),
        (data$post.mean - data$pre.mean) / data$t.value,
        data$delta.mean / data$t.value
      )
      data$delta.se <- ifelse(is.na(data$delta.se), se_from_t, data$delta.se)

      # Make positive
      data$delta.se <- ifelse(data$delta.se < 0, -data$delta.se, data$delta.se)

      # Convert CI to SE (change scores)
      se_from_ci <- (data$delta.CI.upper - data$delta.CI.lower) / 3.92
      data$delta.se <- ifelse(is.na(data$delta.se), se_from_ci, data$delta.se)

      # Convert SE to SD (change scores)
      sd_from_se <- data$delta.se * sqrt(data$n)
      data$delta.sd <- ifelse(is.na(data$delta.sd), sd_from_se, data$delta.sd)

      data
    },

    .calculate_and_impute_correlations = function(data) {
      imputer <- private$.correlation_imputer

      # Calculate pre-post correlations
      data$ri <- sapply(seq_len(nrow(data)), function(i) {
        imputer$calculate_correlation(
          data$pre.sd[i],
          data$post.sd[i],
          data$delta.sd[i]
        )
      })

      # Impute separately for strength and hypertrophy
      strength_data <- data[data$outcome == "Strength", ]
      hypertrophy_data <- data[data$outcome == "Hypertrophy", ]

      imputation_strength <- imputer$impute_via_meta_analysis(strength_data)
      imputation_hypertrophy <- imputer$impute_via_meta_analysis(hypertrophy_data)

      # Apply imputations
      data[data$outcome == "Strength", ]$ri <- ifelse(
        is.na(data[data$outcome == "Strength", ]$ri),
        imputation_strength$imputed_correlation,
        data[data$outcome == "Strength", ]$ri
      )

      data[data$outcome == "Hypertrophy", ]$ri <- ifelse(
        is.na(data[data$outcome == "Hypertrophy", ]$ri),
        imputation_hypertrophy$imputed_correlation,
        data[data$outcome == "Hypertrophy", ]$ri
      )

      # Estimate change score SD where only pre-post data available
      sd_from_prepost <- sqrt(data$pre.sd^2 + data$post.sd^2 -
                               (2 * data$ri * data$pre.sd * data$post.sd))
      data$delta.sd <- ifelse(is.na(data$delta.sd), sd_from_prepost, data$delta.sd)

      list(
        data = data,
        imputation_strength = imputation_strength,
        imputation_hypertrophy = imputation_hypertrophy
      )
    },

    .calculate_effect_sizes = function(data) {
      # SMCR (Standardized Mean Change using Raw score standardization)
      data <- metafor::escalc(
        measure = "SMCR",
        m1i = post.mean,
        m2i = pre.mean,
        sd1i = pre.sd,
        ni = n,
        ri = ri,
        data = data
      )

      # Add weights
      data$weights <- 1 / sqrt(data$vi)

      # RIR for spline models
      data$spline.rir <- ifelse(data$rir.bucket == "Failure", -1, data$avg.rir)

      data
    },

    .subset_by_outcome = function(data, outcome) {
      # Use which() for robust subsetting that handles NAs
      rows_to_keep <- which(data$outcome == outcome)
      subset_data <- data[rows_to_keep, , drop = FALSE]

      if (length(rows_to_keep) == 0) {
        warning(paste("No rows found for outcome:", outcome))
      }

      # Convert factors
      factor_cols <- c("author.year", "study", "group", "obs", "outcome",
                       "train.status", "set.rep.equated")
      for (col in factor_cols) {
        if (col %in% names(subset_data)) {
          subset_data[[col]] <- factor(subset_data[[col]])
        }
      }

      # Relevel set.rep.equated if present
      if ("set.rep.equated" %in% names(subset_data)) {
        subset_data$set.rep.equated <- factor(
          subset_data$set.rep.equated,
          levels = c("set", "rep", "both")
        )
      }

      subset_data
    }
  )
)
