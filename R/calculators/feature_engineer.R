# R/calculators/feature_engineer.R
# Service: Comprehensive Feature Engineering for Anomaly Detection
#
# =============================================================================
# ARCHITECTURE (follows SOLID/CUPID principles):
#
# FeatureEngineer (Facade) - orchestrates feature creation
#   ├── GlobalStatsCalculator - computes population-level statistics
#   ├── SetFeatureCalculator - Group 2: set-level features
#   ├── ParticipantFeatureCalculator - Groups 3, 8, 9: participant features
#   ├── CrossLevelFeatureCalculator - Groups 4, 5, 6: cross-level features
#   └── ObservationFeatureCalculator - Groups 1, 10, 11, 12: observation features
#
# Each calculator has single responsibility (~30 lines per method)
# =============================================================================

box::use(
  R6[R6Class],
  stats[sd, lm, coef, median, var, quantile, IQR],
  dplyr[
    group_by, summarize, mutate, ungroup, n, left_join, .data,
    row_number, first, last, lag, arrange, n_distinct, select, any_of
  ]
)

# =============================================================================
# VALUE OBJECTS
# =============================================================================

#' Feature Engineering Result
#'
#' Immutable container for all engineered features across aggregation levels.
#'
#' @export
FeatureEngineeringResult <- R6Class(

  classname = "FeatureEngineeringResult",
  cloneable = FALSE,

 public = list(
    #' @field observation_features Per-observation features
    observation_features = NULL,
    #' @field set_features Per-set aggregates
    set_features = NULL,
    #' @field participant_features Per-participant aggregates
    participant_features = NULL,
    #' @field participant_set_features Participant × set features
    participant_set_features = NULL,
    #' @field participant_load_features Participant × load features
    participant_load_features = NULL,
    #' @field participant_rir_features Participant × RIR features
    participant_rir_features = NULL,
    #' @field feature_names Feature names by aggregation level
    feature_names = NULL,
    #' @field global_stats Population statistics for normalization
    global_stats = NULL,

    #' @description Create result container
    initialize = function(observation_features, set_features, participant_features,
                          participant_set_features, participant_load_features,
                          participant_rir_features, feature_names, global_stats) {
      self$observation_features <- observation_features
      self$set_features <- set_features
      self$participant_features <- participant_features
      self$participant_set_features <- participant_set_features
      self$participant_load_features <- participant_load_features
      self$participant_rir_features <- participant_rir_features
      self$feature_names <- feature_names
      self$global_stats <- global_stats
    },

    #' @description Summarize feature counts
    summarize = function() {
      list(
        n_observations = nrow(self$observation_features),
        n_sets = nrow(self$set_features),
        n_participants = nrow(self$participant_features),
        total_features = sum(sapply(self$feature_names, length))
      )
    },

    #' @description Extract observation feature matrix for anomaly detection
    get_observation_matrix = function() {
      cols <- intersect(self$feature_names$observation, names(self$observation_features))
      self$observation_features[, cols, drop = FALSE]
    },

    #' @description Extract participant feature matrix for anomaly detection
    get_participant_matrix = function() {
      cols <- intersect(self$feature_names$participant, names(self$participant_features))
      self$participant_features[, cols, drop = FALSE]
    }
  )
)

# =============================================================================
# HELPER: Statistical Utilities (stateless functions)
# =============================================================================

StatUtils <- R6Class(
  classname = "StatUtils",
  cloneable = FALSE,

  public = list(
    #' @description Safe coefficient of variation
    cv = function(x) {
      x <- x[!is.na(x)]
      if (length(x) < 2) return(NA_real_)
      m <- mean(x)
      if (abs(m) < 1e-10) return(NA_real_)
      sd(x) / m
    },

    #' @description Safe linear regression slope
    slope = function(y, x) {
      valid <- !is.na(y) & !is.na(x)
      if (sum(valid) < 2 || length(unique(x[valid])) < 2) return(NA_real_)
      fit <- tryCatch(lm(y[valid] ~ x[valid]), error = function(e) NULL)
      if (is.null(fit)) return(NA_real_)
      coef(fit)[2]
    },

    #' @description Safe linear regression intercept
    intercept = function(y, x) {
      valid <- !is.na(y) & !is.na(x)
      if (sum(valid) < 2) return(NA_real_)
      if (length(unique(x[valid])) < 2) return(mean(y[valid], na.rm = TRUE))
      fit <- tryCatch(lm(y[valid] ~ x[valid]), error = function(e) NULL)
      if (is.null(fit)) return(mean(y[valid], na.rm = TRUE))
      coef(fit)[1]
    },

    #' @description Safe R-squared
    r_squared = function(y, x) {
      valid <- !is.na(y) & !is.na(x)
      if (sum(valid) < 3 || length(unique(x[valid])) < 2) return(NA_real_)
      fit <- tryCatch(lm(y[valid] ~ x[valid]), error = function(e) NULL)
      if (is.null(fit)) return(NA_real_)
      tryCatch(summary(fit)$r.squared, error = function(e) NA_real_)
    },

    #' @description Safe skewness (Fisher's definition)
    skewness = function(x) {
      x <- x[!is.na(x)]
      if (length(x) < 3) return(NA_real_)
      m <- mean(x)
      s <- sd(x)
      if (s < 1e-10) return(NA_real_)
      sum((x - m)^3) / (length(x) * s^3)
    },

    #' @description Safe excess kurtosis
    kurtosis = function(x) {
      x <- x[!is.na(x)]
      if (length(x) < 4) return(NA_real_)
      m <- mean(x)
      s <- sd(x)
      if (s < 1e-10) return(NA_real_)
      sum((x - m)^4) / (length(x) * s^4) - 3
    },

    #' @description Parse load percentage string to numeric
    parse_load = function(load_pct) {
      as.numeric(gsub("%", "", load_pct))
    }
  )
)

# =============================================================================
# CALCULATOR: Global Statistics
# =============================================================================

GlobalStatsCalculator <- R6Class(
  classname = "GlobalStatsCalculator",
  cloneable = FALSE,

  public = list(
    #' @description Calculate population-level statistics
    #' @param data Raw velocity data
    #' @return Named list of global statistics
    calculate = function(data) {
      list(
        mean_velocity = mean(data$mean_velocity, na.rm = TRUE),
        sd_velocity = max(sd(data$mean_velocity, na.rm = TRUE), 1e-10),
        median_velocity = median(data$mean_velocity, na.rm = TRUE),
        mean_rir = mean(data$rir, na.rm = TRUE),
        sd_rir = max(sd(data$rir, na.rm = TRUE), 1e-10),
        mean_cv = NA_real_,
        mean_rir_slope = NA_real_,
        mean_load_slope = NA_real_
      )
    },

    #' @description Update global stats with participant-derived values
    update_from_participants = function(global_stats, participant_features) {
      global_stats$mean_cv <- mean(participant_features$p_velocity_cv, na.rm = TRUE)
      global_stats$mean_rir_slope <- mean(participant_features$p_rir_slope, na.rm = TRUE)
      global_stats$mean_load_slope <- mean(participant_features$p_load_sensitivity, na.rm = TRUE)
      global_stats
    }
  )
)

# =============================================================================
# CALCULATOR: Set-Level Features (Group 2)
# =============================================================================

SetFeatureCalculator <- R6Class(
  classname = "SetFeatureCalculator",
  cloneable = FALSE,

  public = list(
    initialize = function() {
      private$.stats <- StatUtils$new()
    },

    #' @description Calculate set-level aggregate features
    calculate = function(data) {
      set_features <- private$.calculate_base_features(data)
      set_features <- private$.add_sequence_numbers(set_features)
      private$.add_cross_set_ratios(set_features)
    }
  ),

  private = list(
    .stats = NULL,

    .calculate_base_features = function(data) {
      stats <- private$.stats
      data |>
        group_by(.data$id, .data$set_id) |>
        summarize(
          set_mean_velocity = mean(.data$mean_velocity, na.rm = TRUE),
          set_median_velocity = median(.data$mean_velocity, na.rm = TRUE),
          set_velocity_sd = sd(.data$mean_velocity, na.rm = TRUE),
          set_velocity_cv = stats$cv(.data$mean_velocity),
          set_velocity_range = max(.data$mean_velocity, na.rm = TRUE) -
                               min(.data$mean_velocity, na.rm = TRUE),
          set_velocity_iqr = IQR(.data$mean_velocity, na.rm = TRUE),
          set_max_velocity = max(.data$mean_velocity, na.rm = TRUE),
          set_min_velocity = min(.data$mean_velocity, na.rm = TRUE),
          set_first_rep_velocity = first(.data$mean_velocity),
          set_last_rep_velocity = last(.data$mean_velocity),
          set_velocity_drop = first(.data$mean_velocity) - last(.data$mean_velocity),
          set_velocity_drop_pct = private$.velocity_drop_pct(
            first(.data$mean_velocity), last(.data$mean_velocity)
          ),
          set_decay_slope = stats$slope(.data$mean_velocity, .data$rep_number),
          set_decay_r2 = stats$r_squared(.data$mean_velocity, .data$rep_number),
          set_reps_count = n(),
          set_load = first(.data$load_percentage),
          set_rir_start = first(.data$rir),
          set_rir_end = last(.data$rir),
          .groups = "drop"
        )
    },

    .add_sequence_numbers = function(set_features) {
      set_features |>
        group_by(.data$id) |>
        mutate(set_sequence_in_session = row_number()) |>
        ungroup()
    },

    .add_cross_set_ratios = function(set_features) {
      set_features |>
        group_by(.data$id) |>
        arrange(.data$set_sequence_in_session) |>
        mutate(
          consecutive_set_velocity_ratio = private$.safe_ratio(
            .data$set_mean_velocity, lag(.data$set_mean_velocity)
          ),
          consecutive_set_velocity_diff = .data$set_mean_velocity - lag(.data$set_mean_velocity),
          set_to_session_first_ratio = .data$set_mean_velocity / first(.data$set_mean_velocity),
          set_to_session_best_ratio = .data$set_mean_velocity /
                                      max(.data$set_mean_velocity, na.rm = TRUE),
          cumulative_volume_at_set = cumsum(lag(.data$set_reps_count, default = 0)),
          inter_set_recovery = .data$consecutive_set_velocity_diff > 0
        ) |>
        ungroup()
    },

    .velocity_drop_pct = function(first_vel, last_vel) {
      if (is.na(first_vel) || first_vel <= 0) return(NA_real_)
      (first_vel - last_vel) / first_vel * 100
    },

    .safe_ratio = function(numerator, denominator) {
      ifelse(is.na(denominator) | denominator < 1e-10, NA_real_, numerator / denominator)
    }
  )
)

# =============================================================================
# CALCULATOR: Participant-Level Features (Groups 3, 8, 9)
# =============================================================================

ParticipantFeatureCalculator <- R6Class(
  classname = "ParticipantFeatureCalculator",
  cloneable = FALSE,

  public = list(
    initialize = function() {
      private$.stats <- StatUtils$new()
    },

    #' @description Calculate participant-level features
    calculate = function(data, global_stats) {
      p_features <- private$.calculate_base_features(data)
      p_features <- private$.add_global_context(p_features, global_stats)
      private$.add_load_curve_features(data, p_features, global_stats)
    }
  ),

  private = list(
    .stats = NULL,

    .calculate_base_features = function(data) {
      stats <- private$.stats
      data |>
        group_by(.data$id) |>
        summarize(
          # Central tendency
          p_mean_velocity = mean(.data$mean_velocity, na.rm = TRUE),
          p_median_velocity = median(.data$mean_velocity, na.rm = TRUE),
          p_velocity_sd = sd(.data$mean_velocity, na.rm = TRUE),
          p_velocity_cv = stats$cv(.data$mean_velocity),
          p_velocity_iqr = IQR(.data$mean_velocity, na.rm = TRUE),
          p_velocity_var = var(.data$mean_velocity, na.rm = TRUE),

          # Range/extremes
          p_velocity_range = max(.data$mean_velocity, na.rm = TRUE) -
                             min(.data$mean_velocity, na.rm = TRUE),
          p_max_velocity = max(.data$mean_velocity, na.rm = TRUE),
          p_min_velocity = min(.data$mean_velocity, na.rm = TRUE),

          # Percentiles
          p_percentile_10 = quantile(.data$mean_velocity, 0.10, na.rm = TRUE),
          p_percentile_90 = quantile(.data$mean_velocity, 0.90, na.rm = TRUE),

          # Distribution shape
          p_velocity_skewness = stats$skewness(.data$mean_velocity),
          p_velocity_kurtosis = stats$kurtosis(.data$mean_velocity),

          # RIR patterns
          p_rir_slope = stats$slope(.data$mean_velocity, .data$rir),
          p_rir_intercept = stats$intercept(.data$mean_velocity, .data$rir),
          p_rir_r2 = stats$r_squared(.data$mean_velocity, .data$rir),
          p_rir_range = max(.data$rir, na.rm = TRUE) - min(.data$rir, na.rm = TRUE),
          p_mean_rir = mean(.data$rir, na.rm = TRUE),

          # Set patterns
          p_n_sets = n_distinct(.data$set_id),
          p_n_observations = n(),
          p_mean_reps_per_set = mean(unique(.data$reps_to_failure), na.rm = TRUE),
          p_reps_cv = stats$cv(unique(.data$reps_to_failure)),
          p_max_reps_in_set = max(.data$reps_to_failure, na.rm = TRUE),
          p_min_reps_in_set = min(.data$reps_to_failure, na.rm = TRUE),
          p_reps_range = max(.data$reps_to_failure, na.rm = TRUE) -
                         min(.data$reps_to_failure, na.rm = TRUE),

          # Load patterns
          p_load_sensitivity = private$.load_sensitivity(.data$mean_velocity, .data$load_percentage),
          p_load_range = private$.load_range(.data$load_percentage),
          p_mean_load = mean(stats$parse_load(.data$load_percentage), na.rm = TRUE),
          p_n_loads = n_distinct(.data$load_percentage),

          .groups = "drop"
        )
    },

    .add_global_context = function(p_features, global_stats) {
      p_features |>
        mutate(
          p_rir_slope_deviation = .data$p_rir_slope - global_stats$mean_rir_slope,
          p_velocity_ratio_to_population = .data$p_mean_velocity / global_stats$mean_velocity,
          p_velocity_z_global = (.data$p_mean_velocity - global_stats$mean_velocity) /
                                global_stats$sd_velocity,
          p_cv_ratio_to_population = private$.safe_ratio(
            .data$p_velocity_cv, global_stats$mean_cv
          ),
          p_rir_sensitivity_ratio = private$.safe_ratio(
            .data$p_rir_slope, global_stats$mean_rir_slope
          ),
          p_load_sensitivity_ratio = private$.safe_ratio(
            .data$p_load_sensitivity, global_stats$mean_load_slope
          ),
          p_velocity_percentile = rank(.data$p_mean_velocity) / n() * 100,
          p_consistency_percentile = rank(-.data$p_velocity_cv) / n() * 100
        )
    },

    .add_load_curve_features = function(data, p_features, global_stats) {
      stats <- private$.stats
      load_curves <- data |>
        mutate(load_numeric = stats$parse_load(.data$load_percentage)) |>
        group_by(.data$id) |>
        summarize(
          p_load_curve_slope = stats$slope(.data$mean_velocity, .data$load_numeric),
          p_load_curve_intercept = stats$intercept(.data$mean_velocity, .data$load_numeric),
          p_load_curve_r2 = stats$r_squared(.data$mean_velocity, .data$load_numeric),
          .groups = "drop"
        )

      global_slope <- private$.stats$slope(
        data$mean_velocity,
        private$.stats$parse_load(data$load_percentage)
      )

      p_features |>
        left_join(load_curves, by = "id") |>
        mutate(p_load_curve_deviation = .data$p_load_curve_slope - global_slope)
    },

    .load_sensitivity = function(velocity, load_pct) {
      load_numeric <- private$.stats$parse_load(load_pct)
      if (length(unique(load_numeric)) < 2) return(NA_real_)
      slope <- private$.stats$slope(velocity, load_numeric)
      if (is.na(slope)) return(NA_real_)
      slope * 10
    },

    .load_range = function(load_pct) {
      load_numeric <- private$.stats$parse_load(load_pct)
      max(load_numeric, na.rm = TRUE) - min(load_numeric, na.rm = TRUE)
    },

    .safe_ratio = function(num, denom) {
      ifelse(is.na(denom) | abs(denom) < 1e-10, NA_real_, num / denom)
    }
  )
)

# =============================================================================
# CALCULATOR: Cross-Level Features (Groups 4, 5, 6)
# =============================================================================

CrossLevelFeatureCalculator <- R6Class(
  classname = "CrossLevelFeatureCalculator",
  cloneable = FALSE,

  public = list(
    initialize = function() {
      private$.stats <- StatUtils$new()
    },

    #' @description Calculate participant × set features
    calculate_participant_set = function(data, participant_features) {
      stats <- private$.stats
      ps <- data |>
        group_by(.data$id, .data$set_id) |>
        summarize(
          ps_mean_velocity = mean(.data$mean_velocity, na.rm = TRUE),
          ps_velocity_cv = stats$cv(.data$mean_velocity),
          ps_reps_count = n(),
          ps_velocity_range = max(.data$mean_velocity, na.rm = TRUE) -
                              min(.data$mean_velocity, na.rm = TRUE),
          ps_decay_slope = stats$slope(.data$mean_velocity, .data$rep_number),
          .groups = "drop"
        ) |>
        left_join(participant_features[, c("id", "p_mean_velocity", "p_velocity_sd")], by = "id") |>
        mutate(
          ps_deviation_from_p_mean = .data$ps_mean_velocity - .data$p_mean_velocity,
          ps_z_vs_participant = private$.safe_z(.data$ps_mean_velocity, .data$p_mean_velocity, .data$p_velocity_sd)
        ) |>
        group_by(.data$id) |>
        mutate(
          ps_rank_among_sets = rank(-.data$ps_mean_velocity),
          ps_is_best_set = .data$ps_mean_velocity == max(.data$ps_mean_velocity, na.rm = TRUE),
          ps_is_worst_set = .data$ps_mean_velocity == min(.data$ps_mean_velocity, na.rm = TRUE),
          ps_ratio_to_first_set = .data$ps_mean_velocity / first(.data$ps_mean_velocity),
          ps_ratio_to_best_set = .data$ps_mean_velocity / max(.data$ps_mean_velocity, na.rm = TRUE)
        ) |>
        ungroup() |>
        select(-any_of(c("p_mean_velocity", "p_velocity_sd")))

      ps
    },

    #' @description Calculate participant × load features
    calculate_participant_load = function(data) {
      stats <- private$.stats
      global_by_load <- data |>
        group_by(.data$load_percentage) |>
        summarize(
          global_mean_at_load = mean(.data$mean_velocity, na.rm = TRUE),
          global_sd_at_load = sd(.data$mean_velocity, na.rm = TRUE),
          .groups = "drop"
        )

      data |>
        group_by(.data$id, .data$load_percentage) |>
        summarize(
          pl_mean_velocity = mean(.data$mean_velocity, na.rm = TRUE),
          pl_velocity_sd = sd(.data$mean_velocity, na.rm = TRUE),
          pl_velocity_cv = stats$cv(.data$mean_velocity),
          pl_n_observations = n(),
          pl_n_sets = n_distinct(.data$set_id),
          load_numeric = stats$parse_load(first(.data$load_percentage)),
          .groups = "drop"
        ) |>
        left_join(global_by_load, by = "load_percentage") |>
        mutate(
          pl_ratio_to_global_at_load = private$.safe_ratio(.data$pl_mean_velocity, .data$global_mean_at_load),
          pl_z_vs_population_at_load = private$.safe_z(.data$pl_mean_velocity, .data$global_mean_at_load, .data$global_sd_at_load)
        ) |>
        select(-any_of(c("global_mean_at_load", "global_sd_at_load")))
    },

    #' @description Calculate participant × RIR features
    calculate_participant_rir = function(data) {
      stats <- private$.stats
      global_by_rir <- data |>
        group_by(.data$rir) |>
        summarize(
          global_mean_at_rir = mean(.data$mean_velocity, na.rm = TRUE),
          global_sd_at_rir = sd(.data$mean_velocity, na.rm = TRUE),
          .groups = "drop"
        )

      data |>
        group_by(.data$id, .data$rir) |>
        summarize(
          pr_mean_velocity = mean(.data$mean_velocity, na.rm = TRUE),
          pr_velocity_sd = sd(.data$mean_velocity, na.rm = TRUE),
          pr_n_observations = n(),
          .groups = "drop"
        ) |>
        left_join(global_by_rir, by = "rir") |>
        mutate(
          pr_ratio_to_global_at_rir = private$.safe_ratio(.data$pr_mean_velocity, .data$global_mean_at_rir),
          pr_z_vs_population_at_rir = private$.safe_z(.data$pr_mean_velocity, .data$global_mean_at_rir, .data$global_sd_at_rir)
        ) |>
        select(-any_of(c("global_mean_at_rir", "global_sd_at_rir")))
    }
  ),

  private = list(
    .stats = NULL,

    .safe_ratio = function(num, denom) {
      ifelse(is.na(denom) | abs(denom) < 1e-10, NA_real_, num / denom)
    },

    .safe_z = function(x, mean, sd) {
      ifelse(is.na(sd) | sd < 1e-10, 0, (x - mean) / sd)
    }
  )
)

# =============================================================================
# CALCULATOR: Observation-Level Features (Groups 1, 10, 11, 12)
# =============================================================================

ObservationFeatureCalculator <- R6Class(
  classname = "ObservationFeatureCalculator",
  cloneable = FALSE,

  public = list(
    initialize = function() {
      private$.stats <- StatUtils$new()
    },

    #' @description Calculate observation-level features
    calculate = function(data, set_features, participant_features,
                         participant_load_features, participant_rir_features, global_stats) {
      obs <- private$.add_global_z_scores(data, global_stats)
      obs <- private$.add_set_context(obs, set_features)
      obs <- private$.add_participant_context(obs, participant_features)
      obs <- private$.add_session_features(obs)
      obs <- private$.add_derived_features(obs)
      obs <- private$.add_interaction_features(obs)
      obs <- private$.add_cross_level_context(obs, participant_load_features, participant_rir_features)
      obs
    }
  ),

  private = list(
    .stats = NULL,

    .add_global_z_scores = function(data, global_stats) {
      data |>
        mutate(
          velocity_z_global = (.data$mean_velocity - global_stats$mean_velocity) / global_stats$sd_velocity,
          velocity_ratio_to_global = .data$mean_velocity / global_stats$mean_velocity,
          rir_z_global = (.data$rir - global_stats$mean_rir) / global_stats$sd_rir,
          set_position_normalized = .data$rep_number / .data$reps_to_failure,
          is_first_rep = .data$rep_number == 1,
          is_last_rep = .data$rep_number == .data$reps_to_failure,
          load_numeric = private$.stats$parse_load(.data$load_percentage)
        )
    },

    .add_set_context = function(obs, set_features) {
      obs |>
        left_join(
          set_features[, c("id", "set_id", "set_mean_velocity", "set_max_velocity",
                           "set_first_rep_velocity", "set_velocity_sd", "set_sequence_in_session")],
          by = c("id", "set_id")
        ) |>
        mutate(
          velocity_z_within_set = private$.safe_z(.data$mean_velocity, .data$set_mean_velocity, .data$set_velocity_sd),
          velocity_ratio_to_set_max = private$.safe_ratio(.data$mean_velocity, .data$set_max_velocity),
          velocity_ratio_to_set_first = private$.safe_ratio(.data$mean_velocity, .data$set_first_rep_velocity),
          velocity_decay_from_start = .data$mean_velocity - .data$set_first_rep_velocity
        ) |>
        group_by(.data$set_id) |>
        mutate(velocity_rank_in_set = rank(-.data$mean_velocity)) |>
        ungroup()
    },

    .add_participant_context = function(obs, participant_features) {
      obs |>
        left_join(
          participant_features[, c("id", "p_mean_velocity", "p_velocity_sd", "p_rir_slope", "p_rir_intercept")],
          by = "id"
        ) |>
        mutate(
          velocity_z_within_participant = private$.safe_z(.data$mean_velocity, .data$p_mean_velocity, .data$p_velocity_sd),
          velocity_deviation_from_participant = .data$mean_velocity - .data$p_mean_velocity
        )
    },

    .add_session_features = function(obs) {
      obs |>
        group_by(.data$id) |>
        mutate(
          rep_number_in_session = row_number(),
          cumulative_volume = cumsum(1) - 1
        ) |>
        ungroup()
    },

    .add_derived_features = function(obs) {
      obs |>
        mutate(
          velocity_predicted_from_rir = .data$p_rir_intercept + .data$p_rir_slope * .data$rir,
          velocity_residual_from_rir = .data$mean_velocity - .data$velocity_predicted_from_rir,
          standardized_residual_rir = private$.safe_z(
            .data$velocity_residual_from_rir, 0, .data$p_velocity_sd
          ),
          velocity_ratio_to_expected = private$.safe_ratio(
            .data$mean_velocity, .data$velocity_predicted_from_rir
          )
        ) |>
        select(-any_of(c("p_rir_intercept", "velocity_predicted_from_rir")))
    },

    .add_interaction_features = function(obs) {
      obs |>
        mutate(
          rir_x_load_interaction = .data$rir * .data$load_numeric,
          set_position_x_rir = .data$set_position_normalized * .data$rir,
          cumulative_volume_x_velocity = .data$cumulative_volume * .data$mean_velocity
        )
    },

    .add_cross_level_context = function(obs, pl_features, pr_features) {
      obs |>
        left_join(
          pl_features[, c("id", "load_percentage", "pl_z_vs_population_at_load", "pl_ratio_to_global_at_load")],
          by = c("id", "load_percentage")
        ) |>
        left_join(
          pr_features[, c("id", "rir", "pr_z_vs_population_at_rir", "pr_ratio_to_global_at_rir")],
          by = c("id", "rir")
        )
    },

    .safe_ratio = function(num, denom) {
      ifelse(is.na(denom) | abs(denom) < 1e-10, NA_real_, num / denom)
    },

    .safe_z = function(x, mean, sd) {
      ifelse(is.na(sd) | sd < 1e-10, 0, (x - mean) / sd)
    }
  )
)

# =============================================================================
# FACADE: Feature Engineer
# =============================================================================

#' Comprehensive Feature Engineer
#'
#' Orchestrates feature creation across 12 feature groups using specialized calculators.
#' Follows facade pattern for simple interface, delegation for implementation.
#'
#' @export
FeatureEngineer <- R6Class(
  classname = "FeatureEngineer",
  cloneable = FALSE,

  public = list(
    #' @description Create feature engineer with all calculators
    initialize = function() {
      private$.global_calc <- GlobalStatsCalculator$new()
      private$.set_calc <- SetFeatureCalculator$new()
      private$.participant_calc <- ParticipantFeatureCalculator$new()
      private$.cross_level_calc <- CrossLevelFeatureCalculator$new()
      private$.observation_calc <- ObservationFeatureCalculator$new()
    },

    #' @description Engineer all features
    #' @param data Data frame with: id, set_id, rep_number, mean_velocity, rir, load_percentage, reps_to_failure
    #' @return FeatureEngineeringResult
    engineer_features = function(data) {
      private$.validate_columns(data)

      # Calculate global statistics
      global_stats <- private$.global_calc$calculate(data)

      # Calculate set features
      set_features <- private$.set_calc$calculate(data)

      # Calculate participant features
      participant_features <- private$.participant_calc$calculate(data, global_stats)

      # Update global stats with participant-derived values
      global_stats <- private$.global_calc$update_from_participants(global_stats, participant_features)

      # Calculate cross-level features
      ps_features <- private$.cross_level_calc$calculate_participant_set(data, participant_features)
      pl_features <- private$.cross_level_calc$calculate_participant_load(data)
      pr_features <- private$.cross_level_calc$calculate_participant_rir(data)

      # Calculate observation features
      obs_features <- private$.observation_calc$calculate(
        data, set_features, participant_features, pl_features, pr_features, global_stats
      )

      FeatureEngineeringResult$new(
        observation_features = obs_features,
        set_features = set_features,
        participant_features = participant_features,
        participant_set_features = ps_features,
        participant_load_features = pl_features,
        participant_rir_features = pr_features,
        feature_names = private$.get_feature_names(),
        global_stats = global_stats
      )
    },

    #' @description Get recommended features for observation-level anomaly detection
    get_observation_feature_names = function() {
      c("velocity_z_global", "velocity_z_within_participant", "velocity_z_within_set",
        "velocity_ratio_to_set_max", "velocity_decay_from_start", "set_position_normalized",
        "velocity_residual_from_rir", "standardized_residual_rir",
        "rir_x_load_interaction", "set_position_x_rir")
    },

    #' @description Get recommended features for participant-level anomaly detection
    get_participant_feature_names = function() {
      c("p_mean_velocity", "p_velocity_cv", "p_velocity_range", "p_velocity_iqr",
        "p_rir_slope", "p_rir_r2", "p_mean_reps_per_set", "p_reps_cv",
        "p_load_sensitivity", "p_velocity_z_global", "p_cv_ratio_to_population",
        "p_load_curve_slope", "p_load_curve_r2")
    }
  ),

  private = list(
    .global_calc = NULL,
    .set_calc = NULL,
    .participant_calc = NULL,
    .cross_level_calc = NULL,
    .observation_calc = NULL,

    .validate_columns = function(data) {
      required <- c("id", "set_id", "rep_number", "mean_velocity", "rir",
                    "load_percentage", "reps_to_failure")
      missing <- setdiff(required, names(data))
      if (length(missing) > 0) {
        stop("Missing required columns: ", paste(missing, collapse = ", "))
      }
    },

    .get_feature_names = function() {
      list(
        observation = c(
          "velocity_z_global", "velocity_ratio_to_global", "rir_z_global",
          "set_position_normalized", "is_first_rep", "is_last_rep",
          "velocity_z_within_set", "velocity_ratio_to_set_max",
          "velocity_ratio_to_set_first", "velocity_decay_from_start",
          "velocity_rank_in_set", "velocity_z_within_participant",
          "velocity_deviation_from_participant",
          "rep_number_in_session", "cumulative_volume",
          "velocity_residual_from_rir", "standardized_residual_rir",
          "velocity_ratio_to_expected",
          "rir_x_load_interaction", "set_position_x_rir", "cumulative_volume_x_velocity",
          "pl_z_vs_population_at_load", "pl_ratio_to_global_at_load",
          "pr_z_vs_population_at_rir", "pr_ratio_to_global_at_rir"
        ),
        set = c(
          "set_mean_velocity", "set_median_velocity", "set_velocity_sd", "set_velocity_cv",
          "set_velocity_range", "set_velocity_iqr", "set_max_velocity", "set_min_velocity",
          "set_first_rep_velocity", "set_last_rep_velocity",
          "set_velocity_drop", "set_velocity_drop_pct",
          "set_decay_slope", "set_decay_r2", "set_reps_count",
          "consecutive_set_velocity_ratio", "consecutive_set_velocity_diff",
          "set_to_session_first_ratio", "set_to_session_best_ratio",
          "cumulative_volume_at_set", "inter_set_recovery"
        ),
        participant = c(
          "p_mean_velocity", "p_median_velocity", "p_velocity_sd",
          "p_velocity_cv", "p_velocity_iqr", "p_velocity_var",
          "p_velocity_range", "p_max_velocity", "p_min_velocity",
          "p_percentile_10", "p_percentile_90",
          "p_velocity_skewness", "p_velocity_kurtosis",
          "p_rir_slope", "p_rir_intercept", "p_rir_r2", "p_rir_range", "p_mean_rir",
          "p_rir_slope_deviation",
          "p_n_sets", "p_n_observations", "p_mean_reps_per_set",
          "p_reps_cv", "p_max_reps_in_set", "p_min_reps_in_set", "p_reps_range",
          "p_load_sensitivity", "p_load_range", "p_mean_load", "p_n_loads",
          "p_velocity_ratio_to_population", "p_velocity_z_global",
          "p_cv_ratio_to_population", "p_rir_sensitivity_ratio",
          "p_load_sensitivity_ratio", "p_velocity_percentile", "p_consistency_percentile",
          "p_load_curve_slope", "p_load_curve_intercept", "p_load_curve_r2", "p_load_curve_deviation"
        ),
        participant_set = c(
          "ps_mean_velocity", "ps_velocity_cv", "ps_reps_count", "ps_velocity_range",
          "ps_decay_slope", "ps_deviation_from_p_mean", "ps_z_vs_participant",
          "ps_rank_among_sets", "ps_is_best_set", "ps_is_worst_set",
          "ps_ratio_to_first_set", "ps_ratio_to_best_set"
        ),
        participant_load = c(
          "pl_mean_velocity", "pl_velocity_sd", "pl_velocity_cv",
          "pl_n_observations", "pl_n_sets",
          "pl_ratio_to_global_at_load", "pl_z_vs_population_at_load"
        ),
        participant_rir = c(
          "pr_mean_velocity", "pr_velocity_sd", "pr_n_observations",
          "pr_ratio_to_global_at_rir", "pr_z_vs_population_at_rir"
        )
      )
    }
  )
)
