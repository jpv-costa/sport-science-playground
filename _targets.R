# _targets.R
# Pipeline Definition for Deadlift Study Meta-Analysis
#
# SOLID Principles Applied:
# - SRP: Each target has a single responsibility
# - DIP: Targets depend on abstractions (R6 classes) not implementations
#
# Usage:
#   targets::tar_make()           # Run pipeline
#   targets::tar_visnetwork()     # Visualize DAG
#   targets::tar_progress()       # Check status

library(targets)
library(tarchetypes)

# ==============================================================================
# Configuration
# ==============================================================================

# Set box module path (required for box::use to find R/ modules)
options(box.path = getwd())

# Load project modules via box
box::use(
  R/utils/config[load_config, get_param],
  R/utils/seed[set_reproducible_seed, get_default_seed],
  R/utils/logging[setup_logging],
  R/domain/study[Study],
  R/domain/treatment_group[TreatmentGroup],
  R/domain/effect_size[EffectSize],
  R/calculators/effect_size_calculator[SMCRCalculator, ROMCCalculator]
)

# Set global target options
tar_option_set(
  packages = c(
    "data.table",
    "arrow",
    "readxl",
    "metafor",
    "R6"
  ),
  format = "qs",                    # Fast serialization
  memory = "transient",             # Memory efficient
  garbage_collection = TRUE,
  error = "continue"                # Continue on error, mark as errored
)

# ==============================================================================
# Helper Functions (used by targets)
# ==============================================================================

#' Load and parse raw Excel data
#' @param file_path Path to Excel file
#' @param config Configuration object
#' @return data.table with raw study data
load_raw_data <- function(file_path, config) {
  data <- readxl::read_excel(file_path)
  data.table::as.data.table(data)
}

#' Create TreatmentGroup objects from raw data
#' @param raw_data data.table with study data
#' @return List of TreatmentGroup objects (valid only)
create_treatment_groups <- function(raw_data) {
  # Create groups from data with correct column mapping
  groups <- lapply(seq_len(nrow(raw_data)), function(i) {
    row <- raw_data[i]
    # Use avg.rir for RIR, pre.sd for SD, etc.
    tryCatch({
      TreatmentGroup$new(
        group_id = paste(row$study, row$group, sep = "_"),
        repetitions_in_reserve = row$avg.rir,
        mean_pre = row$pre.mean,
        mean_post = row$post.mean,
        standard_deviation_pre = row$pre.sd,
        sample_size = row$n,
        pre_post_correlation = if (is.na(row$pre.post.correlation)) NULL else row$pre.post.correlation
      )
    }, error = function(e) NULL)
  })

  # Filter out NULL (failed) groups
  Filter(Negate(is.null), groups)
}

#' Calculate effect sizes for all treatment groups
#' @param treatment_groups List of TreatmentGroup objects
#' @param effect_type "smcr" or "romc"
#' @param config Configuration object
#' @return List of EffectSize objects (valid only)
calculate_all_effect_sizes <- function(treatment_groups, effect_type, config) {
  calculator <- if (effect_type == "smcr") {
    SMCRCalculator$new(
      default_correlation = get_param(config, "effect_size", "default_correlation"),
      apply_bias_correction = get_param(config, "effect_size", "hedges_correction")
    )
  } else {
    ROMCCalculator$new(
      default_correlation = get_param(config, "effect_size", "default_correlation")
    )
  }

  # Filter to valid groups first
  valid_groups <- Filter(function(g) g$validate()$is_valid, treatment_groups)
  message(sprintf("Processing %d valid groups out of %d total", length(valid_groups), length(treatment_groups)))

  calculator$calculate_batch(valid_groups)
}

#' Convert EffectSize objects to metafor-compatible data.table
#' @param effect_sizes List of EffectSize objects
#' @param raw_data Original data for metadata
#' @return data.table ready for metafor
prepare_meta_data <- function(effect_sizes, raw_data) {
  effect_lists <- lapply(effect_sizes, function(es) es$to_list())
  effect_dt <- data.table::rbindlist(effect_lists)

  # Combine with study metadata (using correct column names from data)
  # Filter to only valid rows (matching number of effect sizes)
  n_effects <- nrow(effect_dt)
  metadata <- raw_data[1:n_effects, .(
    study,
    group,
    outcome_type = outcome,
    rir = avg.rir,
    training_status = train.status
  )]

  cbind(metadata, effect_dt)
}

# ==============================================================================
# Pipeline Definition
# ==============================================================================

list(
  # --------------------------------------------------------------------------
  # Stage 1: Configuration
  # --------------------------------------------------------------------------

  tar_target(
    name = config_file,
    command = "config/default.yml",
    format = "file"
  ),

  tar_target(
    name = analysis_config,
    command = load_config(config_file)
  ),

  tar_target(
    name = random_seed,
    command = get_default_seed(analysis_config)
  ),

  # --------------------------------------------------------------------------
  # Stage 2: Data Loading
  # --------------------------------------------------------------------------

  tar_target(
    name = raw_data_file,
    command = file.path(
      get_param(analysis_config, "data", "raw_path"),
      get_param(analysis_config, "data", "primary_dataset")
    ),
    format = "file"
  ),

  tar_target(
    name = raw_study_data,
    command = load_raw_data(raw_data_file, analysis_config)
  ),

  # --------------------------------------------------------------------------
  # Stage 3: Data Transformation
  # --------------------------------------------------------------------------

  tar_target(
    name = treatment_groups,
    command = create_treatment_groups(raw_study_data)
  ),

  tar_target(
    name = effect_sizes,
    command = calculate_all_effect_sizes(
      treatment_groups,
      effect_type = get_param(analysis_config, "effect_size", "type"),
      config = analysis_config
    )
  ),

  tar_target(
    name = meta_analysis_data,
    command = prepare_meta_data(effect_sizes, raw_study_data)
  ),

  # --------------------------------------------------------------------------
  # Stage 4: Save Processed Data
  # --------------------------------------------------------------------------

  tar_target(
    name = processed_data_file,
    command = {
      output_path <- file.path(
        get_param(analysis_config, "data", "processed_path"),
        "meta_analysis_data.parquet"
      )
      arrow::write_parquet(meta_analysis_data, output_path)
      output_path
    },
    format = "file"
  ),

  # --------------------------------------------------------------------------
  # Stage 5: Meta-Analysis Models (to be expanded)
  # --------------------------------------------------------------------------

  # Placeholder for model fitting targets
  # These will be added as we implement the model classes

  tar_target(
    name = pipeline_complete,
    command = {
      message("Pipeline completed successfully!")
      list(
        n_studies = nrow(raw_study_data),
        n_effect_sizes = length(effect_sizes),
        output_file = processed_data_file
      )
    }
  )
)
