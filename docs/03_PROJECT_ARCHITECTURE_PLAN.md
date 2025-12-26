# Deadlift Study: Project Architecture & Infrastructure Plan

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [Guiding Principles](#2-guiding-principles)
3. [Technology Stack Decisions](#3-technology-stack-decisions)
4. [Project Structure](#4-project-structure)
5. [Naming Conventions](#5-naming-conventions)
6. [R6 Class Architecture](#6-r6-class-architecture)
7. [Pipeline Orchestration](#7-pipeline-orchestration)
8. [Configuration Management](#8-configuration-management)
9. [Data Validation Framework](#9-data-validation-framework)
10. [Testing Strategy](#10-testing-strategy)
11. [Environment Reproducibility](#11-environment-reproducibility)
12. [Reporting & Documentation](#12-reporting--documentation)
13. [Decision Charts](#13-decision-charts)
14. [Implementation Roadmap](#14-implementation-roadmap)
15. [Sources & References](#15-sources--references)

---

## 1. Executive Summary

This document defines the architecture, tooling, and principles for the Deadlift Study project - a meta-analysis examining the relationship between Repetitions in Reserve (RIR) and resistance training outcomes (strength/hypertrophy).

**Key Design Goals:**
- **Reproducibility**: Exact replication of analyses across environments and time
- **Modularity**: Composable components following CUPID principles
- **Testability**: Comprehensive testing including property-based tests for statistical code
- **Scalability**: Handle multiple studies, hypotheses, and large datasets
- **Maintainability**: Clear structure following OOP and clean architecture principles

---

## 2. Guiding Principles

### 2.1 CUPID Principles (Adapted for R)

| Principle | Application |
|-----------|-------------|
| **Composable** | R6 classes with <10 public methods, dependency injection via constructor |
| **Unix** | Each class has one clear responsibility (e.g., `DataLoader`, `EffectSizeCalculator`) |
| **Predictable** | Random seeds parameterized, immutable configs, deterministic outputs |
| **Idiomatic** | Follow tidyverse style guide, use `box::use()` for imports |
| **Domain-Based** | Class names reflect meta-analysis domain, not implementation |

### 2.2 SOLID Principles Applied

| Principle | R Implementation |
|-----------|------------------|
| **SRP** | One actor per R6 class (e.g., `EffectSizeCalculator` serves statisticians) |
| **OCP** | Use R6 inheritance and composition for extension |
| **LSP** | Subclasses maintain parent contracts |
| **ISP** | Small, focused interfaces via R6 active bindings |
| **DIP** | Depend on abstractions (pass collaborators via constructor) |

### 2.3 Scientific Validity Tenets

1. **Reproducibility**: All random operations accept `seed` parameter
2. **SHAP Additivity**: Validate `shap_values.sum() + expected_value ≈ predictions`
3. **Numerical Stability**: Safe division, clipped values, log-space for extremes
4. **Data Leakage Prevention**: Fit on train, transform on test - never fit on test
5. **Statistical Correctness**: Use validated methods from established packages

---

## 3. Technology Stack Decisions

### 3.1 Decision Matrix

| Category | Choice | Rationale |
|----------|--------|-----------|
| **OOP Framework** | R6 | Encapsulated OOP, mutable objects, method chaining, aligns with CLAUDE.md OOP-first |
| **Pipeline** | targets | Automatic caching, dependency tracking, dynamic branching for multi-study |
| **Package Management** | renv + pak | renv for lockfiles, pak for fast installation |
| **R Version** | rig | Cross-platform, Rust-based, fast |
| **Data Validation** | pointblank | VALID-I/II/III workflows, Quarto integration |
| **Module System** | box | Python-like imports, explicit dependencies |
| **Logging** | lgr | Structured logging, multiple appenders |
| **Testing** | testthat + hedgehog | Unit tests + property-based testing |
| **Linting** | lintr + styler | PEP8-like consistency via pre-commit |
| **Config** | config + yaml | YAML-based, environment inheritance |
| **Reporting** | Quarto | Multi-language, parameterized reports |
| **Data Format** | Arrow/Parquet | Columnar, fast, cross-language |
| **Data Versioning** | DVC | Large file versioning, remote storage |
| **Containerization** | Rocker | Versioned R images, reproducible builds |

### 3.2 Package Categories

```r
# Core Analysis
c("metafor", "brms", "bayesmeta", "metaBMA", "RoBMA")

# Data Manipulation
c("data.table", "arrow", "dplyr", "tidyr")

# OOP & Architecture
c("R6", "box", "config", "yaml")

# Pipeline & Orchestration
c("targets", "tarchetypes", "crew")

# Validation & Testing
c("pointblank", "testthat", "hedgehog", "quickcheck", "assertr")

# Reporting
c("quarto", "ggplot2", "gt", "gtsummary")

# Dev Tools
c("renv", "pak", "lintr", "styler", "precommit", "lgr")
```

---

## 4. Project Structure

```
deadlift-study/
├── .R-version                    # rig: R version lock (e.g., "4.4.2")
├── .Rprofile                     # Project-specific R settings
├── .Renviron                     # Environment variables (secrets excluded from git)
├── .lintr                        # Linting configuration
├── .pre-commit-config.yaml       # Pre-commit hooks
├── .gitattributes                # Git LFS patterns
├── .gitignore                    # Standard ignores
├── .dvc/                         # DVC configuration
├── .dvcignore                    # DVC ignores
├── renv.lock                     # Package lockfile
├── renv/                         # renv library
├── Makefile                      # Automation commands
├── _targets.R                    # Pipeline definition
├── _quarto.yml                   # Quarto project config
│
├── config/                       # Configuration files
│   ├── default.yml               # Base configuration
│   ├── development.yml           # Dev overrides
│   ├── production.yml            # Production settings
│   └── studies/                  # Study-specific configs
│       ├── strength_rir.yml
│       └── hypertrophy_rir.yml
│
├── R/                            # Source code (box modules)
│   ├── __init__.R                # Package-level exports
│   ├── domain/                   # Domain entities
│   │   ├── study.R               # Study R6 class
│   │   ├── effect_size.R         # EffectSize R6 class
│   │   └── meta_result.R         # MetaResult R6 class
│   ├── loaders/                  # Data loading
│   │   ├── excel_loader.R        # ExcelLoader R6 class
│   │   ├── csv_loader.R          # CsvLoader R6 class
│   │   └── parquet_loader.R      # ParquetLoader R6 class
│   ├── calculators/              # Effect size calculation
│   │   ├── smcr_calculator.R     # SMCRCalculator R6 class
│   │   └── romc_calculator.R     # ROMCCalculator R6 class
│   ├── models/                   # Meta-analysis models
│   │   ├── base_model.R          # BaseMetaModel R6 class
│   │   ├── frequentist_model.R   # FrequentistModel R6 class
│   │   └── bayesian_model.R      # BayesianModel R6 class
│   ├── validators/               # Data validation
│   │   ├── study_validator.R     # StudyValidator R6 class
│   │   └── effect_validator.R    # EffectValidator R6 class
│   ├── reporters/                # Output generation
│   │   ├── forest_plot.R         # ForestPlotReporter R6 class
│   │   └── summary_table.R       # SummaryTableReporter R6 class
│   └── utils/                    # Shared utilities
│       ├── logging.R             # Logging setup
│       ├── config.R              # Config loading
│       └── seed.R                # Seed management
│
├── tests/                        # Test suite
│   ├── testthat.R                # Test runner
│   ├── testthat/
│   │   ├── setup.R               # Test fixtures
│   │   ├── test-smcr_calculator.R
│   │   ├── test-romc_calculator.R
│   │   ├── test-frequentist_model.R
│   │   └── properties/           # Property-based tests
│   │       ├── test-effect_size_properties.R
│   │       └── test-model_properties.R
│   └── fixtures/                 # Test data
│       └── sample_study.rds
│
├── data/                         # Data directory (DVC tracked)
│   ├── raw/                      # Original data (immutable)
│   │   └── pelland_dataset.xlsx
│   ├── processed/                # Cleaned data
│   │   └── studies.parquet
│   └── external/                 # External datasets
│
├── analyses/                     # Analysis scripts (targets entry points)
│   ├── 01_primary_analysis.qmd
│   ├── 02_sensitivity_analysis.qmd
│   ├── 03_publication_bias.qmd
│   └── 04_bayesian_analysis.qmd
│
├── reports/                      # Generated reports (gitignored)
│   ├── primary/
│   └── sensitivity/
│
├── docs/                         # Documentation
│   ├── 01_CURRENT_EXPERIMENT_PLAN.md
│   ├── 02_METHODOLOGICAL_IMPROVEMENTS.md
│   └── 03_PROJECT_ARCHITECTURE_PLAN.md
│
├── docker/                       # Docker configuration
│   ├── Dockerfile
│   └── docker-compose.yml
│
└── scripts/                      # Utility scripts
    ├── setup.R                   # Initial setup
    └── check.R                   # CI checks
```

---

## 5. Naming Conventions

### 5.1 Files

| Type | Convention | Example |
|------|------------|---------|
| R6 Classes | `snake_case.R`, singular noun | `effect_size_calculator.R` |
| Functions | `snake_case.R`, verb phrase | `calculate_effect_size.R` |
| Tests | `test-{module}.R` | `test-smcr_calculator.R` |
| Configs | `snake_case.yml` | `strength_rir.yml` |
| Data | `snake_case.{ext}` | `processed_studies.parquet` |
| Analyses | `##_{description}.qmd` | `01_primary_analysis.qmd` |

### 5.2 R6 Classes

| Element | Convention | Example |
|---------|------------|---------|
| Class Name | `PascalCase`, singular noun | `EffectSizeCalculator` |
| Public Methods | `snake_case`, verb phrase | `calculate()`, `validate()` |
| Private Methods | `snake_case`, prefixed `._` | `._compute_variance()` |
| Active Bindings | `snake_case`, noun | `result`, `is_valid` |
| Constructor Args | `snake_case`, descriptive | `raw_data`, `config_path` |

### 5.3 Variables

| Type | Convention | Example |
|------|------------|---------|
| Data Frames | `snake_case`, plural noun | `studies`, `effect_sizes` |
| Single Values | `snake_case`, noun | `sample_size`, `mean_diff` |
| Booleans | `is_` or `has_` prefix | `is_valid`, `has_moderators` |
| Counts | `n_` prefix or `_count` suffix | `n_studies`, `observation_count` |
| Lists | `snake_case`, plural or `_list` | `models`, `result_list` |
| Matrices | `snake_case`, `_matrix` suffix | `correlation_matrix` |

### 5.4 Functions (Non-OOP)

| Category | Convention | Example |
|----------|------------|---------|
| Actions | verb + noun | `calculate_effect_size()` |
| Queries | verb + noun | `fetch_study()`, `find_outliers()` |
| Transformations | verb + noun | `transform_to_log()` |
| Predicates | `is_` prefix | `is_valid_effect_size()` |
| Factories | `create_` prefix | `create_meta_model()` |

### 5.5 Configuration Keys

```yaml
# Use snake_case for all keys
analysis:
  effect_size_type: "smcr"      # Not: effectSizeType
  min_sample_size: 10           # Not: minSampleSize
  random_state: 42              # Not: randomState

models:
  frequentist:
    method: "REML"
    test_type: "knha"           # Knapp-Hartung adjustment
```

### 5.6 Column Names (Data)

| Type | Convention | Example |
|------|------------|---------|
| IDs | `{entity}_id` | `study_id`, `group_id` |
| Measurements | descriptive, with units | `mean_pre`, `sd_post` |
| Counts | `n_` prefix | `n_participants`, `n_sets` |
| Categories | descriptive noun | `training_status`, `exercise_type` |
| Booleans | `is_` or `has_` prefix | `is_failure_trained` |
| Timestamps | `_at` suffix | `created_at`, `updated_at` |

---

## 6. R6 Class Architecture

### 6.1 Base Classes

```r
# R/domain/base_entity.R
box::use(R6[R6Class])

#' Base Entity
#'
#' Abstract base class for all domain entities.
#' Provides common validation and serialization.
BaseEntity <- R6Class(
  classname = "BaseEntity",
  public = list(
    #' @description Create a new entity
    #' @param id Unique identifier
    initialize = function(id = NULL) {
      private$.id <- id %||% private$.generate_id()
      private$.created_at <- Sys.time()
    },

    #' @description Validate entity state
    #' @return Logical, TRUE if valid
    validate = function() {
      stop("Subclass must implement validate()")
    },

    #' @description Convert to list for serialization
    #' @return Named list
    to_list = function() {
      stop("Subclass must implement to_list()")
    }
  ),

  active = list(
    #' @field id Entity identifier
    id = function() private$.id,

    #' @field created_at Creation timestamp
    created_at = function() private$.created_at
  ),

  private = list(
    .id = NULL,
    .created_at = NULL,

    .generate_id = function() {
      paste0(class(self)[1], "_", format(Sys.time(), "%Y%m%d%H%M%S"))
    }
  )
)
```

### 6.2 Domain Entity: Study

```r
# R/domain/study.R
box::use(
  R6[R6Class],
  ./base_entity[BaseEntity]
)

#' Study Entity
#'
#' Represents a single study in the meta-analysis.
Study <- R6Class(
  classname = "Study",
  inherit = BaseEntity,

  public = list(
    #' @description Create a new Study
    #' @param study_id Study identifier
    #' @param author First author name
    #' @param year Publication year
    #' @param outcome_type "strength" or "hypertrophy"
    #' @param groups List of treatment groups
    initialize = function(study_id, author, year, outcome_type, groups = list()) {
      super$initialize(id = study_id)
      private$.author <- author
      private$.year <- year
      private$.outcome_type <- outcome_type
      private$.groups <- groups
    },

    #' @description Validate study data
    validate = function() {
      errors <- character()

      if (is.na(private$.year) || private$.year < 1900 || private$.year > 2100) {
        errors <- c(errors, "Year must be between 1900 and 2100")
      }

      if (!private$.outcome_type %in% c("strength", "hypertrophy")) {
        errors <- c(errors, "Outcome type must be 'strength' or 'hypertrophy'")
      }

      if (length(errors) > 0) {
        stop(paste(errors, collapse = "; "))
      }

      invisible(TRUE)
    },

    #' @description Add a treatment group
    #' @param group TreatmentGroup object
    #' @return self (for chaining)
    add_group = function(group) {
      private$.groups <- c(private$.groups, list(group))
      invisible(self)
    },

    #' @description Convert to list
    to_list = function() {
      list(
        study_id = self$id,
        author = private$.author,
        year = private$.year,
        outcome_type = private$.outcome_type,
        n_groups = length(private$.groups)
      )
    }
  ),

  active = list(
    author = function() private$.author,
    year = function() private$.year,
    outcome_type = function() private$.outcome_type,
    groups = function() private$.groups,
    n_groups = function() length(private$.groups)
  ),

  private = list(
    .author = NULL,
    .year = NULL,
    .outcome_type = NULL,
    .groups = NULL
  )
)
```

### 6.3 Calculator: SMCR Effect Size

```r
# R/calculators/smcr_calculator.R
box::use(
  R6[R6Class],
  stats[...],
  ../utils/logging[get_logger]
)

#' SMCR Effect Size Calculator
#'
#' Calculates Standardized Mean Change using Raw score standardization.
#' Uses Hedge's g correction for small sample bias.
SMCRCalculator <- R6Class(
  classname = "SMCRCalculator",

  public = list(
    #' @description Create calculator with configuration
    #' @param correlation_method Method for pre-post correlation imputation
    initialize = function(correlation_method = "fisher_z") {
      private$.correlation_method <- correlation_method
      private$.log <- get_logger("SMCRCalculator")
    },

    #' @description Calculate SMCR effect size
    #' @param mean_pre Pre-intervention mean
    #' @param mean_post Post-intervention mean
    #' @param sd_pre Pre-intervention SD
    #' @param n Sample size
    #' @param correlation Pre-post correlation (optional)
    #' @return Named list with effect_size and variance
    calculate = function(mean_pre, mean_post, sd_pre, n, correlation = NULL) {
      private$.validate_inputs(mean_pre, mean_post, sd_pre, n)

      if (is.null(correlation)) {
        correlation <- private$.impute_correlation()
      }

      # Raw score standardization
      mean_change <- mean_post - mean_pre
      effect_size <- mean_change / sd_pre

      # Hedge's g correction (small sample bias)
      j <- private$.hedges_correction(n)
      effect_size_corrected <- effect_size * j

      # Variance calculation
      variance <- private$.calculate_variance(n, effect_size, correlation)
      variance_corrected <- variance * (j^2)

      private$.log$info(
        "Calculated SMCR: ES = {round(effect_size_corrected, 3)}, Var = {round(variance_corrected, 4)}"
      )

      list(
        effect_size = effect_size_corrected,
        variance = variance_corrected,
        se = sqrt(variance_corrected),
        n = n,
        correlation = correlation,
        method = "smcr_hedges_g"
      )
    }
  ),

  private = list(
    .correlation_method = NULL,
    .log = NULL,

    .validate_inputs = function(mean_pre, mean_post, sd_pre, n) {
      if (sd_pre <= 0) {
        stop("Standard deviation must be positive")
      }
      if (n < 2) {
        stop("Sample size must be at least 2")
      }
    },

    .impute_correlation = function() {
      # Using Fisher's z transformation pooled estimate
      # Default: r = 0.5 (moderate correlation)
      0.5
    },

    .hedges_correction = function(n) {
      # Small sample bias correction
      # J ≈ 1 - 3/(4*df - 1) where df = n - 1
      df <- n - 1
      1 - 3 / (4 * df - 1)
    },

    .calculate_variance = function(n, effect_size, correlation) {
      # Variance for repeated measures standardized mean change
      # Var(d) = (2(1-r)/n) + (d^2 / 2n)
      term1 <- 2 * (1 - correlation) / n
      term2 <- effect_size^2 / (2 * n)
      term1 + term2
    }
  )
)
```

### 6.4 Model: Frequentist Meta-Regression

```r
# R/models/frequentist_model.R
box::use(
  R6[R6Class],
  metafor[rma.mv, robust, predict.rma],
  ../utils/logging[get_logger]
)

#' Frequentist Meta-Regression Model
#'
#' Multi-level meta-regression using metafor::rma.mv
#' Supports nested random effects and robust variance estimation.
FrequentistModel <- R6Class(
  classname = "FrequentistModel",

  public = list(
    #' @description Create model with configuration
    #' @param random_structure Formula for random effects
    #' @param method Estimation method (default: "REML")
    #' @param robust Use robust variance estimation
    initialize = function(
      random_structure = ~ 1 | study/group/obs,
      method = "REML",
      robust = TRUE
    ) {
      private$.random_structure <- random_structure
      private$.method <- method
      private$.robust <- robust
      private$.log <- get_logger("FrequentistModel")
    },

    #' @description Fit the meta-regression model
    #' @param data Data frame with effect sizes and variances
    #' @param moderators Optional moderator formula
    #' @return self (for chaining)
    fit = function(data, moderators = NULL) {
      private$.validate_data(data)

      formula <- if (is.null(moderators)) ~ 1 else moderators

      private$.log$info("Fitting meta-regression with {nrow(data)} observations")

      private$.model <- rma.mv(
        yi = yi,
        V = vi,
        mods = formula,
        random = private$.random_structure,
        data = data,
        method = private$.method,
        sparse = TRUE
      )

      if (private$.robust) {
        private$.robust_results <- robust(
          private$.model,
          cluster = data$study
        )
      }

      private$.fitted <- TRUE
      private$.log$info("Model fitted successfully")

      invisible(self)
    },

    #' @description Get model summary
    #' @return List with coefficients, heterogeneity stats
    summary = function() {
      private$.check_fitted()

      model <- if (private$.robust) private$.robust_results else private$.model

      list(
        coefficients = coef(model),
        se = model$se,
        ci_lower = model$ci.lb,
        ci_upper = model$ci.ub,
        pval = model$pval,
        tau2 = private$.model$sigma2,
        i2 = private$.calculate_i2(),
        q_statistic = private$.model$QE,
        q_pval = private$.model$QEp,
        n_studies = private$.model$k
      )
    },

    #' @description Generate predictions
    #' @param newdata Optional new data for prediction
    #' @param level Confidence level (default: 0.95)
    #' @return Predictions with intervals
    predict = function(newdata = NULL, level = 0.95) {
      private$.check_fitted()

      predict.rma(
        private$.model,
        newmods = newdata,
        level = level,
        pi.type = "default"
      )
    }
  ),

  active = list(
    is_fitted = function() private$.fitted,
    model = function() private$.model,
    coefficients = function() if (private$.fitted) coef(private$.model) else NULL
  ),

  private = list(
    .random_structure = NULL,
    .method = NULL,
    .robust = NULL,
    .model = NULL,
    .robust_results = NULL,
    .fitted = FALSE,
    .log = NULL,

    .validate_data = function(data) {
      required <- c("yi", "vi", "study")
      missing <- setdiff(required, names(data))
      if (length(missing) > 0) {
        stop("Missing required columns: ", paste(missing, collapse = ", "))
      }
    },

    .check_fitted = function() {
      if (!private$.fitted) {
        stop("Model must be fitted before calling this method")
      }
    },

    .calculate_i2 = function() {
      # I² = (Q - df) / Q × 100%
      q <- private$.model$QE
      df <- private$.model$k - 1
      if (q <= df) return(0)
      100 * (q - df) / q
    }
  )
)
```

---

## 7. Pipeline Orchestration

### 7.1 targets Configuration

```r
# _targets.R
library(targets)
library(tarchetypes)

# Source R modules
box::use(
  R/loaders/excel_loader[ExcelLoader],
  R/calculators/smcr_calculator[SMCRCalculator],
  R/calculators/romc_calculator[ROMCCalculator],
  R/models/frequentist_model[FrequentistModel],
  R/validators/study_validator[StudyValidator],
  R/utils/config[load_config]
)

# Load configuration
config <- load_config()

# Set global options
tar_option_set(
  packages = c("data.table", "arrow", "metafor"),
  format = "qs",  # Fast serialization
  memory = "transient",  # Memory efficient
  garbage_collection = TRUE
)

# Define targets
list(
  # Configuration
  tar_target(
    config_file,
    "config/default.yml",
    format = "file"
  ),

  tar_target(
    analysis_config,
    load_config(config_file)
  ),

  # Data Loading
  tar_target(
    raw_data_file,
    "data/raw/pelland_dataset.xlsx",
    format = "file"
  ),

  tar_target(
    raw_studies,
    {
      loader <- ExcelLoader$new()
      loader$load(raw_data_file)
    }
  ),

  # Data Validation
  tar_target(
    validated_studies,
    {
      validator <- StudyValidator$new()
      validator$validate(raw_studies)
    }
  ),

  # Effect Size Calculation - Dynamic Branching
  tar_target(
    effect_sizes_smcr,
    {
      calculator <- SMCRCalculator$new()
      validated_studies |>
        split(by = "outcome_type") |>
        lapply(function(df) calculator$calculate_batch(df))
    },
    pattern = map(validated_studies)
  ),

  # Primary Analysis - Strength
  tar_target(
    model_strength,
    {
      model <- FrequentistModel$new(
        random_structure = ~ 1 | study/group/obs,
        robust = TRUE
      )
      strength_data <- effect_sizes_smcr[outcome_type == "strength"]
      model$fit(strength_data, moderators = ~ rir)
    }
  ),

  # Primary Analysis - Hypertrophy
  tar_target(
    model_hypertrophy,
    {
      model <- FrequentistModel$new(
        random_structure = ~ 1 | study/group/obs,
        robust = TRUE
      )
      hypertrophy_data <- effect_sizes_smcr[outcome_type == "hypertrophy"]
      model$fit(hypertrophy_data, moderators = ~ rir)
    }
  ),

  # Reports - Parameterized Quarto
  tar_quarto(
    primary_report,
    path = "analyses/01_primary_analysis.qmd",
    extra_files = c("_quarto.yml")
  )
)
```

### 7.2 Pipeline Commands (Makefile)

```makefile
# Makefile

.PHONY: all pipeline test lint check clean

# Default target
all: check pipeline

# Run full pipeline
pipeline:
	Rscript -e "targets::tar_make()"

# Run pipeline in parallel
pipeline-parallel:
	Rscript -e "targets::tar_make_future(workers = 4)"

# Visualize pipeline
pipeline-viz:
	Rscript -e "targets::tar_visnetwork()"

# Check pipeline status
pipeline-status:
	Rscript -e "targets::tar_progress()"

# Invalidate and rerun
pipeline-clean:
	Rscript -e "targets::tar_invalidate(everything())"

# Test suite
test:
	Rscript -e "testthat::test_local()"

test-coverage:
	Rscript -e "covr::package_coverage()"

# Linting
lint:
	Rscript -e "lintr::lint_package()"

style:
	Rscript -e "styler::style_pkg()"

# Pre-commit
check: lint test
	pre-commit run --all-files

# Clean artifacts
clean:
	rm -rf _targets/
	rm -rf reports/
```

---

## 8. Configuration Management

### 8.1 Base Configuration

```yaml
# config/default.yml
default:
  project:
    name: "Deadlift Study Meta-Analysis"
    version: "1.0.0"
    authors:
      - "Research Team"

  analysis:
    effect_size_type: "smcr"
    correlation_imputation_method: "fisher_z"
    default_correlation: 0.5
    min_sample_size: 5
    outlier_threshold_z: 3.0

  models:
    frequentist:
      method: "REML"
      test_type: "knha"
      random_structure: "~1|study/group/obs"
      robust_variance: true
    bayesian:
      chains: 4
      iterations: 4000
      warmup: 1000
      seed: 42

  validation:
    required_columns:
      - study_id
      - author
      - year
      - rir
      - mean_pre
      - mean_post
      - sd_pre
      - n
    value_ranges:
      rir:
        min: 0
        max: 20
      year:
        min: 1980
        max: 2025

  reporting:
    confidence_level: 0.95
    prediction_interval: true
    forest_plot_theme: "minimal"

  logging:
    level: "INFO"
    file: "logs/analysis.log"

  random:
    seed: 42

development:
  inherits: default
  logging:
    level: "DEBUG"

production:
  inherits: default
  models:
    bayesian:
      chains: 8
      iterations: 10000
```

### 8.2 Config Loader

```r
# R/utils/config.R
box::use(
  config[get],
  yaml[read_yaml],
  ../utils/logging[get_logger]
)

#' Load Configuration
#'
#' Loads configuration from YAML file with environment override.
#'
#' @param config_file Path to config file
#' @param environment Environment name (default: from R_CONFIG_ACTIVE)
#' @return Validated configuration list
load_config <- function(config_file = "config/default.yml",
                        environment = Sys.getenv("R_CONFIG_ACTIVE", "default")) {
  log <- get_logger("Config")

  log$info("Loading configuration from {config_file} for environment: {environment}")

  cfg <- config::get(file = config_file, config = environment)

  # Validate required keys
  required <- c("analysis", "models", "validation")
  missing <- setdiff(required, names(cfg))
  if (length(missing) > 0) {
    stop("Missing required config sections: ", paste(missing, collapse = ", "))
  }

  log$debug("Configuration loaded successfully")
  cfg
}

#' Get Analysis Parameter
#'
#' @param cfg Configuration object
#' @param ... Path to parameter (e.g., "models", "frequentist", "method")
#' @return Parameter value
get_param <- function(cfg, ...) {
  path <- c(...)
  value <- cfg
  for (key in path) {
    value <- value[[key]]
    if (is.null(value)) {
      stop("Config key not found: ", paste(path, collapse = "."))
    }
  }
  value
}
```

---

## 9. Data Validation Framework

### 9.1 pointblank Validation Agent

```r
# R/validators/study_validator.R
box::use(
  R6[R6Class],
  pointblank[...],
  ../utils/logging[get_logger]
)

#' Study Data Validator
#'
#' Validates study data using pointblank VALID-II workflow.
#' Produces validation reports and fails on critical errors.
StudyValidator <- R6Class(
  classname = "StudyValidator",

  public = list(
    #' @description Create validator with configuration
    #' @param config Validation configuration list
    initialize = function(config = NULL) {
      private$.config <- config %||% private$.default_config()
      private$.log <- get_logger("StudyValidator")
    },

    #' @description Validate study data frame
    #' @param data Data frame to validate
    #' @param stop_on_fail Stop execution on validation failure
    #' @return Validated data frame
    validate = function(data, stop_on_fail = TRUE) {
      private$.log$info("Starting validation for {nrow(data)} rows")

      agent <- create_agent(
        tbl = data,
        tbl_name = "study_data",
        label = "Study Data Validation"
      ) |>
        # Required columns exist
        col_exists(
          columns = private$.config$required_columns
        ) |>

        # No NA in critical columns
        col_vals_not_null(
          columns = c("study_id", "author", "year", "rir")
        ) |>

        # Value range checks
        col_vals_between(
          columns = "rir",
          left = private$.config$value_ranges$rir$min,
          right = private$.config$value_ranges$rir$max
        ) |>
        col_vals_between(
          columns = "year",
          left = private$.config$value_ranges$year$min,
          right = private$.config$value_ranges$year$max
        ) |>

        # Positive values where required
        col_vals_gt(
          columns = c("n", "sd_pre"),
          value = 0
        ) |>

        # Sample size minimum
        col_vals_gte(
          columns = "n",
          value = private$.config$min_sample_size
        ) |>

        # Execute validation
        interrogate()

      # Log results
      summary <- get_agent_report(agent)
      private$.log$info("Validation complete: {sum(summary$pass)} passed, {sum(summary$fail)} failed")

      # Check for critical failures
      if (all_passed(agent)) {
        private$.log$info("All validations passed")
      } else if (stop_on_fail) {
        stop("Validation failed. See report for details.")
      } else {
        private$.log$warn("Some validations failed but continuing")
      }

      # Store agent for reporting
      private$.agent <- agent

      invisible(data)
    },

    #' @description Generate validation report
    #' @param output_file Path for HTML report
    export_report = function(output_file = "reports/validation_report.html") {
      if (is.null(private$.agent)) {
        stop("Must run validate() before exporting report")
      }

      export_report(private$.agent, filename = output_file)
      private$.log$info("Validation report exported to {output_file}")
    }
  ),

  active = list(
    is_valid = function() {
      if (is.null(private$.agent)) return(NA)
      all_passed(private$.agent)
    },

    failure_count = function() {
      if (is.null(private$.agent)) return(NA)
      sum(get_agent_report(private$.agent)$fail)
    }
  ),

  private = list(
    .config = NULL,
    .agent = NULL,
    .log = NULL,

    .default_config = function() {
      list(
        required_columns = c("study_id", "author", "year", "rir",
                            "mean_pre", "mean_post", "sd_pre", "n"),
        value_ranges = list(
          rir = list(min = 0, max = 20),
          year = list(min = 1980, max = 2025)
        ),
        min_sample_size = 5
      )
    }
  )
)
```

---

## 10. Testing Strategy

### 10.1 Test Structure

```
tests/
├── testthat.R                    # Test runner
├── testthat/
│   ├── setup.R                   # Global fixtures
│   │
│   ├── test-smcr_calculator.R    # Unit tests
│   ├── test-romc_calculator.R
│   ├── test-frequentist_model.R
│   ├── test-study_validator.R
│   │
│   ├── integration/              # Integration tests
│   │   ├── test-pipeline.R
│   │   └── test-end_to_end.R
│   │
│   └── properties/               # Property-based tests
│       ├── test-effect_size_properties.R
│       └── test-model_invariants.R
│
└── fixtures/                     # Test data
    ├── sample_study.rds
    └── expected_results.rds
```

### 10.2 Unit Test Example

```r
# tests/testthat/test-smcr_calculator.R
box::use(
  testthat[...],
  ../../R/calculators/smcr_calculator[SMCRCalculator]
)

describe("SMCRCalculator", {
  calculator <- SMCRCalculator$new()

  describe("calculate()", {
    it("returns correct effect size for known values", {
      result <- calculator$calculate(
        mean_pre = 100,
        mean_post = 110,
        sd_pre = 15,
        n = 30,
        correlation = 0.5
      )

      # Expected: (110-100)/15 * J ≈ 0.667 * 0.974 ≈ 0.650
      expect_equal(result$effect_size, 0.650, tolerance = 0.01)
    })

    it("applies Hedge's g correction for small samples", {
      small_n <- calculator$calculate(
        mean_pre = 100, mean_post = 110, sd_pre = 15, n = 10
      )
      large_n <- calculator$calculate(
        mean_pre = 100, mean_post = 110, sd_pre = 15, n = 100
      )

      # Small sample effect size should be smaller (more correction)
      expect_lt(small_n$effect_size, large_n$effect_size)
    })

    it("rejects non-positive standard deviation", {
      expect_error(
        calculator$calculate(mean_pre = 100, mean_post = 110, sd_pre = 0, n = 30),
        "Standard deviation must be positive"
      )
    })

    it("rejects sample size less than 2", {
      expect_error(
        calculator$calculate(mean_pre = 100, mean_post = 110, sd_pre = 15, n = 1),
        "Sample size must be at least 2"
      )
    })
  })
})
```

### 10.3 Property-Based Test Example

```r
# tests/testthat/properties/test-effect_size_properties.R
box::use(
  testthat[...],
  hedgehog[...],
  ../../R/calculators/smcr_calculator[SMCRCalculator]
)

describe("Effect Size Properties", {
  calculator <- SMCRCalculator$new()

  # Generator for valid effect size inputs
  gen_valid_inputs <- gen.and_then(
    gen.list(
      mean_pre = gen.double(min = 1, max = 1000),
      mean_post = gen.double(min = 1, max = 1000),
      sd_pre = gen.double(min = 0.1, max = 100),
      n = gen.int(min = 5, max = 500),
      correlation = gen.double(min = 0.01, max = 0.99)
    )
  )

  it("effect size is positive when post > pre", {
    forall(gen_valid_inputs, function(inputs) {
      # Ensure post > pre for this property
      inputs$mean_post <- inputs$mean_pre + abs(inputs$mean_post - inputs$mean_pre) + 0.01

      result <- calculator$calculate(
        mean_pre = inputs$mean_pre,
        mean_post = inputs$mean_post,
        sd_pre = inputs$sd_pre,
        n = inputs$n,
        correlation = inputs$correlation
      )

      expect_true(result$effect_size > 0)
    })
  })

  it("variance is always positive", {
    forall(gen_valid_inputs, function(inputs) {
      result <- calculator$calculate(
        mean_pre = inputs$mean_pre,
        mean_post = inputs$mean_post,
        sd_pre = inputs$sd_pre,
        n = inputs$n,
        correlation = inputs$correlation
      )

      expect_true(result$variance > 0)
    })
  })

  it("effect size scales linearly with mean difference", {
    forall(gen_valid_inputs, function(inputs) {
      result1 <- calculator$calculate(
        mean_pre = inputs$mean_pre,
        mean_post = inputs$mean_pre + 10,
        sd_pre = inputs$sd_pre,
        n = inputs$n,
        correlation = inputs$correlation
      )

      result2 <- calculator$calculate(
        mean_pre = inputs$mean_pre,
        mean_post = inputs$mean_pre + 20,
        sd_pre = inputs$sd_pre,
        n = inputs$n,
        correlation = inputs$correlation
      )

      # Effect size should approximately double
      ratio <- result2$effect_size / result1$effect_size
      expect_equal(ratio, 2, tolerance = 0.01)
    })
  })

  it("Hedge's g correction decreases with larger n", {
    forall(gen_valid_inputs, function(inputs) {
      small_n <- calculator$calculate(
        mean_pre = inputs$mean_pre,
        mean_post = inputs$mean_post,
        sd_pre = inputs$sd_pre,
        n = 10,
        correlation = inputs$correlation
      )

      large_n <- calculator$calculate(
        mean_pre = inputs$mean_pre,
        mean_post = inputs$mean_post,
        sd_pre = inputs$sd_pre,
        n = 100,
        correlation = inputs$correlation
      )

      # Correction reduces absolute effect size for small n
      # So |small_n$effect_size| <= |large_n$effect_size|
      expect_lte(abs(small_n$effect_size), abs(large_n$effect_size) + 0.01)
    })
  })
})
```

### 10.4 Numerical Precision Tests

```r
# tests/testthat/test-numerical_stability.R
box::use(
  testthat[...],
  ../../R/calculators/smcr_calculator[SMCRCalculator]
)

describe("Numerical Stability", {
  calculator <- SMCRCalculator$new()

  it("handles very small effect sizes without underflow", {
    result <- calculator$calculate(
      mean_pre = 100,
      mean_post = 100.001,
      sd_pre = 15,
      n = 30
    )

    expect_true(is.finite(result$effect_size))
    expect_true(is.finite(result$variance))
    expect_true(result$variance > 0)
  })

  it("handles very large values without overflow", {
    result <- calculator$calculate(
      mean_pre = 1e6,
      mean_post = 1.1e6,
      sd_pre = 1e5,
      n = 30
    )

    expect_true(is.finite(result$effect_size))
    expect_true(is.finite(result$variance))
  })

  it("produces consistent results across platforms (reproducibility)", {
    # Known input with expected output (verified manually)
    result <- calculator$calculate(
      mean_pre = 100,
      mean_post = 115,
      sd_pre = 20,
      n = 50,
      correlation = 0.5
    )

    # Expected: (115-100)/20 * J(50) = 0.75 * 0.985 ≈ 0.739
    expect_equal(result$effect_size, 0.739, tolerance = 0.001)
  })
})
```

---

## 11. Environment Reproducibility

### 11.1 R Version Lock

```
# .R-version
4.4.2
```

### 11.2 .Rprofile

```r
# .Rprofile
# Project-specific R configuration

# Set CRAN mirror to Posit Public Package Manager for reproducibility
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/2024-12-01"))

# Load renv for package management
source("renv/activate.R")

# Box module options
options(box.path = "R")

# Logging defaults
options(lgr.default_threshold = "info")

# Reproducibility
options(
  warnPartialMatchArgs = TRUE,
  warnPartialMatchAttr = TRUE,
  warnPartialMatchDollar = TRUE
)

# Development helpers (only in interactive sessions)
if (interactive()) {
  # Load devtools for development
  if (requireNamespace("devtools", quietly = TRUE)) {
    suppressMessages(require(devtools))
  }

  # Set editor
  options(editor = "vim")

  # Prompt with project name
  options(prompt = "[deadlift-study] > ")
}

cat("Loaded .Rprofile for deadlift-study project\n")
```

### 11.3 Dockerfile

```dockerfile
# docker/Dockerfile
FROM rocker/r-ver:4.4.2

# System dependencies for R packages
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# Install pak for fast package installation
RUN R -e "install.packages('pak', repos = 'https://r-lib.github.io/p/pak/stable/')"

# Install Quarto
RUN curl -LO https://quarto.org/download/latest/quarto-linux-amd64.deb \
    && dpkg -i quarto-linux-amd64.deb \
    && rm quarto-linux-amd64.deb

# Set working directory
WORKDIR /project

# Copy renv files first (for layer caching)
COPY renv.lock renv.lock
COPY renv/activate.R renv/activate.R
COPY .Rprofile .Rprofile

# Initialize renv and restore packages
RUN R -e "source('renv/activate.R'); renv::restore()"

# Copy project files
COPY . .

# Default command: run targets pipeline
CMD ["Rscript", "-e", "targets::tar_make()"]
```

### 11.4 docker-compose.yml

```yaml
# docker/docker-compose.yml
version: '3.8'

services:
  analysis:
    build:
      context: ..
      dockerfile: docker/Dockerfile
    volumes:
      - ../data:/project/data:ro
      - ../reports:/project/reports
    environment:
      - R_CONFIG_ACTIVE=production
    command: ["Rscript", "-e", "targets::tar_make()"]

  dev:
    build:
      context: ..
      dockerfile: docker/Dockerfile
    volumes:
      - ..:/project
    ports:
      - "8787:8787"
    environment:
      - R_CONFIG_ACTIVE=development
      - PASSWORD=rstudio
    command: ["/init"]
    # Uses rocker/rstudio features for development
```

---

## 12. Reporting & Documentation

### 12.1 Quarto Configuration

```yaml
# _quarto.yml
project:
  type: website
  output-dir: reports

website:
  title: "Deadlift Study: Meta-Analysis"
  navbar:
    left:
      - text: "Primary Analysis"
        file: analyses/01_primary_analysis.qmd
      - text: "Sensitivity"
        file: analyses/02_sensitivity_analysis.qmd
      - text: "Publication Bias"
        file: analyses/03_publication_bias.qmd
      - text: "Bayesian"
        file: analyses/04_bayesian_analysis.qmd

format:
  html:
    theme: cosmo
    toc: true
    toc-depth: 3
    code-fold: true
    code-summary: "Show code"
    fig-width: 8
    fig-height: 6

execute:
  freeze: auto
  cache: true

params:
  outcome: "strength"
  seed: 42
  confidence_level: 0.95
```

### 12.2 Parameterized Analysis Template

```qmd
---
title: "Primary Meta-Analysis: `{r} params$outcome`"
params:
  outcome: "strength"
  seed: 42
  config_file: "config/default.yml"
---

```{r setup}
#| include: false
box::use(
  R/utils/config[load_config],
  R/models/frequentist_model[FrequentistModel],
  targets[tar_read]
)

config <- load_config(params$config_file)
set.seed(params$seed)
```

## Methods

This analysis examines `{r} params$outcome` outcomes using multi-level meta-regression.

```{r model-results}
# Load pre-computed results from targets
model <- tar_read(paste0("model_", params$outcome))
summary <- model$summary()
```

## Results

### Overall Effect

The overall effect size was **`{r} round(summary$coefficients[1], 3)`**
(95% CI: `{r} round(summary$ci_lower[1], 3)` to `{r} round(summary$ci_upper[1], 3)`).

### Heterogeneity

- τ² = `{r} round(summary$tau2, 4)`
- I² = `{r} round(summary$i2, 1)`%
- Q = `{r} round(summary$q_statistic, 2)` (p = `{r} format.pval(summary$q_pval)`)
```

---

## 13. Decision Charts

### 13.1 When to Use Which Effect Size

```
                    ┌────────────────────────────────┐
                    │ What type of outcome measure?  │
                    └───────────────┬────────────────┘
                                    │
              ┌─────────────────────┼─────────────────────┐
              │                     │                     │
              ▼                     ▼                     ▼
    ┌─────────────────┐   ┌─────────────────┐   ┌─────────────────┐
    │ Continuous      │   │ Ratio/Percentage│   │ Dichotomous     │
    │ (strength, CSA) │   │ (% change)      │   │ (success/fail)  │
    └────────┬────────┘   └────────┬────────┘   └────────┬────────┘
             │                     │                     │
             ▼                     ▼                     ▼
    ┌─────────────────┐   ┌─────────────────┐   ┌─────────────────┐
    │ Use SMCR        │   │ Use ROMC        │   │ Use OR/RR       │
    │ (Hedge's g)     │   │ (log response   │   │ (odds/risk      │
    │                 │   │  ratio)         │   │  ratio)         │
    └─────────────────┘   └─────────────────┘   └─────────────────┘
```

### 13.2 Model Selection Workflow

```
                    ┌────────────────────────────────┐
                    │ Start: Choose analysis method  │
                    └───────────────┬────────────────┘
                                    │
                    ┌───────────────┼───────────────┐
                    │               │               │
                    ▼               ▼               ▼
           ┌──────────────┐ ┌──────────────┐ ┌──────────────┐
           │ Need prior   │ │ Standard     │ │ Publication  │
           │ information? │ │ frequentist? │ │ bias focus?  │
           └──────┬───────┘ └──────┬───────┘ └──────┬───────┘
                  │                │                │
                  ▼                ▼                ▼
           ┌──────────────┐ ┌──────────────┐ ┌──────────────┐
           │ Bayesian     │ │ rma.mv +     │ │ RoBMA /      │
           │ (brms)       │ │ robust VE    │ │ Selection    │
           └──────────────┘ └──────────────┘ │ Models       │
                                             └──────────────┘
```

### 13.3 Data Validation Decision Tree

```
                    ┌────────────────────────────────┐
                    │ Data Quality Check             │
                    └───────────────┬────────────────┘
                                    │
              ┌─────────────────────┼─────────────────────┐
              │                     │                     │
              ▼                     ▼                     ▼
    ┌─────────────────┐   ┌─────────────────┐   ┌─────────────────┐
    │ Missing values? │   │ Out of range?   │   │ Duplicates?     │
    └────────┬────────┘   └────────┬────────┘   └────────┬────────┘
             │                     │                     │
     ┌───────┼───────┐     ┌───────┼───────┐             │
     ▼       ▼       ▼     ▼       ▼       ▼             ▼
   ┌───┐   ┌───┐   ┌───┐ ┌───┐   ┌───┐   ┌───┐   ┌─────────────┐
   │<5%│   │5-20%│ │>20%│ │RIR│   │Year│  │SD │   │Remove exact │
   │   │   │    │  │   │  │>20│   │<1980│ │<=0│   │duplicates   │
   └─┬─┘   └─┬──┘  └─┬─┘ └─┬─┘   └──┬──┘ └─┬─┘   └─────────────┘
     │       │       │      │       │      │
     ▼       ▼       ▼      ▼       ▼      ▼
  ┌─────┐ ┌─────┐ ┌─────┐ ┌────────────────────┐
  │Impute││Review│ │Exclude│ │    FLAG & REVIEW   │
  │      ││each  │ │study  │ │                    │
  └─────┘└─────┘ └─────┘ └────────────────────┘
```

### 13.4 Testing Strategy Decision

```
                    ┌────────────────────────────────┐
                    │ What are you testing?          │
                    └───────────────┬────────────────┘
                                    │
        ┌───────────────────────────┼───────────────────────────┐
        │                           │                           │
        ▼                           ▼                           ▼
┌───────────────┐         ┌───────────────┐         ┌───────────────┐
│ Single unit   │         │ Statistical   │         │ Integration   │
│ of code       │         │ properties    │         │ behavior      │
└───────┬───────┘         └───────┬───────┘         └───────┬───────┘
        │                         │                         │
        ▼                         ▼                         ▼
┌───────────────┐         ┌───────────────┐         ┌───────────────┐
│ testthat      │         │ hedgehog +    │         │ testthat      │
│ unit tests    │         │ quickcheck    │         │ integration   │
│               │         │ property-based│         │ tests         │
│ - expect_*()  │         │               │         │               │
│ - Known I/O   │         │ - Invariants  │         │ - End-to-end  │
│ - Edge cases  │         │ - Boundaries  │         │ - Pipeline    │
└───────────────┘         │ - Fuzz testing│         │ - Multi-class │
                          └───────────────┘         └───────────────┘
```

---

## 14. Implementation Roadmap

### Phase 1: Foundation (Week 1)
- [ ] Create directory structure
- [ ] Set up .R-version, .Rprofile, .Renviron
- [ ] Initialize renv with core packages
- [ ] Configure pre-commit hooks
- [ ] Set up logging infrastructure

### Phase 2: Core Classes (Week 2)
- [ ] Implement BaseEntity R6 class
- [ ] Implement Study domain entity
- [ ] Implement SMCRCalculator
- [ ] Implement ROMCCalculator
- [ ] Write unit tests for calculators

### Phase 3: Models (Week 3)
- [ ] Implement FrequentistModel
- [ ] Implement BayesianModel (brms wrapper)
- [ ] Write property-based tests
- [ ] Validate against original script results

### Phase 4: Pipeline (Week 4)
- [ ] Set up targets pipeline
- [ ] Implement data loading targets
- [ ] Implement validation targets
- [ ] Implement model fitting targets

### Phase 5: Validation & Reporting (Week 5)
- [ ] Implement pointblank validation
- [ ] Create Quarto report templates
- [ ] Parameterize reports
- [ ] Generate validation reports

### Phase 6: Reproducibility (Week 6)
- [ ] Create Dockerfile
- [ ] Test container builds
- [ ] Set up DVC for data versioning
- [ ] Document setup process

---

## 15. Sources & References

### Package Documentation
- [R6 OOP](https://r6.r-lib.org/)
- [targets Pipeline](https://docs.ropensci.org/targets/)
- [renv Package Management](https://rstudio.github.io/renv/)
- [pak Package Installation](https://pak.r-lib.org/)
- [pointblank Data Validation](https://rstudio.github.io/pointblank/)
- [box Module System](https://klmr.me/box/)
- [lgr Logging](https://s-fleck.github.io/lgr/)
- [testthat Testing](https://testthat.r-lib.org/)
- [hedgehog Property Testing](https://cran.r-project.org/web/packages/hedgehog/)
- [precommit Hooks](https://lorenzwalthert.github.io/precommit/)

### R Version Management
- [rig R Installation Manager](https://github.com/r-lib/rig)

### Literate Programming
- [Quarto](https://quarto.org/)
- [Parameterized Reports](https://quarto.org/docs/computations/parameters.html)

### Data Tools
- [Arrow R Package](https://arrow.apache.org/docs/r/)
- [data.table](https://rdatatable.gitlab.io/data.table/)
- [polars](https://rpolars.r-universe.dev/)

### Containerization
- [Rocker Project](https://rocker-project.org/)
- [Rocker Reproducibility Essay (2024)](https://seanfobbe.com/posts/2024-12-16_reproducibility-limits-rocker-docker-images-r/)

### Meta-Analysis
- [metafor Package](https://www.metafor-project.org/)
- [brms Bayesian](https://paul-buerkner.github.io/brms/)
- [Doing Meta-Analysis in R](https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/)

### Research Compendia
- [rcompendium](https://frbcesab.github.io/rcompendium/)
- [The Turing Way - Research Compendia](https://book.the-turing-way.org/reproducible-research/compendia.html)

### Data Version Control
- [DVC](https://dvc.org/)
- [DVC R Integration Guide](https://dvc.org/doc/user-guide)

### Configuration Management
- [config Package](https://rstudio.github.io/config/)
- [Python Hydra (Inspiration)](https://hydra.cc/)

### Experiment Tracking
- [MLflow R API](https://mlflow.org/docs/latest/R-api.html)

---

*Document Version: 1.0.0*
*Last Updated: 2024-12-25*
*Authors: Research Team*
