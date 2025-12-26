# R/calculators/correlation_imputer.R
# Service: Pre-Post Correlation Imputation
#
# SOLID Principles Applied:
# - SRP: Single responsibility - impute missing correlations (actor: statistician)
# - OCP: Open for extension via different imputation strategies
# - DIP: Depends on metafor abstractions for meta-analysis

box::use(
  R6[R6Class]
)

#' Correlation Imputation Result
#'
#' Value object containing the result of correlation imputation.
#'
#' @export
CorrelationImputationResult <- R6Class(


  classname = "CorrelationImputationResult",
  cloneable = FALSE,

  public = list(

    #' @description Create imputation result
    #' @param imputed_correlation The imputed r value
    #' @param fisher_z The Fisher's z transformation
    #' @param standard_error Standard error of the estimate
    #' @param n_studies Number of studies used for imputation
    #' @param method Imputation method used
    initialize = function(imputed_correlation,
                          fisher_z,
                          standard_error,
                          n_studies,
                          method = "meta-analytic") {

      private$.imputed_correlation <- imputed_correlation
      private$.fisher_z <- fisher_z
      private$.standard_error <- standard_error
      private$.n_studies <- n_studies
      private$.method <- method
    },

    #' @description Convert to list
    to_list = function() {
      list(
        imputed_correlation = private$.imputed_correlation,
        fisher_z = private$.fisher_z,
        standard_error = private$.standard_error,
        n_studies = private$.n_studies,
        method = private$.method
      )
    }
  ),

  active = list(
    imputed_correlation = function() private$.imputed_correlation,
    fisher_z = function() private$.fisher_z,
    standard_error = function() private$.standard_error,
    n_studies = function() private$.n_studies,
    method = function() private$.method
  ),

  private = list(
    .imputed_correlation = NULL,
    .fisher_z = NULL,
    .standard_error = NULL,
    .n_studies = NULL,
    .method = NULL
  )
)


#' Pre-Post Correlation Imputer
#'
#' Imputes missing pre-post correlations using meta-analysis of
#' available correlations within the dataset.
#'
#' @description
#' Service class for imputing missing pre-post correlations using
#' a multi-level meta-analysis approach. Follows Pelland et al. methodology.
#'
#' @export
CorrelationImputer <- R6Class(

  classname = "CorrelationImputer",

  public = list(

    # =========================================================================
    # Constructor
    # =========================================================================

    #' @description Create imputer
    #' @param random_effects Formula for random effects structure
    initialize = function(random_effects = "~ 1 | study/group/obs") {
      private$.random_effects <- random_effects
    },

    # =========================================================================
    # Public Methods
    # =========================================================================

    #' @description Calculate pre-post correlation from SDs
    #' @param pre_sd Pre-intervention standard deviation
    #' @param post_sd Post-intervention standard deviation
    #' @param delta_sd Change score standard deviation
    #' @return Numeric correlation or NA if invalid
    calculate_correlation = function(pre_sd, post_sd, delta_sd) {
      # Formula: ri = (sd_pre^2 + sd_post^2 - sd_delta^2) / (2 * sd_pre * sd_post)
      numerator <- pre_sd^2 + post_sd^2 - delta_sd^2
      denominator <- 2 * pre_sd * post_sd

      if (denominator == 0) return(NA_real_)

      r <- numerator / denominator

      # Correlation must be between -1 and 1
      if (is.na(r) || r < -1 || r > 1) {
        return(NA_real_)
      }

      r
    },

    #' @description Impute correlation via meta-analysis
    #' @param data Data frame with ri (correlations), n (sample sizes),
    #'             study, group, obs columns
    #' @param use_robust Use robust variance estimation
    #' @return CorrelationImputationResult object
    impute_via_meta_analysis = function(data, use_robust = TRUE) {
      if (!requireNamespace("metafor", quietly = TRUE)) {
        stop("Package 'metafor' is required for meta-analytic imputation")
      }
      if (!requireNamespace("psych", quietly = TRUE)) {
        stop("Package 'psych' is required for Fisher z transformation")
      }

      # Filter to rows with valid correlations
      valid_data <- data[!is.na(data$ri), ]
      n_studies <- length(unique(valid_data$study))

      if (nrow(valid_data) < 2) {
        stop("At least 2 valid correlations required for meta-analysis")
      }

      # Convert to Fisher's z
      escalc_data <- metafor::escalc(
        measure = "ZCOR",
        ri = ri,
        ni = n,
        data = valid_data
      )

      # Fit multi-level model
      model <- metafor::rma.mv(
        yi = yi,
        V = vi,
        data = escalc_data,
        random = list(~ 1 | study, ~ 1 | group, ~ 1 | obs),
        method = "REML",
        test = "t"
      )

      # Apply robust variance estimation if requested
      if (use_robust) {
        model <- metafor::robust(model, cluster = escalc_data$study)
      }

      # Convert back to r
      fisher_z <- model$b[1, 1]
      imputed_r <- psych::fisherz2r(fisher_z)

      CorrelationImputationResult$new(
        imputed_correlation = imputed_r,
        fisher_z = fisher_z,
        standard_error = model$se[1],
        n_studies = n_studies,
        method = ifelse(use_robust, "robust-meta-analysis", "meta-analysis")
      )
    },

    #' @description Apply imputation to fill missing correlations
    #' @param data Data frame with ri column
    #' @param imputed_value Value to use for imputation
    #' @return Data frame with filled correlations
    apply_imputation = function(data, imputed_value) {
      data$ri <- ifelse(is.na(data$ri), imputed_value, data$ri)
      data
    }
  ),

  private = list(
    .random_effects = NULL
  )
)
