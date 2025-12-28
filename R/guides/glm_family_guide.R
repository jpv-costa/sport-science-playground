# R/guides/glm_family_guide.R
# Service: GLM Family Selection Guide for Sport Science Research
#
# =============================================================================
# PRINCIPLES APPLIED (from CLAUDE.md)
# =============================================================================
#
# Based on "The Practitioner's Guide to Generalized Linear Models"
#
# NAMING PRINCIPLES:
# - Understandability: Domain terms (family, link, outcome_type)
# - Consistency: All public methods verb-based (suggest, describe, validate)
# - Distinguishability: suggest_family vs describe_family vs validate_choice
#
# SOLID PRINCIPLES:
# - SRP: Single responsibility - GLM family guidance only
# - OCP: Open for extension with new family types
# - DIP: Depends on data characteristics, not specific implementations
#
# SCIENTIFIC VALIDITY:
# - Family recommendations based on outcome variable characteristics
# - Link functions matched to scientific interpretation needs
# - References to statistical literature (Burnham & Anderson, Koo & Li)
#
# =============================================================================
# GLM FAMILY REFERENCE (from Practitioner's Guide)
# =============================================================================
# Outcome Type          | Family           | Link    | Example
# ----------------------|------------------|---------|-------------------------
# Binary                | Binomial         | logit   | Survived vs died
# Skewed Continuous     | Gamma            | inverse | Velocity at failure
# Count Data            | Poisson          | log     | Number of reps
# Count (overdispersed) | Negative Binomial| log     | Alternative to Poisson
# Ordered Discrete      | Binomial         | logit   | Likert scales
# Unordered Categorical | Multinomial      | logit   | Exercise type choice
# Zero-Inflated         | Two-part model   | varies  | Injury counts
# =============================================================================

box::use(
  R6[R6Class]
)

#' GLM Family Information
#'
#' Value object containing information about a GLM family.
#' @export
GlmFamilyInfo <- R6Class(
  classname = "GlmFamilyInfo",
  cloneable = FALSE,

  public = list(
    #' @field name Family name (e.g., "binomial", "gamma", "poisson")
    name = NULL,

    #' @field link Default link function
    link = NULL,

    #' @field outcome_type Description of suitable outcome type
    outcome_type = NULL,

    #' @field r_function R function to use (e.g., "glm", "glmer", "glmmTMB")
    r_function = NULL,

    #' @field interpretation How to interpret coefficients
    interpretation = NULL,

    #' @field assumptions Key assumptions to check
    assumptions = NULL,

    #' @field example Example use case
    example = NULL,

    #' @description Create GLM family info
    initialize = function(name, link, outcome_type, r_function,
                          interpretation, assumptions, example) {
      self$name <- name
      self$link <- link
      self$outcome_type <- outcome_type
      self$r_function <- r_function
      self$interpretation <- interpretation
      self$assumptions <- assumptions
      self$example <- example
    },

    #' @description Print family information
    print_info = function() {
      cat("GLM Family:", self$name, "\n")
      cat("  Link function:", self$link, "\n")
      cat("  Outcome type:", self$outcome_type, "\n")
      cat("  R function:", self$r_function, "\n")
      cat("  Interpretation:", self$interpretation, "\n")
      cat("  Assumptions:", paste(self$assumptions, collapse = "; "), "\n")
      cat("  Example:", self$example, "\n")
      invisible(self)
    },

    #' @description Convert to list
    to_list = function() {
      list(
        name = self$name,
        link = self$link,
        outcome_type = self$outcome_type,
        r_function = self$r_function,
        interpretation = self$interpretation,
        assumptions = self$assumptions,
        example = self$example
      )
    }
  )
)

#' GLM Family Guide
#'
#' R6 class providing guidance on selecting appropriate GLM families
#' for different outcome variable types in sport science research.
#'
#' @section Usage:
#' ```r
#' guide <- GlmFamilyGuide$new()
#' guide$suggest_family(outcome_type = "continuous_skewed")
#' guide$describe_family("gamma")
#' guide$print_decision_tree()
#' ```
#'
#' @export
GlmFamilyGuide <- R6Class(
  classname = "GlmFamilyGuide",
  cloneable = FALSE,

  public = list(

    #' @description Create GLM family guide
    initialize = function() {
      private$.init_families()
    },

    #' @description Suggest appropriate GLM family based on outcome characteristics
    #'
    #' @param outcome_type Character: "binary", "count", "continuous_normal",
    #'   "continuous_skewed", "ordinal", "nominal", "zero_inflated_count"
    #' @param mean_equals_variance Logical: For count data, is mean ~ variance?
    #' @return GlmFamilyInfo object with recommendation
    suggest_family = function(outcome_type, mean_equals_variance = NULL) {
      outcome_type <- tolower(outcome_type)

      recommendation <- switch(
        outcome_type,
        "binary" = private$.families$binomial,
        "count" = {
          if (is.null(mean_equals_variance) || mean_equals_variance) {
            private$.families$poisson
          } else {
            private$.families$negative_binomial
          }
        },
        "continuous_normal" = private$.families$gaussian,
        "continuous_skewed" = private$.families$gamma,
        "ordinal" = private$.families$ordinal,
        "nominal" = private$.families$multinomial,
        "zero_inflated_count" = private$.families$zero_inflated,
        stop("Unknown outcome type: ", outcome_type)
      )

      recommendation
    },

    #' @description Get detailed information about a specific family
    #' @param family_name Character: family name
    #' @return GlmFamilyInfo object
    describe_family = function(family_name) {
      family_name <- tolower(family_name)
      family_info <- private$.families[[family_name]]

      if (is.null(family_info)) {
        available <- paste(names(private$.families), collapse = ", ")
        stop("Unknown family: ", family_name, ". Available: ", available)
      }

      family_info
    },

    #' @description List all available families
    #' @return Character vector of family names
    list_families = function() {
      names(private$.families)
    },

    #' @description Print decision tree for family selection
    print_decision_tree = function() {
      cat("GLM Family Selection Decision Tree\n")
      cat("===================================\n\n")
      cat("What is your outcome variable type?\n\n")

      cat("1. BINARY (yes/no, survived/died)\n")
      cat("   -> Use: Binomial family with logit link (logistic regression)\n\n")

      cat("2. COUNT (integers >= 0)\n")
      cat("   Is mean approximately equal to variance?\n")
      cat("   -> YES: Use Poisson family with log link\n")
      cat("   -> NO:  Use Negative Binomial family with log link\n")
      cat("   -> Many zeros? Consider Zero-Inflated model\n\n")

      cat("3. CONTINUOUS\n")
      cat("   Is distribution approximately normal?\n")
      cat("   -> YES: Use Gaussian family with identity link (linear model)\n")
      cat("   -> NO (positively skewed, > 0): Use Gamma family\n\n")

      cat("4. ORDINAL (ordered categories, Likert scales)\n")
      cat("   -> Use: Ordinal logistic regression (cumulative link)\n\n")

      cat("5. NOMINAL (unordered categories)\n")
      cat("   -> Use: Multinomial logistic regression\n\n")

      cat("VBT-Specific Guidance:\n")
      cat("----------------------\n")
      cat("- Mean velocity (m/s): Gaussian or Gamma if skewed\n")
      cat("- RIR counts: Poisson or Negative Binomial\n")
      cat("- Set completion (yes/no): Binomial\n")
      cat("- RPE ratings (1-10): Ordinal logistic\n")

      invisible(self)
    },

    #' @description Analyze outcome variable and suggest family
    #'
    #' @param y Vector of outcome values
    #' @return List with suggested family and diagnostic info
    analyze_outcome = function(y) {
      stopifnot("y must be a vector" = is.vector(y) || is.factor(y))

      diagnostics <- list()

      # Determine type
      if (is.factor(y)) {
        n_levels <- nlevels(y)
        if (n_levels == 2) {
          diagnostics$detected_type <- "binary"
          diagnostics$suggestion <- self$suggest_family("binary")
        } else if (is.ordered(y)) {
          diagnostics$detected_type <- "ordinal"
          diagnostics$suggestion <- self$suggest_family("ordinal")
        } else {
          diagnostics$detected_type <- "nominal"
          diagnostics$suggestion <- self$suggest_family("nominal")
        }
      } else if (is.numeric(y)) {
        # Check if binary
        unique_vals <- unique(y[!is.na(y)])
        if (length(unique_vals) == 2 && all(unique_vals %in% c(0, 1))) {
          diagnostics$detected_type <- "binary"
          diagnostics$suggestion <- self$suggest_family("binary")
        }
        # Check if count (integers >= 0)
        else if (all(y >= 0, na.rm = TRUE) && all(y == floor(y), na.rm = TRUE)) {
          diagnostics$detected_type <- "count"
          diagnostics$mean <- mean(y, na.rm = TRUE)
          diagnostics$variance <- var(y, na.rm = TRUE)
          diagnostics$dispersion_ratio <- diagnostics$variance / diagnostics$mean
          diagnostics$zero_proportion <- mean(y == 0, na.rm = TRUE)

          if (diagnostics$zero_proportion > 0.5) {
            diagnostics$note <- "High proportion of zeros - consider zero-inflated model"
            diagnostics$suggestion <- self$suggest_family("zero_inflated_count")
          } else if (diagnostics$dispersion_ratio > 1.5) {
            diagnostics$note <- "Overdispersion detected (var/mean > 1.5)"
            diagnostics$suggestion <- self$suggest_family("count",
                                                          mean_equals_variance = FALSE)
          } else {
            diagnostics$suggestion <- self$suggest_family("count",
                                                          mean_equals_variance = TRUE)
          }
        }
        # Continuous
        else {
          diagnostics$mean <- mean(y, na.rm = TRUE)
          diagnostics$sd <- sd(y, na.rm = TRUE)
          diagnostics$min <- min(y, na.rm = TRUE)
          diagnostics$max <- max(y, na.rm = TRUE)
          diagnostics$skewness <- private$.calculate_skewness(y)

          if (diagnostics$min > 0 && abs(diagnostics$skewness) > 1) {
            diagnostics$detected_type <- "continuous_skewed"
            diagnostics$note <- paste0("Positively skewed (skewness = ",
                                       round(diagnostics$skewness, 2), ")")
            diagnostics$suggestion <- self$suggest_family("continuous_skewed")
          } else {
            diagnostics$detected_type <- "continuous_normal"
            diagnostics$suggestion <- self$suggest_family("continuous_normal")
          }
        }
      }

      diagnostics
    },

    #' @description Get coefficient interpretation guide
    #'
    #' @param family_name Character: family name
    #' @return Character string with interpretation guidance
    get_interpretation_guide = function(family_name) {
      family_name <- tolower(family_name)

      guides <- list(
        gaussian = paste(
          "GAUSSIAN (Identity Link):",
          "Coefficients represent the change in Y for a 1-unit change in X.",
          "Example: b = 0.04 means velocity increases by 0.04 m/s per RIR.",
          sep = "\n"
        ),

        binomial = paste(
          "BINOMIAL (Logit Link):",
          "Raw coefficients are log-odds.",
          "Exponentiate (exp(b)) to get odds ratios.",
          "Example: exp(0.5) = 1.65 means 65% higher odds per unit X increase.",
          sep = "\n"
        ),

        poisson = paste(
          "POISSON (Log Link):",
          "Raw coefficients are log(rate ratios).",
          "Exponentiate (exp(b)) to get multiplicative effect.",
          "Example: exp(-0.62) = 0.54 means 46% fewer counts per unit X increase.",
          sep = "\n"
        ),

        negative_binomial = paste(
          "NEGATIVE BINOMIAL (Log Link):",
          "Same interpretation as Poisson.",
          "Exponentiate (exp(b)) to get rate ratio.",
          "Accounts for overdispersion in count data.",
          sep = "\n"
        ),

        gamma = paste(
          "GAMMA (Inverse or Log Link):",
          "For log link: exponentiate for multiplicative effect.",
          "For inverse link: interpretation is more complex.",
          "Use predicted values for practical interpretation.",
          sep = "\n"
        )
      )

      guide <- guides[[family_name]]
      if (is.null(guide)) {
        return(paste("No interpretation guide available for:", family_name))
      }
      guide
    }
  ),

  private = list(
    .families = NULL,

    .init_families = function() {
      private$.families <- list(
        gaussian = GlmFamilyInfo$new(
          name = "Gaussian",
          link = "identity",
          outcome_type = "Continuous, approximately normally distributed",
          r_function = "glm(y ~ x, family = gaussian())",
          interpretation = "Coefficients are direct effects in original units",
          assumptions = c("Normality of residuals", "Homoscedasticity",
                          "Independence", "Linearity"),
          example = "Mean velocity (m/s) with symmetric distribution"
        ),

        binomial = GlmFamilyInfo$new(
          name = "Binomial",
          link = "logit",
          outcome_type = "Binary (0/1, yes/no, success/failure)",
          r_function = "glm(y ~ x, family = binomial())",
          interpretation = "Exponentiate coefficients for odds ratios",
          assumptions = c("Binary outcome", "Independence",
                          "No perfect separation"),
          example = "Set completed to failure (yes/no)"
        ),

        poisson = GlmFamilyInfo$new(
          name = "Poisson",
          link = "log",
          outcome_type = "Count data (integers >= 0)",
          r_function = "glm(y ~ x, family = poisson())",
          interpretation = "Exponentiate coefficients for rate ratios",
          assumptions = c("Mean equals variance", "Independence",
                          "Log-linear relationship"),
          example = "Number of repetitions completed"
        ),

        negative_binomial = GlmFamilyInfo$new(
          name = "Negative Binomial",
          link = "log",
          outcome_type = "Overdispersed count data (variance > mean)",
          r_function = "MASS::glm.nb(y ~ x)",
          interpretation = "Same as Poisson - exponentiate for rate ratios",
          assumptions = c("Overdispersion allowed", "Independence"),
          example = "Number of training sessions with high variance"
        ),

        gamma = GlmFamilyInfo$new(
          name = "Gamma",
          link = "inverse (or log)",
          outcome_type = "Positive continuous, positively skewed",
          r_function = "glm(y ~ x, family = Gamma(link = 'log'))",
          interpretation = "For log link: exponentiate for multiplicative effect",
          assumptions = c("Positive values only", "Right-skewed distribution"),
          example = "Minimum velocity threshold (MVT) at failure"
        ),

        ordinal = GlmFamilyInfo$new(
          name = "Ordinal (Cumulative Link)",
          link = "logit (cumulative)",
          outcome_type = "Ordered categorical (Likert scales, rankings)",
          r_function = "ordinal::clm(ordered_y ~ x)",
          interpretation = "Cumulative odds ratios for category transitions",
          assumptions = c("Proportional odds", "Ordinal outcome"),
          example = "RPE ratings (1-10 scale)"
        ),

        multinomial = GlmFamilyInfo$new(
          name = "Multinomial",
          link = "logit (baseline category)",
          outcome_type = "Unordered categorical (multiple categories)",
          r_function = "nnet::multinom(y ~ x)",
          interpretation = "Odds ratios relative to reference category",
          assumptions = c("Independence of irrelevant alternatives"),
          example = "Exercise type selection (squat/deadlift/bench)"
        ),

        zero_inflated = GlmFamilyInfo$new(
          name = "Zero-Inflated",
          link = "Two-part: logit + log",
          outcome_type = "Count data with excess zeros",
          r_function = "pscl::zeroinfl(y ~ x | z, dist = 'poisson')",
          interpretation = "Two sets of coefficients: presence vs. count",
          assumptions = c("Two generating processes: zeros and counts"),
          example = "Injury counts (many athletes with zero injuries)"
        )
      )
    },

    .calculate_skewness = function(x) {
      x <- x[!is.na(x)]
      n <- length(x)
      m <- mean(x)
      s <- sd(x)
      sum((x - m)^3) / (n * s^3)
    }
  )
)
