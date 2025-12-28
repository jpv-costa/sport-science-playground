# R/guides/statistical_methods_guide.R
# Educational Guide: Statistical Methods for Mixed Effects Models
#
# =============================================================================
# PRINCIPLES APPLIED (from CLAUDE.md)
# =============================================================================
#
# NAMING PRINCIPLES:
# - Understandability: Educational terms (explain_*, interpret_*, learn_*)
# - Consistency: All public methods follow verb_noun pattern
# - Distinguishability: Clear method names for different concepts
# - Conciseness: BF = Bayes Factor, LMM = Linear Mixed Model
#
# SOLID PRINCIPLES:
# - SRP: Single responsibility - only educational content generation
# - OCP: Open for extension via new educational topics
# - DIP: Depends on abstractions, not specific implementations
#
# CUPID PRINCIPLES:
# - Composable: Each explanation method is independent
# - Unix: Each method explains one concept
# - Predictable: Same inputs -> same explanations
# - Idiomatic: Follows R conventions
# - Domain-based: Names reflect statistical education concepts
# =============================================================================

box::use(
  R6[R6Class]
)

#' Statistical Methods Guide
#'
#' R6 class providing educational content about statistical methods
#' used in mixed effects modeling. Designed for practitioners who want
#' to understand the "why" behind the methods.
#'
#' @section Philosophy:
#' This guide follows the principle that understanding comes before
#' application. Each method explains not just "what" but "why" and "when".
#'
#' @export
StatisticalMethodsGuide <- R6Class(
  classname = "StatisticalMethodsGuide",
  cloneable = FALSE,

  public = list(
    #' @description Create a new guide instance
    initialize = function() {
      # No initialization needed
    },

    # =========================================================================
    # MIXED EFFECTS MODELS
    # =========================================================================

    #' @description Explain why mixed models are needed
    #' @return Character vector with explanation
    explain_why_mixed_models = function() {
      c(
        "WHY MIXED MODELS?",
        "================",
        "",
        "The Problem: Non-Independent Data",
        "---------------------------------",
        "When the same person provides multiple measurements, those measurements",
        "are correlated - they share that person's characteristics. Standard",
        "regression assumes each observation is independent, which is violated.",
        "",
        "The Consequence of Ignoring Clustering:",
        "1. Standard errors are too small (overconfident)",
        "2. P-values are too small (false positives)",
        "3. Confidence intervals are too narrow",
        "",
        "The Solution: Mixed Effects Models",
        "-----------------------------------",
        "Mixed models include RANDOM EFFECTS that account for:",
        "- Random intercepts: Each person has their own baseline",
        "- Random slopes: Each person has their own relationship",
        "",
        "Key Insight:",
        "Mixed models borrow strength across individuals while respecting",
        "that each person is unique. They're not 'pooled' (ignoring individuals)",
        "nor 'unpooled' (treating each person completely separately).",
        "",
        "When to Use:",
        "- Repeated measures on same subjects",
        "- Clustered data (students in schools, patients in clinics)",
        "- Longitudinal studies",
        "- Any hierarchical/nested data structure"
      )
    },

    #' @description Explain random effects structure decisions
    #' @return Character vector with explanation
    explain_random_effects_structure = function() {
      c(
        "CHOOSING RANDOM EFFECTS STRUCTURE",
        "==================================",
        "",
        "Three-Step Decision Framework:",
        "",
        "Step 1: Do individuals differ in their BASELINE?",
        "------------------------------------------------",
        "Question: At any given predictor value, do some individuals",
        "consistently score higher/lower than others?",
        "",
        "If YES -> Include RANDOM INTERCEPTS",
        "Each individual gets their own starting point.",
        "",
        "Step 2: Do individuals differ in their RESPONSE?",
        "------------------------------------------------",
        "Question: Does the effect of the predictor vary across individuals?",
        "Some might show strong effects, others weak or opposite effects.",
        "",
        "If YES -> Include RANDOM SLOPES",
        "Each individual gets their own relationship with the predictor.",
        "",
        "Step 3: Does the DESIGN support random slopes?",
        "----------------------------------------------",
        "Random slopes require:",
        "- Multiple observations per individual",
        "- Variation in the predictor within each individual",
        "- Enough individuals to estimate slope variance",
        "",
        "Rule of thumb: Need at least 5-10 observations per individual",
        "at different predictor values to estimate random slopes reliably.",
        "",
        "Common Mistake:",
        "Adding random slopes without enough data leads to convergence",
        "problems or unreliable estimates. Start simple, add complexity",
        "only when justified by data AND theory."
      )
    },

    # =========================================================================
    # BAYES FACTORS
    # =========================================================================

    #' @description Explain Bayes factors for model comparison
    #' @return Character vector with explanation
    explain_bayes_factors = function() {
      c(
        "BAYES FACTORS: QUANTIFYING EVIDENCE",
        "====================================",
        "",
        "The Intuition:",
        "--------------",
        "A Bayes factor (BF) answers: 'How many times more likely is the data",
        "under Model A than Model B?'",
        "",
        "BF = 10 means: 'The data is 10x more likely under Model A'",
        "BF = 0.1 means: 'The data is 10x more likely under Model B'",
        "",
        "Why Better Than P-Values for Model Comparison:",
        "-----------------------------------------------",
        "1. Directly compares models (not testing against null)",
        "2. Quantifies evidence strength (not just yes/no)",
        "3. Can support simpler models (not just reject)",
        "4. Intuitive interpretation (likelihood ratios)",
        "",
        "Jeffreys Scale for Interpretation:",
        "-----------------------------------",
        "BF > 100     : Decisive evidence",
        "BF 30-100    : Very strong evidence",
        "BF 10-30     : Strong evidence",
        "BF 3-10      : Moderate evidence",
        "BF 1-3       : Weak/anecdotal evidence",
        "BF ~ 1       : No evidence either way",
        "BF < 1       : Evidence for the OTHER model",
        "",
        "BIC Approximation:",
        "------------------",
        "For nested models, BF can be approximated from BIC:",
        "BF ≈ exp((BIC_simpler - BIC_complex) / 2)",
        "",
        "This approximation works well when:",
        "- Sample size is reasonably large (n > 30)",
        "- Models are properly specified",
        "- Data are not too sparse",
        "",
        "Practical Advice:",
        "-----------------",
        "- BF > 10: Strong enough for most decisions",
        "- BF 3-10: Consider other evidence too",
        "- BF < 3: Insufficient to distinguish models"
      )
    },

    # =========================================================================
    # VARIANCE COMPONENTS
    # =========================================================================

    #' @description Explain ICC and variance components
    #' @return Character vector with explanation
    explain_variance_components = function() {
      c(
        "VARIANCE COMPONENTS: UNDERSTANDING CLUSTERING",
        "==============================================",
        "",
        "Intraclass Correlation (ICC):",
        "-----------------------------",
        "ICC = Between-group variance / Total variance",
        "",
        "Interpretation:",
        "- ICC = 0.30 means 30% of variance is due to differences",
        "  BETWEEN individuals (stable individual differences)",
        "- The remaining 70% is WITHIN individuals (measurement-to-measurement)",
        "",
        "ICC Guidelines:",
        "- ICC < 0.05: Negligible clustering, maybe ignore",
        "- ICC 0.05-0.15: Small clustering, consider mixed models",
        "- ICC 0.15-0.25: Moderate clustering, use mixed models",
        "- ICC > 0.25: Substantial clustering, definitely use mixed models",
        "",
        "Design Effect:",
        "--------------",
        "Design effect = 1 + (average cluster size - 1) * ICC",
        "",
        "This tells you how much your effective sample size is reduced",
        "due to clustering. If design effect = 5, your 500 observations",
        "only provide as much information as 100 independent observations.",
        "",
        "Effective Sample Size:",
        "----------------------",
        "n_effective = n_total / design_effect",
        "",
        "Why This Matters:",
        "-----------------",
        "1. Affects power calculations: Need more participants, not just",
        "   more observations per participant",
        "2. Explains why adding more measurements per person has",
        "   diminishing returns",
        "3. Guides study design decisions",
        "",
        "Variance Decomposition Table:",
        "-----------------------------",
        "Report both between and within variance components.",
        "The ratio helps readers understand your data structure."
      )
    },

    # =========================================================================
    # MODEL DIAGNOSTICS
    # =========================================================================

    #' @description Explain visual-first diagnostics approach
    #' @return Character vector with explanation
    explain_visual_diagnostics = function() {
      c(
        "VISUAL-FIRST MODEL DIAGNOSTICS",
        "===============================",
        "",
        "Philosophy: Look First, Test Never",
        "-----------------------------------",
        "Statistical tests of assumptions (Shapiro-Wilk, Levene's test)",
        "are problematic because:",
        "",
        "1. With small samples: Tests have low power, miss violations",
        "2. With large samples: Tests flag trivial violations",
        "3. Binary outcome: Tests give yes/no, not 'how much'",
        "",
        "Better Approach: Visual Assessment",
        "-----------------------------------",
        "",
        "QQ Plot (Normality of Residuals):",
        "- Points should follow diagonal line",
        "- Slight deviations at tails are usually OK",
        "- Look for systematic curves (non-normality)",
        "- LMMs are robust to moderate non-normality with n > 30",
        "",
        "Residuals vs Fitted (Homoscedasticity):",
        "- Points should scatter randomly around zero",
        "- No fan/funnel shape (heteroscedasticity)",
        "- No curves (non-linearity)",
        "- Use cluster-robust SEs if concerned",
        "",
        "Random Effects QQ Plot:",
        "- Checks normality of random effects",
        "- Important for proper shrinkage",
        "- More sensitive with few clusters (< 30)",
        "",
        "What to Report:",
        "---------------",
        "1. Show the plots",
        "2. Describe what you see (not p-values)",
        "3. State whether violations are concerning",
        "4. Run robustness checks if borderline",
        "",
        "Key Insight:",
        "------------",
        "The question isn't 'Are assumptions perfectly met?'",
        "The question is 'Are violations severe enough to matter?'"
      )
    },

    # =========================================================================
    # ROBUSTNESS CHECKS
    # =========================================================================

    #' @description Explain robustness checks (pit stop philosophy)
    #' @return Character vector with explanation
    explain_robustness_checks = function() {
      c(
        "ROBUSTNESS CHECKS: THE PIT STOP PHILOSOPHY",
        "==========================================",
        "",
        "The Mindset Shift:",
        "------------------",
        "DON'T think: 'My model might be broken, let me check'",
        "DO think: 'Let me confirm my conclusions don't depend on",
        "          specific technical choices'",
        "",
        "Robustness checks are planned validation stops, not emergency repairs.",
        "",
        "Three Key Pit Stops:",
        "--------------------",
        "",
        "1. CLUSTER-ROBUST STANDARD ERRORS",
        "   What: SEs that don't assume homoscedasticity",
        "   Why: Handles unequal variance across clusters",
        "   Interpretation: If SE ratio (robust/Wald) ~ 1.0,",
        "   heteroscedasticity isn't a problem",
        "",
        "2. BOOTSTRAP CONFIDENCE INTERVALS",
        "   What: CIs from resampling, not formulas",
        "   Why: Makes no distributional assumptions",
        "   Interpretation: If bootstrap CIs match parametric CIs,",
        "   non-normality isn't a problem",
        "",
        "3. SENSITIVITY ANALYSIS",
        "   What: Re-run with different model specifications",
        "   Why: Check if conclusions depend on specific choices",
        "   Interpretation: If results consistent across specifications,",
        "   conclusions are robust to modeling choices",
        "",
        "How to Report:",
        "--------------",
        "'Standard and robust methods agree, confirming that our",
        "conclusions don't depend on specific technical assumptions.'",
        "",
        "When Results Disagree:",
        "----------------------",
        "1. Report both sets of results",
        "2. Investigate the source of disagreement",
        "3. Be more conservative in conclusions",
        "4. Consider it a finding, not a failure"
      )
    },

    # =========================================================================
    # EFFECT SIZES
    # =========================================================================

    #' @description Explain effect sizes for mixed models
    #' @return Character vector with explanation
    explain_effect_sizes = function() {
      c(
        "EFFECT SIZES: BEYOND STATISTICAL SIGNIFICANCE",
        "==============================================",
        "",
        "Why Effect Sizes Matter:",
        "------------------------",
        "P-values tell you IF an effect exists.",
        "Effect sizes tell you HOW BIG it is.",
        "",
        "A tiny effect can be 'significant' with enough data.",
        "Effect sizes let readers judge practical importance.",
        "",
        "Key Effect Sizes for Mixed Models:",
        "-----------------------------------",
        "",
        "1. COHEN'S d (Standardized Effect)",
        "   Formula: d = (Mean_high - Mean_low) / SD_pooled",
        "   ",
        "   Interpretation (Cohen's guidelines):",
        "   - d = 0.2: Small effect",
        "   - d = 0.5: Medium effect",
        "   - d = 0.8: Large effect",
        "   - d > 1.0: Very large effect",
        "",
        "2. R-SQUARED (Marginal vs Conditional)",
        "   ",
        "   R²_marginal: Variance explained by FIXED effects only",
        "   R²_conditional: Variance explained by FIXED + RANDOM effects",
        "   ",
        "   The gap between them shows how much individual differences",
        "   contribute beyond the main effects.",
        "",
        "3. PREDICTION DIFFERENCE",
        "   Raw-unit effect of a standardized change in predictor.",
        "   e.g., 'A 1-SD increase in X is associated with a 0.15 m/s",
        "   change in velocity'",
        "",
        "Practical Significance:",
        "-----------------------",
        "Always ask: 'Is this effect big enough to matter in practice?'",
        "",
        "Consider:",
        "- Measurement precision (can we detect a 0.01 m/s difference?)",
        "- Decision thresholds (does 0.05 m/s change training decisions?)",
        "- Cost-benefit (is the effect worth the intervention cost?)"
      )
    },

    # =========================================================================
    # CONVENIENCE METHODS
    # =========================================================================

    #' @description Get all topics as a list
    #' @return Named list of all explanations
    get_all_topics = function() {
      list(
        why_mixed_models = self$explain_why_mixed_models(),
        random_effects_structure = self$explain_random_effects_structure(),
        bayes_factors = self$explain_bayes_factors(),
        variance_components = self$explain_variance_components(),
        visual_diagnostics = self$explain_visual_diagnostics(),
        robustness_checks = self$explain_robustness_checks(),
        effect_sizes = self$explain_effect_sizes()
      )
    },

    #' @description Print a specific topic
    #' @param topic Character name of topic
    print_topic = function(topic) {
      topics <- self$get_all_topics()
      if (topic %in% names(topics)) {
        cat(paste(topics[[topic]], collapse = "\n"))
        cat("\n")
      } else {
        cat("Available topics:\n")
        cat(paste("-", names(topics), collapse = "\n"))
        cat("\n")
      }
    },

    #' @description Get quick reference card
    #' @return Character vector with quick reference
    get_quick_reference = function() {
      c(
        "STATISTICAL METHODS QUICK REFERENCE",
        "====================================",
        "",
        "WHEN TO USE MIXED MODELS:",
        "- Repeated measures / longitudinal data",
        "- Nested / hierarchical data",
        "- ICC > 0.05",
        "",
        "RANDOM EFFECTS DECISION:",
        "1. Baseline varies? -> Random intercepts",
        "2. Slopes vary? -> Random slopes",
        "3. Design supports? -> Include if yes",
        "",
        "BAYES FACTOR INTERPRETATION:",
        "BF > 10: Strong evidence",
        "BF 3-10: Moderate evidence",
        "BF 1-3: Weak evidence",
        "BF < 1: Evidence for other model",
        "",
        "ICC INTERPRETATION:",
        "< 0.05: Negligible",
        "0.05-0.15: Small",
        "0.15-0.25: Moderate",
        "> 0.25: Substantial",
        "",
        "EFFECT SIZE (Cohen's d):",
        "0.2: Small",
        "0.5: Medium",
        "0.8: Large",
        "",
        "DIAGNOSTICS APPROACH:",
        "1. Look at plots (not tests)",
        "2. Assess severity of violations",
        "3. Run robustness checks if concerned",
        "",
        "ROBUSTNESS PIT STOPS:",
        "1. Cluster-robust SEs (heteroscedasticity)",
        "2. Bootstrap CIs (non-normality)",
        "3. Sensitivity analysis (model choices)"
      )
    }
  )
)

#' Convenience function to print guide topic
#'
#' @param topic Character name of topic (or NULL for all)
#' @export
print_statistical_guide <- function(topic = NULL) {
  guide <- StatisticalMethodsGuide$new()
  if (is.null(topic)) {
    cat(paste(guide$get_quick_reference(), collapse = "\n"))
  } else {
    guide$print_topic(topic)
  }
}

#' Get available statistical guide topics
#'
#' @return Character vector of topic names
#' @export
get_guide_topics <- function() {
  c(
    "why_mixed_models",
    "random_effects_structure",
    "bayes_factors",
    "variance_components",
    "visual_diagnostics",
    "robustness_checks",
    "effect_sizes"
  )
}
