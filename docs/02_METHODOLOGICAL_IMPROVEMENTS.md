# Methodological Improvements for the PTF Meta-Analysis

> **Document Version:** 1.0
> **Purpose:** Advanced statistical methods and enhancements to strengthen the meta-analysis
> **Last Updated:** December 2025

---

## Table of Contents

1. [Publication Bias Assessment](#1-publication-bias-assessment)
2. [Assumption Testing & Diagnostics](#2-assumption-testing--diagnostics)
3. [Heterogeneity Assessment Improvements](#3-heterogeneity-assessment-improvements)
4. [Bayesian Meta-Analysis](#4-bayesian-meta-analysis)
5. [Causal Inference with DAGs](#5-causal-inference-with-dags)
6. [Cross-Validation & Model Validation](#6-cross-validation--model-validation)
7. [Dose-Response Modeling Enhancements](#7-dose-response-modeling-enhancements)
8. [Sensitivity Analyses](#8-sensitivity-analyses)
9. [Individual Participant Data Meta-Analysis](#9-individual-participant-data-ipd-meta-analysis)
10. [Robust Variance Estimation Enhancements](#10-robust-variance-estimation-enhancements)
11. [Structural Equation Modeling](#11-structural-equation-modeling-sem--path-analysis)
12. [Conformal Prediction Methods](#12-conformal-prediction-methods)
13. [Implementation Priority & Roadmap](#13-implementation-priority--roadmap)
14. [Additional R Packages Required](#14-additional-r-packages-required)
15. [References](#15-references)

---

## 1. Publication Bias Assessment

### 1.1 Why This Matters

**Current Gap:** The existing script does not assess publication bias.

**Evidence of Problem:**
- A 2024 BMJ review found **55% of meta-analyses in high-impact journals don't assess publication bias**
- P-hacking and selective reporting are common in exercise science
- Smaller studies with non-significant results are less likely to be published

### 1.2 Methods to Implement

#### 1.2.1 Funnel Plot Visualization

```r
# Basic funnel plot
funnel(primary.str.g.best.model,
       main = "Funnel Plot: Strength Outcomes",
       xlab = "Hedge's g")

# Enhanced funnel with contours
funnel(primary.str.g.best.model,
       level = c(90, 95, 99),
       shade = c("white", "gray75", "gray60"),
       refline = 0,
       legend = TRUE)
```

**Interpretation:**
- Symmetric funnel → No evidence of publication bias
- Asymmetric funnel → Potential bias (small studies missing on one side)

#### 1.2.2 Egger's Regression Test

```r
# Egger's test for funnel plot asymmetry
regtest(primary.str.g.best.model, model = "lm")

# For multilevel models, use modified approach
regtest(primary.str.g.best.model,
        predictor = "sei",    # Standard error as predictor
        model = "rma")
```

**Interpretation:**
- p < 0.10 suggests significant asymmetry
- Intercept significantly different from 0 indicates bias

#### 1.2.3 PET-PEESE Method

**Precision-Effect Test (PET):**
```r
# PET: regress effect on SE
pet_model <- rma.mv(yi, vi,
                    mods = ~ sqrt(vi),  # SE as moderator
                    random = list(~1|study, ~1|group, ~1|obs),
                    data = primary.data.str)

# If PET intercept is significant, use PEESE
peese_model <- rma.mv(yi, vi,
                      mods = ~ vi,  # Variance as moderator
                      random = list(~1|study, ~1|group, ~1|obs),
                      data = primary.data.str)
```

**Decision Rule:**
1. If PET p > 0.10: Use unadjusted estimate
2. If PET p ≤ 0.10 and PET estimate significant: Use PEESE estimate
3. If PET p ≤ 0.10 and PET estimate non-significant: Use PET estimate

#### 1.2.4 Selection Models

```r
# Three-parameter selection model (3PSM)
selmodel(primary.str.g.best.model,
         type = "stepfun",
         steps = c(0.025, 0.5))  # Two-tailed significance thresholds

# Weight function selection model
selmodel(primary.str.g.best.model,
         type = "beta")
```

**What it does:** Models probability of publication as a function of p-value

#### 1.2.5 Trim and Fill

```r
# Trim and fill method
trimfill(primary.str.g.best.model, side = "left")
trimfill(primary.str.g.best.model, side = "right")

# Visualize
funnel(trimfill(primary.str.g.best.model))
```

**Limitation:** Assumes bias is the only source of asymmetry

#### 1.2.6 Robust Bayesian Meta-Analysis (RoBMA)

```r
library(RoBMA)

# Fit RoBMA model
robma_model <- RoBMA(
  d = primary.data.str$yi,
  se = sqrt(primary.data.str$vi),
  study_names = primary.data.str$author.year,
  seed = 42,
  parallel = TRUE
)

summary(robma_model)
```

**Advantages:**
- Model-averages over publication bias models
- Produces posterior probabilities for effect and heterogeneity
- Handles both publication bias and p-hacking

### 1.3 Reporting Template

```markdown
## Publication Bias Assessment

### Visual Assessment
- Funnel plot [Figure X] shows [symmetric/asymmetric] distribution

### Statistical Tests
| Test | Statistic | p-value | Interpretation |
|------|-----------|---------|----------------|
| Egger's | z = X.XX | p = X.XX | [No evidence/Evidence] of asymmetry |
| PET | b0 = X.XX | p = X.XX | [Adjusted/Unadjusted] estimate appropriate |
| Trim & Fill | k imputed = X | - | [X] studies potentially missing |

### Sensitivity Analysis
| Method | Original Estimate | Adjusted Estimate | Change |
|--------|------------------|-------------------|--------|
| PET-PEESE | X.XX [CI] | X.XX [CI] | X% |
| Selection Model | X.XX [CI] | X.XX [CI] | X% |
| RoBMA | X.XX [CrI] | - | - |
```

---

## 2. Assumption Testing & Diagnostics

### 2.1 Why This Matters

**Current Gap:** No formal assumption testing in the script.

**Evidence of Problem:**
- A 2023 BMC Medicine study found **15-26% of Cochrane meta-analyses violate normality assumptions**
- Violations don't bias point estimates but **invalidate confidence intervals**

### 2.2 Normality of Random Effects

#### 2.2.1 Shapiro-Wilk Test

```r
# Extract standardized residuals
resid_std <- rstandard(primary.str.g.best.model)

# Shapiro-Wilk test
shapiro.test(resid_std$z)

# Interpretation
# p > 0.05: No evidence against normality
# p < 0.05: Evidence of non-normality
```

#### 2.2.2 Q-Q Plot

```r
# Q-Q plot of standardized residuals
qqnorm(resid_std$z, main = "Q-Q Plot: Standardized Residuals")
qqline(resid_std$z, col = "red", lwd = 2)

# Add confidence envelope
library(car)
qqPlot(resid_std$z,
       main = "Q-Q Plot with 95% Confidence Envelope",
       ylab = "Standardized Residuals")
```

#### 2.2.3 Density Plot

```r
# Density plot of residuals vs. normal
ggplot(data.frame(resid = resid_std$z), aes(x = resid)) +
  geom_histogram(aes(y = ..density..), bins = 30,
                 fill = "steelblue", alpha = 0.7) +
  stat_function(fun = dnorm, color = "red", size = 1) +
  labs(title = "Distribution of Standardized Residuals",
       x = "Standardized Residual", y = "Density") +
  theme_minimal()
```

### 2.3 Homoscedasticity

```r
# Residuals vs. fitted values plot
fitted_vals <- fitted(primary.str.g.best.model)

ggplot(data.frame(fitted = fitted_vals, resid = resid_std$z),
       aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  labs(title = "Residuals vs. Fitted Values",
       x = "Fitted Values", y = "Standardized Residuals") +
  theme_minimal()

# Breusch-Pagan test (approximate for meta-analysis)
# Check if residual variance relates to predictors
bp_test <- lm(resid_std$z^2 ~ primary.data.str$avg.rir +
                              primary.data.str$load.set)
summary(bp_test)
```

### 2.4 Influence Diagnostics

#### 2.4.1 Leave-One-Out Analysis

```r
# Leave-one-out diagnostics
loo_results <- leave1out(primary.str.g.best.model)

# Visualize
ggplot(loo_results, aes(x = reorder(slab, estimate), y = estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), width = 0.2) +
  geom_hline(yintercept = coef(primary.str.g.best.model)[1],
             linetype = "dashed", color = "red") +
  coord_flip() +
  labs(title = "Leave-One-Out Analysis",
       x = "Study Removed", y = "Pooled Estimate") +
  theme_minimal()
```

#### 2.4.2 Influence Statistics

```r
# Calculate influence statistics
inf <- influence(primary.str.g.best.model)

# Plot all influence diagnostics
plot(inf, layout = c(4, 2))

# Key statistics:
# - Cook's distance: Overall influence
# - DFFITS: Influence on fitted values
# - Covariance ratio: Influence on precision
# - Hat values: Leverage
```

#### 2.4.3 Outlier Detection

```r
# Identify outliers using studentized residuals
outliers <- which(abs(rstudent(primary.str.g.best.model)$z) > 2)

# Iterative outlier detection (Meng et al., 2024)
find_outliers <- function(model, threshold = 2) {
  resid <- rstudent(model)$z
  outliers <- which(abs(resid) > threshold)
  return(outliers)
}

# Remove outliers and refit
primary.data.str.clean <- primary.data.str[-outliers, ]
```

### 2.5 Multicollinearity Check

```r
# Variance Inflation Factors for moderators
library(car)

# Create linear model version for VIF
lm_version <- lm(yi ~ avg.rir + load.set + weeks +
                      as.numeric(train.status) + as.numeric(set.rep.equated),
                 data = primary.data.str)

vif(lm_version)

# VIF > 5 indicates concerning multicollinearity
# VIF > 10 indicates severe multicollinearity
```

### 2.6 Reporting Template

```markdown
## Model Diagnostics

### Normality Assessment
| Test | Statistic | p-value | Conclusion |
|------|-----------|---------|------------|
| Shapiro-Wilk | W = X.XX | p = X.XX | [Normal/Non-normal] |

Visual inspection of Q-Q plot [Figure X] suggests [approximate normality/departures from normality].

### Influential Studies
[X] studies identified as potentially influential based on Cook's distance > [threshold].

| Study | Cook's D | Effect on Estimate | Recommendation |
|-------|----------|-------------------|----------------|
| Author (Year) | X.XX | ±X.XX | [Retain/Sensitivity] |

### Multicollinearity
| Predictor | VIF | Conclusion |
|-----------|-----|------------|
| avg.rir | X.XX | [OK/Concerning] |
| load.set | X.XX | [OK/Concerning] |
```

---

## 3. Heterogeneity Assessment Improvements

### 3.1 Prediction Intervals (Already Implemented - Enhancement)

**Current Status:** The script calculates prediction intervals.

**Enhancement:** Ensure proper reporting and interpretation.

```r
# Calculate prediction interval
pi_lower <- estimate - qt(0.975, df) * sqrt(se^2 + sum(tau^2))
pi_upper <- estimate + qt(0.975, df) * sqrt(se^2 + sum(tau^2))

# Create forest plot with PI
forest(primary.str.g.best.model,
       addpred = TRUE,           # Add prediction interval
       header = TRUE,
       xlim = c(-2, 4),
       alim = c(-1, 3))
```

**Key Interpretation Points:**
- **72%** of significant meta-analyses have PIs that cross null
- If PI includes clinically important effects in both directions, practical recommendations are limited

### 3.2 Tau² Decomposition

```r
# Extract variance components
var_components <- data.frame(
  Level = c("Study", "Group", "Observation", "Total"),
  Tau2 = c(primary.str.g.best.model$sigma2[1],
           primary.str.g.best.model$sigma2[2],
           primary.str.g.best.model$sigma2[3],
           sum(primary.str.g.best.model$sigma2)),
  Proportion = c(primary.str.g.best.model$sigma2[1] / sum(primary.str.g.best.model$sigma2),
                 primary.str.g.best.model$sigma2[2] / sum(primary.str.g.best.model$sigma2),
                 primary.str.g.best.model$sigma2[3] / sum(primary.str.g.best.model$sigma2),
                 1.00)
)

print(var_components)
```

### 3.3 I² Variants for Multilevel Models

```r
# Multilevel I² calculation
i2_ml <- function(model) {
  W <- diag(1/model$vi)
  X <- model.matrix(model)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

  # Total variance
  total_var <- sum(model$sigma2) + (model$k - model$p) / sum(diag(P))

  # I² for each level
  i2_study <- model$sigma2[1] / total_var * 100
  i2_group <- model$sigma2[2] / total_var * 100
  i2_obs <- model$sigma2[3] / total_var * 100
  i2_total <- sum(model$sigma2) / total_var * 100

  return(list(study = i2_study, group = i2_group,
              obs = i2_obs, total = i2_total))
}

# Use orchaRd package for I² in multilevel models
library(orchaRd)
i2_results <- i2_ml(primary.str.g.best.model)
```

### 3.4 Outlier Handling Strategies

```r
# Strategy 1: Winsorization
winsorize_effects <- function(data, threshold = 2) {
  z <- (data$yi - mean(data$yi)) / sd(data$yi)
  data$yi_winsor <- ifelse(abs(z) > threshold,
                           sign(z) * threshold * sd(data$yi) + mean(data$yi),
                           data$yi)
  return(data)
}

# Strategy 2: Robust estimation
library(robumeta)
robust_model <- robu(yi ~ avg.rir + load.set,
                     data = primary.data.str,
                     studynum = study,
                     var.eff.size = vi)

# Strategy 3: Sensitivity analysis excluding outliers
outlier_ids <- which(abs(rstudent(primary.str.g.best.model)$z) > 3)
sensitivity_model <- update(primary.str.g.best.model,
                            subset = -outlier_ids)
```

---

## 4. Bayesian Meta-Analysis

### 4.1 Why Switch to Bayesian?

| Advantage | Frequentist | Bayesian |
|-----------|-------------|----------|
| **Output** | Point estimate + CI | Full posterior distribution |
| **Few studies** | Poor τ² estimation | Better τ² via priors |
| **Interpretation** | "95% of CIs contain true value" | "95% probability true value is in CrI" |
| **Prior knowledge** | Cannot incorporate | Formally incorporated |
| **Model comparison** | AIC/BIC | WAIC/LOO-CV with SE |
| **Probability statements** | Not directly | P(effect > threshold) |

### 4.2 Implementation with brms

#### 4.2.1 Basic Bayesian Meta-Regression

```r
library(brms)

# Prepare data
primary.data.str$sei <- sqrt(primary.data.str$vi)

# Define priors
priors <- c(
  prior(normal(0, 1), class = Intercept),      # Weakly informative
  prior(normal(0, 0.5), class = b),            # Skeptical prior for effects
  prior(cauchy(0, 0.5), class = sd)            # Half-Cauchy for tau
)

# Fit Bayesian model
bayes_str_g <- brm(
  yi | se(sei) ~ log1p(avg.rir) + load.set + set.rep.equated + weeks + train.status +
    (1 | study) + (1 | study:group) + (1 | obs),
  data = primary.data.str,
  prior = priors,
  iter = 4000,
  warmup = 1000,
  chains = 4,
  cores = 4,
  seed = 42,
  control = list(adapt_delta = 0.95)
)

# Model summary
summary(bayes_str_g)
```

#### 4.2.2 Prior Sensitivity Analysis

```r
# Skeptical prior (effect likely near zero)
prior_skeptical <- c(
  prior(normal(0, 0.5), class = Intercept),
  prior(normal(0, 0.2), class = b),
  prior(cauchy(0, 0.3), class = sd)
)

# Enthusiastic prior (effect likely positive)
prior_enthusiastic <- c(
  prior(normal(0.5, 1), class = Intercept),
  prior(normal(0, 0.5), class = b),
  prior(cauchy(0, 0.5), class = sd)
)

# Fit with different priors
bayes_skeptical <- update(bayes_str_g, prior = prior_skeptical)
bayes_enthusiastic <- update(bayes_str_g, prior = prior_enthusiastic)

# Compare posteriors
posterior_samples <- as.data.frame(rbind(
  cbind(Model = "Weakly Informative", as_draws_df(bayes_str_g)),
  cbind(Model = "Skeptical", as_draws_df(bayes_skeptical)),
  cbind(Model = "Enthusiastic", as_draws_df(bayes_enthusiastic))
))
```

#### 4.2.3 Posterior Inference

```r
# Posterior summaries
posterior_summary(bayes_str_g)

# Probability of effect > 0
hypothesis(bayes_str_g, "log1pavg.rir < 0")

# Probability of meaningful effect (e.g., > 0.2 SMD per RIR)
hypothesis(bayes_str_g, "log1pavg.rir < -0.2")

# Region of Practical Equivalence (ROPE)
library(bayestestR)
rope(bayes_str_g, range = c(-0.1, 0.1), ci = 0.95)

# Bayes Factors for parameters
bayesfactor_parameters(bayes_str_g)
```

#### 4.2.4 Model Comparison

```r
# Leave-one-out cross-validation
loo_bayes <- loo(bayes_str_g)
print(loo_bayes)

# Compare models
loo_compare(loo(bayes_linear), loo(bayes_log), loo(bayes_spline))

# WAIC
waic(bayes_str_g)
```

#### 4.2.5 Posterior Predictive Checks

```r
# Posterior predictive check
pp_check(bayes_str_g, type = "dens_overlay", ndraws = 100)
pp_check(bayes_str_g, type = "stat", stat = "mean")
pp_check(bayes_str_g, type = "scatter_avg")
```

### 4.3 Reporting Bayesian Results

```markdown
## Bayesian Meta-Regression Results

### Prior Specification
| Parameter | Prior | Justification |
|-----------|-------|---------------|
| Intercept | Normal(0, 1) | Weakly informative, centered at null |
| β coefficients | Normal(0, 0.5) | Skeptical, small effects expected |
| τ (between-study SD) | Half-Cauchy(0, 0.5) | Weakly informative, allows large values |

### Posterior Summaries
| Parameter | Mean | 95% CrI | P(effect > 0) | ROPE |
|-----------|------|---------|---------------|------|
| log1p(RIR) | X.XX | [X.XX, X.XX] | X.X% | X% in ROPE |

### Model Fit
- WAIC: X.XX (SE = X.XX)
- LOO-IC: X.XX (SE = X.XX)
- Posterior predictive checks indicate [adequate/inadequate] model fit
```

---

## 5. Causal Inference with DAGs

### 5.1 Why Use DAGs?

**Benefits for Meta-Analysis:**
1. **Clarify confounders** vs. mediators vs. colliders
2. **Guide covariate selection** — "more covariates" isn't always better
3. **Identify selection bias** sources
4. **Transparently communicate assumptions**
5. **Compare adjustment strategies** across included studies

### 5.2 Creating a DAG for PTF Meta-Analysis

#### 5.2.1 Conceptual DAG

```r
library(ggdag)
library(dagitty)

# Define the DAG
ptf_dag <- dagify(
  # Outcome relationships
  Strength ~ RIR + Load + Volume + Weeks + TrainStatus + Genetics + Nutrition,
  Hypertrophy ~ RIR + Load + Volume + Weeks + TrainStatus + Genetics + Nutrition,

  # RIR relationships (what influences RIR assignment)
  RIR ~ StudyDesign + Progression + FailureDef,

  # Load relationships
  Load ~ StudyDesign + TrainStatus,

  # Volume relationships
  Volume ~ StudyDesign + Load,

  # Publication bias (selection)
  Published ~ EffectSize + SampleSize + Journal,

  # Unmeasured confounders
  TrainStatus ~ Age + Genetics,

  # Define exposure and outcome

  exposure = "RIR",
  outcome = "Strength",

  # Define latent/unmeasured variables
  latent = c("Genetics", "Nutrition", "Journal")
)

# Visualize
ggdag(ptf_dag, layout = "sugiyama") +
  theme_dag() +
  labs(title = "Directed Acyclic Graph for PTF Meta-Analysis")
```

#### 5.2.2 Identify Adjustment Sets

```r
# What should we adjust for?
adjustmentSets(ptf_dag, exposure = "RIR", outcome = "Strength",
               type = "minimal")

# What are the backdoor paths?
paths(ptf_dag, from = "RIR", to = "Strength")

# Identify colliders (should NOT adjust for)
# Volume may be a collider if RIR → Volume and Load → Volume
```

#### 5.2.3 Using DAGitty Online Tool

1. Go to [https://www.dagitty.net/](https://www.dagitty.net/)
2. Draw nodes: RIR, Strength, Load, Volume, Weeks, TrainStatus, etc.
3. Draw arrows representing causal relationships
4. Mark exposure (RIR) and outcome (Strength/Hypertrophy)
5. Request minimal adjustment set

### 5.3 DAG-Informed Analysis

```r
# Based on DAG, minimum sufficient adjustment set might be:
# {TrainStatus, StudyDesign}

# Compare with current model (may be over-adjusted)
dag_informed_model <- rma.mv(yi, vi,
                              mods = ~ log1p(avg.rir) + train.status + within.between.design,
                              random = list(~1|study, ~1|group, ~1|obs),
                              data = primary.data.str)

# Check if results differ from over-adjusted model
# If Load is a mediator (RIR → Load → Strength), adjusting for it
# blocks part of the causal effect
```

### 5.4 Reporting DAG Analysis

```markdown
## Causal Assumptions

### Directed Acyclic Graph
[Figure X] presents our assumed causal structure for the RIR-Strength relationship.

### Key Assumptions
1. RIR directly affects strength outcomes
2. Training status confounds the RIR-Strength relationship
3. Load may be a **mediator** (RIR → Load → Strength), not a confounder
4. Volume may be a **collider** (RIR → Volume ← Load)

### Implications for Analysis
| Current Adjustment | DAG Recommendation | Implication |
|--------------------|-------------------|-------------|
| Load | Potential mediator | May over-adjust |
| Volume | Potential collider | Should not adjust |
| Training status | Confounder | Must adjust |

### Sensitivity to Unmeasured Confounding
Using E-value analysis, an unmeasured confounder would need an RR of X.XX with both RIR and Strength to explain away the observed association.
```

---

## 6. Cross-Validation & Model Validation

### 6.1 Leave-One-Out Cross-Validation

```r
# Leave-one-out at study level
studies <- unique(primary.data.str$study)

loo_cv <- sapply(studies, function(s) {
  # Fit model excluding study s
  train_data <- primary.data.str[primary.data.str$study != s, ]
  test_data <- primary.data.str[primary.data.str$study == s, ]

  model <- rma.mv(yi, vi,
                  mods = ~ log1p(avg.rir) + load.set + set.rep.equated +
                           weeks + train.status,
                  random = list(~1|study, ~1|group, ~1|obs),
                  data = train_data, method = "REML")

  # Predict for held-out study
  newmods <- model.matrix(~ log1p(avg.rir) + load.set + set.rep.equated +
                           weeks + train.status, data = test_data)[, -1]
  pred <- predict(model, newmods = newmods)

  # Calculate squared error
  se <- (pred$pred - test_data$yi)^2
  return(mean(se))
})

# Cross-validation MSE
cv_mse <- mean(loo_cv)
cv_rmse <- sqrt(cv_mse)

cat("LOO-CV RMSE:", round(cv_rmse, 3), "\n")
```

### 6.2 K-Fold Cross-Validation

```r
# 5-fold cross-validation by study
set.seed(42)
studies <- unique(primary.data.str$study)
folds <- sample(rep(1:5, length.out = length(studies)))
names(folds) <- studies

kfold_cv <- sapply(1:5, function(k) {
  test_studies <- studies[folds == k]
  train_data <- primary.data.str[!primary.data.str$study %in% test_studies, ]
  test_data <- primary.data.str[primary.data.str$study %in% test_studies, ]

  # Skip if test set is empty
  if(nrow(test_data) == 0) return(NA)

  model <- rma.mv(yi, vi,
                  mods = ~ log1p(avg.rir) + load.set + set.rep.equated +
                           weeks + train.status,
                  random = list(~1|study, ~1|group, ~1|obs),
                  data = train_data, method = "REML")

  newmods <- model.matrix(~ log1p(avg.rir) + load.set + set.rep.equated +
                           weeks + train.status, data = test_data)[, -1]
  pred <- predict(model, newmods = newmods)

  mse <- mean((pred$pred - test_data$yi)^2)
  return(mse)
})

# Report
cv_results <- data.frame(
  Fold = 1:5,
  MSE = kfold_cv,
  RMSE = sqrt(kfold_cv)
)
cv_results$Overall <- c(rep("", 4), paste("Mean RMSE:", round(mean(sqrt(kfold_cv), na.rm = TRUE), 3)))
print(cv_results)
```

### 6.3 Statistical Validity (Vn Statistic)

```r
# Willis (2017) cross-validation statistic
calculate_vn <- function(model, data) {
  n <- nrow(data)
  pred_errors <- numeric(n)

  for(i in 1:n) {
    # Fit model without observation i
    model_i <- update(model, subset = -i)

    # Predict observation i
    newmods <- model.matrix(model)[i, -1]
    pred <- predict(model_i, newmods = matrix(newmods, nrow = 1))

    pred_errors[i] <- (data$yi[i] - pred$pred)^2 / pred$se^2
  }

  # Vn statistic
  Vn <- mean(pred_errors)

  # Under null (model is valid), Vn ~ 1
  # Vn >> 1 suggests poor generalization

  return(list(Vn = Vn, pred_errors = pred_errors))
}
```

### 6.4 Model Comparison via CV

```r
# Compare functional forms using CV
models <- list(
  linear = ~ avg.rir + load.set + weeks + train.status,
  log = ~ log1p(avg.rir) + load.set + weeks + train.status,
  quadratic = ~ poly(avg.rir, 2) + load.set + weeks + train.status,
  spline = ~ bs(spline.rir, degree = 1, knots = 0) + load.set + weeks + train.status
)

cv_comparison <- sapply(models, function(formula) {
  model <- rma.mv(yi, vi, mods = formula,
                  random = list(~1|study, ~1|group, ~1|obs),
                  data = primary.data.str)

  # 5-fold CV
  # ... (implementation as above)
  # Return RMSE
})

cv_comparison_df <- data.frame(
  Model = names(models),
  CV_RMSE = cv_comparison
) %>% arrange(CV_RMSE)
```

---

## 7. Dose-Response Modeling Enhancements

### 7.1 Fractional Polynomials

```r
library(mfp)

# First-degree fractional polynomial
fp1_model <- rma.mv(yi, vi,
                    mods = ~ I(avg.rir^0.5) + load.set + weeks + train.status,
                    random = list(~1|study, ~1|group, ~1|obs),
                    data = primary.data.str)

# Second-degree fractional polynomial (FP2)
# Common power pairs: (-2,-1), (-1,-0.5), (-0.5,0), (0,0.5), (0.5,1), (1,2)
fp2_powers <- list(
  c(-2, -1), c(-1, -0.5), c(-0.5, 0), c(0, 0.5), c(0.5, 1), c(1, 2), c(2, 3)
)

fp2_results <- lapply(fp2_powers, function(p) {
  # Handle log transformation when power = 0
  x1 <- if(p[1] == 0) log(avg.rir + 1) else (avg.rir + 1)^p[1]
  x2 <- if(p[2] == 0) log(avg.rir + 1) else (avg.rir + 1)^p[2]

  model <- rma.mv(yi, vi,
                  mods = ~ x1 + x2 + load.set + weeks + train.status,
                  random = list(~1|study, ~1|group, ~1|obs),
                  data = primary.data.str)

  list(powers = p, AIC = AIC(model), model = model)
})

# Select best FP2 by AIC
best_fp2 <- fp2_results[[which.min(sapply(fp2_results, function(x) x$AIC))]]
```

### 7.2 Gaussian Process Regression

```r
library(brms)

# Gaussian process for smooth non-parametric dose-response
gp_model <- brm(
  yi | se(sei) ~ gp(avg.rir, k = 10, c = 5/4) + load.set + weeks + train.status +
    (1 | study) + (1 | study:group),
  data = primary.data.str %>% mutate(sei = sqrt(vi)),
  prior = c(
    prior(normal(0, 1), class = Intercept),
    prior(normal(0, 0.5), class = b),
    prior(inv_gamma(2, 1), class = lscale),
    prior(cauchy(0, 1), class = sdgp)
  ),
  iter = 4000, warmup = 1000, chains = 4,
  control = list(adapt_delta = 0.95)
)

# Plot smooth dose-response
conditional_effects(gp_model, effects = "avg.rir")
```

### 7.3 Penalized Splines (P-Splines)

```r
library(mgcv)

# For visualization and comparison
# Note: mgcv doesn't directly handle meta-analysis, but concepts apply

# Using metafor with many knots (approximates penalized spline)
ns_model <- rma.mv(yi, vi,
                   mods = ~ ns(avg.rir, df = 5) + load.set + weeks + train.status,
                   random = list(~1|study, ~1|group, ~1|obs),
                   data = primary.data.str)
```

### 7.4 Model Comparison for Dose-Response

```r
# Compare all dose-response models
dose_response_models <- list(
  linear = rma.mv(yi, vi, mods = ~ avg.rir + load.set + weeks + train.status, ...),
  log = rma.mv(yi, vi, mods = ~ log1p(avg.rir) + load.set + weeks + train.status, ...),
  quadratic = rma.mv(yi, vi, mods = ~ poly(avg.rir, 2) + load.set + weeks + train.status, ...),
  cubic = rma.mv(yi, vi, mods = ~ poly(avg.rir, 3) + load.set + weeks + train.status, ...),
  spline_linear = rma.mv(yi, vi, mods = ~ bs(spline.rir, degree = 1, knots = 0) + ..., ...),
  spline_cubic = rma.mv(yi, vi, mods = ~ bs(spline.rir, degree = 3, knots = 0) + ..., ...),
  rcs = rma.mv(yi, vi, mods = ~ rcs(avg.rir, 4) + load.set + weeks + train.status, ...),
  fp2_best = best_fp2$model
)

# Comparison table
comparison_table <- data.frame(
  Model = names(dose_response_models),
  AIC = sapply(dose_response_models, AIC),
  BIC = sapply(dose_response_models, BIC),
  logLik = sapply(dose_response_models, logLik)
) %>% arrange(AIC)
```

---

## 8. Sensitivity Analyses

### 8.1 Comprehensive Sensitivity Analysis Framework

| Analysis | Purpose | Implementation |
|----------|---------|----------------|
| Exclude high RoB studies | Risk of bias | Filter by quality |
| Exclude small studies | Small-study effects | Filter n < threshold |
| Alternative correlation | ri imputation | Re-run with ri = 0.5, 0.7, 0.9 |
| Different effect sizes | Metric robustness | Compare SMCR vs SMD vs ROM |
| One-step vs two-step | Pooling approach | Compare estimates |
| Trim outliers | Outlier influence | Exclude |z| > 3 |
| Cluster at sample level | Overlapping samples | Adjust clustering |
| Different knot locations | Spline specification | Vary knots |
| Exclude specific subgroups | Subgroup influence | Sequential exclusion |

### 8.2 Correlation Sensitivity

```r
# Test sensitivity to imputed correlations
correlations <- c(0.3, 0.5, 0.7, 0.9)

corr_sensitivity <- lapply(correlations, function(r) {
  # Recalculate effect sizes with different r
  data_temp <- primary.data.str
  data_temp$ri_test <- r
  data_temp$delta.sd_test <- sqrt(data_temp$pre.sd^2 + data_temp$post.sd^2 -
                                   2 * r * data_temp$pre.sd * data_temp$post.sd)

  # Recalculate SMCR
  data_temp <- escalc(measure = "SMCR",
                      m1i = post.mean, m2i = pre.mean,
                      sd1i = pre.sd, ni = n, ri = ri_test,
                      data = data_temp,
                      var.names = c("yi_test", "vi_test"))

  # Refit model
  model <- rma.mv(yi_test, vi_test,
                  mods = ~ log1p(avg.rir) + load.set + set.rep.equated +
                           weeks + train.status,
                  random = list(~1|study, ~1|group, ~1|obs),
                  data = data_temp)

  return(list(r = r, estimate = coef(model)[2],
              ci.lb = model$ci.lb[2], ci.ub = model$ci.ub[2]))
})

# Plot sensitivity
corr_sensitivity_df <- do.call(rbind, lapply(corr_sensitivity, as.data.frame))

ggplot(corr_sensitivity_df, aes(x = factor(r), y = estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Sensitivity to Pre-Post Correlation Assumption",
       x = "Assumed Correlation (r)", y = "Effect of log(RIR+1)") +
  theme_minimal()
```

### 8.3 Effect Size Metric Sensitivity

```r
# Compare SMCR vs ROM
smcr_model <- primary.str.g.best.model
rom_model <- primary.str.rom.best.model

# Also try raw mean difference
rmd_data <- escalc(measure = "MC", m1i = post.mean, m2i = pre.mean,
                   sd1i = delta.sd, ni = n, ri = ri, data = primary.data.str)

rmd_model <- rma.mv(yi, vi,
                    mods = ~ log1p(avg.rir) + load.set + set.rep.equated +
                             weeks + train.status,
                    random = list(~1|study, ~1|group, ~1|obs),
                    data = rmd_data)

# Comparison
effect_size_comparison <- data.frame(
  Metric = c("SMCR (Hedge's g)", "ROMC (log RR)", "Raw Mean Change"),
  Estimate = c(coef(smcr_model)[2], coef(rom_model)[2], coef(rmd_model)[2]),
  SE = c(smcr_model$se[2], rom_model$se[2], rmd_model$se[2]),
  p = c(smcr_model$pval[2], rom_model$pval[2], rmd_model$pval[2])
)
```

### 8.4 Outlier Sensitivity

```r
# Identify potential outliers
resid <- rstudent(primary.str.g.best.model)
outliers <- which(abs(resid$z) > 2)

# Fit model excluding outliers
sensitivity_no_outliers <- update(primary.str.g.best.model,
                                   subset = -outliers)

# Compare
cat("Original estimate:", round(coef(primary.str.g.best.model)[2], 4), "\n")
cat("Excluding outliers:", round(coef(sensitivity_no_outliers)[2], 4), "\n")
cat("Change:", round((coef(sensitivity_no_outliers)[2] - coef(primary.str.g.best.model)[2]) /
                      coef(primary.str.g.best.model)[2] * 100, 1), "%\n")
```

### 8.5 Training Status Sensitivity

```r
# Trained only
trained_only <- rma.mv(yi, vi,
                       mods = ~ log1p(avg.rir) + load.set + set.rep.equated + weeks,
                       random = list(~1|study, ~1|group, ~1|obs),
                       data = primary.data.str[primary.data.str$train.status == "trained", ])

# Untrained only
untrained_only <- rma.mv(yi, vi,
                         mods = ~ log1p(avg.rir) + load.set + set.rep.equated + weeks,
                         random = list(~1|study, ~1|group, ~1|obs),
                         data = primary.data.str[primary.data.str$train.status == "untrained", ])
```

### 8.6 Reporting Sensitivity Analyses

```markdown
## Sensitivity Analyses

### Summary Table
| Analysis | Original Estimate | Sensitivity Estimate | Δ% | Conclusion |
|----------|-------------------|---------------------|-----|------------|
| Exclude high RoB | X.XX [CI] | X.XX [CI] | X% | Robust |
| r = 0.5 | X.XX [CI] | X.XX [CI] | X% | Sensitive |
| r = 0.9 | X.XX [CI] | X.XX [CI] | X% | Robust |
| Exclude outliers | X.XX [CI] | X.XX [CI] | X% | Robust |
| Trained only | X.XX [CI] | X.XX [CI] | X% | Consistent |
| RIR ≤ 10 | X.XX [CI] | X.XX [CI] | X% | Robust |

### Interpretation
Results were [robust/sensitive] to analytical choices. The direction and
magnitude of the RIR-Strength relationship remained [consistent/variable]
across sensitivity analyses, with estimates ranging from X.XX to X.XX.
```

---

## 9. Individual Participant Data (IPD) Meta-Analysis

### 9.1 Why IPD is the "Gold Standard"

| Advantage | Aggregate Data | IPD |
|-----------|----------------|-----|
| **Effect modifiers** | Study-level only (ecological fallacy) | Individual-level |
| **Subgroup analyses** | Low power | High power |
| **Missing data** | Limited options | Multiple imputation |
| **Standardization** | Post-hoc | A priori |
| **Bias assessment** | Study-level | Individual-level |
| **Covariate adjustment** | Aggregate | Precise |

### 9.2 Data Requirements for IPD

```r
# Required individual-level variables
ipd_schema <- data.frame(
  Variable = c("participant_id", "study_id", "age", "sex", "training_experience",
               "pre_strength", "post_strength", "pre_muscle_size", "post_muscle_size",
               "avg_rir", "load_percent", "sets_per_week", "frequency",
               "intervention_weeks"),
  Type = c("ID", "ID", "continuous", "binary", "continuous",
           "continuous", "continuous", "continuous", "continuous",
           "continuous", "continuous", "continuous", "continuous",
           "continuous"),
  Description = c("Unique participant identifier", "Study identifier",
                  "Age in years", "Male=1, Female=0", "Years of RT experience",
                  "Baseline 1RM (kg)", "Post-intervention 1RM (kg)",
                  "Baseline muscle CSA/thickness", "Post muscle CSA/thickness",
                  "Average RIR across training", "Average load (% 1RM)",
                  "Weekly set volume", "Sessions per week",
                  "Duration of intervention")
)
```

### 9.3 One-Stage IPD Analysis

```r
library(lme4)
library(lmerTest)

# One-stage mixed model (if IPD available)
ipd_model <- lmer(
  post_strength ~ pre_strength + log1p(avg_rir) + load_percent +
                  sets_per_week + frequency + intervention_weeks +
                  age + sex + training_experience +
                  (1 + log1p(avg_rir) | study_id) + (1 | participant_id),
  data = ipd_data,
  REML = TRUE
)

summary(ipd_model)
```

### 9.4 Two-Stage IPD Analysis

```r
# Stage 1: Estimate within-study effects
stage1_results <- ipd_data %>%
  group_by(study_id) %>%
  do({
    model <- lm(post_strength ~ pre_strength + log1p(avg_rir) +
                load_percent + age + sex, data = .)
    data.frame(
      estimate = coef(model)["log1p(avg_rir)"],
      se = summary(model)$coefficients["log1p(avg_rir)", "Std. Error"],
      n = nrow(.)
    )
  })

# Stage 2: Meta-analyze study effects
stage2_model <- rma(yi = estimate, sei = se, data = stage1_results,
                    method = "REML")
```

### 9.5 Practical Considerations

**Challenges:**
1. **Time-consuming data acquisition** - requires author contact and data sharing
2. **Resource intensive** - data cleaning, harmonization, validation
3. **Incomplete IPD** - rarely available from all studies
4. **Variable heterogeneity** - different measures across studies

**Recommendations:**
1. Start with collaborative consortium (e.g., prospective IPD-MA)
2. Use established data sharing platforms
3. Pre-register IPD-MA protocol
4. Plan for missing covariates

---

## 10. Robust Variance Estimation Enhancements

### 10.1 Current Implementation

```r
# Current approach in script
primary.str.g.best.model <- robust.rma.mv(
  primary.str.g.best.model.prep,
  cluster = primary.data.str$study
)
```

### 10.2 Small-Sample Corrections

```r
library(clubSandwich)

# CR2 estimator with Satterthwaite degrees of freedom
robust_cr2 <- coef_test(primary.str.g.best.model.prep,
                        vcov = "CR2",
                        cluster = primary.data.str$study,
                        test = "Satterthwaite")

print(robust_cr2)
```

### 10.3 Wild Bootstrap for Few Clusters

```r
# When number of clusters (studies) is small (< 20)
library(clubSandwich)

# Wild cluster bootstrap
boot_results <- Wald_test(primary.str.g.best.model.prep,
                          constraints = "log1p(avg.rir)",
                          vcov = "CR2",
                          cluster = primary.data.str$study,
                          test = "naive-tp")

# Bootstrap confidence intervals
wild_bootstrap_ci <- function(model, cluster_var, n_boot = 1000) {
  # Implementation of wild cluster bootstrap
  # ... (complex implementation)
}
```

### 10.4 Comparison of Robust Estimators

```r
# Compare different robust estimators
robust_comparison <- data.frame(
  Estimator = c("Naive", "CR0", "CR1", "CR2"),
  Estimate = rep(coef(primary.str.g.best.model.prep)[2], 4),
  SE = c(
    primary.str.g.best.model.prep$se[2],
    sqrt(vcovCR(primary.str.g.best.model.prep, cluster = primary.data.str$study,
                type = "CR0")[2,2]),
    sqrt(vcovCR(primary.str.g.best.model.prep, cluster = primary.data.str$study,
                type = "CR1")[2,2]),
    sqrt(vcovCR(primary.str.g.best.model.prep, cluster = primary.data.str$study,
                type = "CR2")[2,2])
  )
)

robust_comparison$t <- robust_comparison$Estimate / robust_comparison$SE
print(robust_comparison)
```

---

## 11. Structural Equation Modeling (SEM) / Path Analysis

### 11.1 When to Use SEM in Meta-Analysis

- **Mediation hypotheses** (e.g., RIR → Load → Strength)
- **Latent variable modeling** (e.g., "training quality" as latent)
- **Complex path models** with multiple mediators
- **Meta-analytic structural equation modeling (MASEM)**

### 11.2 Testing Mediation

```r
library(lavaan)

# If Load mediates RIR → Strength relationship
# (Would require correlation matrix from each study)

mediation_model <- '
  # Direct effects
  strength ~ c*rir + b*load + weeks + train_status

  # Mediation path
  load ~ a*rir

  # Indirect effect
  indirect := a*b
  total := c + a*b
'

# For meta-analytic SEM, use metaSEM package
library(metaSEM)

# Would require correlation matrices from each study
# masem_model <- tssem1(Cov = correlation_matrices, n = sample_sizes)
```

### 11.3 Latent Variable Approach

```r
# Conceptual: "Training Quality" as latent construct
latent_model <- '
  # Latent variable
  training_quality =~ load + volume + frequency + rir

  # Structural model
  strength ~ training_quality + weeks + train_status
'
```

---

## 12. Conformal Prediction Methods

### 12.1 Why Conformal Prediction for Meta-Analysis?

**Current Gap:** Prediction intervals in meta-analysis rely on parametric assumptions (normal random effects) that are frequently violated.

**Evidence of Problem:**
- A 2023 BMC Medicine study found **15-26% of Cochrane meta-analyses violate the normality assumption**
- With few studies, parametric prediction intervals have **poor finite-sample coverage**
- The Higgins-Thompson-Spiegelhalter prediction interval can have coverage well below nominal

**What Conformal Prediction Offers:**
- **Distribution-free guarantees** — no normality assumption required
- **Finite-sample validity** — works even with few studies
- **Model-agnostic** — works with any underlying predictor
- **Explicit coverage guarantee** — P(true value ∈ PI) ≥ 1 - α

### 12.2 Key Concepts

#### 12.2.1 The Core Idea

Conformal prediction constructs prediction sets by:
1. Computing "nonconformity scores" on calibration data (how unusual is each observation?)
2. Using the distribution of these scores to calibrate prediction sets for new observations
3. Guaranteeing that P(Y_new ∈ prediction set) ≥ 1 - α under exchangeability

> "Conformal prediction is a user-friendly paradigm for creating statistically rigorous uncertainty sets/intervals for the predictions of such models. Critically, the sets are valid in a distribution-free sense: they possess explicit, non-asymptotic guarantees even without distributional assumptions."
> — Angelopoulos & Bates (2022)

#### 12.2.2 Comparison with Parametric Prediction Intervals

| Aspect | Parametric PI (Current) | Conformal PI |
|--------|------------------------|--------------|
| **Formula** | est ± t × √(SE² + τ²) | est ± q̂(nonconformity scores) |
| **Assumptions** | Normal random effects, known τ² | Only exchangeability |
| **Coverage guarantee** | Asymptotic (approximate) | Finite-sample (exact) |
| **With model misspecification** | Coverage degrades | Still valid |
| **With few studies** | Poor coverage | Valid coverage |
| **Computational cost** | Low | Medium (split) to High (full) |

#### 12.2.3 Types of Conformal Prediction

| Method | Description | Trade-off |
|--------|-------------|-----------|
| **Full Conformal** | Refit model for each candidate y | Most accurate, computationally expensive |
| **Split Conformal** | Split data into train/calibration | Fast, slight efficiency loss |
| **Jackknife+** | Leave-one-out with + correction | Good balance, no data splitting |
| **CV+** | Cross-validation based | Uses all data efficiently |
| **CQR** | Conformalized Quantile Regression | Adaptive interval width |

### 12.3 The Hierarchical Challenge

Standard conformal prediction assumes **exchangeability** of observations. In meta-analysis with structure `~1|study/group/obs`, observations within studies are **not exchangeable** with observations from other studies.

**Solution:** Hierarchical Conformal Prediction (Dunn et al., 2022)

> "The validity of conformal prediction hinges on the exchangeability of the data, which does not hold when groups of observations come from distinct distributions... We extend conformal methods to this hierarchical setting."

#### 12.3.1 Approaches for Hierarchical Data

| Method | Description | Reference |
|--------|-------------|-----------|
| **CDF Pooling** | Pool empirical CDFs across groups | Dunn et al. (2022) |
| **Single Subsampling** | Subsample one observation per group | Dunn et al. (2022) |
| **Repeated Subsampling** | Average over multiple subsamples | Dunn et al. (2022) |
| **Multi-scale CP** | Separate conformity at each level | Gasparin et al. (2025) |

### 12.4 Implementation

#### 12.4.1 Basic Split Conformal for Meta-Regression

```r
#' Split Conformal Prediction for Meta-Analysis
#'
#' @param data Full dataset
#' @param formula Model formula
#' @param alpha Miscoverage rate (default 0.05 for 95% PI)
#' @param cal_prop Proportion of data for calibration (default 0.2)
#'
#' @return List with model, quantile, and prediction function

conformal_meta <- function(data, formula, alpha = 0.05, cal_prop = 0.2) {

  # 1. Split data into training and calibration sets
  set.seed(42)
  n <- nrow(data)
  cal_idx <- sample(1:n, size = floor(cal_prop * n))
  train_data <- data[-cal_idx, ]
  cal_data <- data[cal_idx, ]

  # 2. Fit model on training set
  model <- rma.mv(yi, vi,
                  mods = formula,
                  random = list(~1|study, ~1|group, ~1|obs),
                  data = train_data,
                  method = "REML")

  # 3. Get predictions on calibration set
  X_cal <- model.matrix(formula, data = cal_data)[, -1]
  cal_preds <- predict(model, newmods = X_cal)

  # 4. Compute nonconformity scores (absolute residuals)
  #    Can also use: |residual| / SE for studentized version
  nonconf_scores <- abs(cal_data$yi - cal_preds$pred)

  # 5. Compute quantile for prediction interval
  #    Use ceiling((n_cal + 1)(1 - alpha)) / n_cal for finite-sample correction
  n_cal <- length(nonconf_scores)
  q_level <- ceiling((n_cal + 1) * (1 - alpha)) / n_cal
  q_hat <- quantile(nonconf_scores, probs = min(q_level, 1))

  # Return prediction function
  predict_conformal <- function(newdata) {
    X_new <- model.matrix(formula, data = newdata)[, -1]
    preds <- predict(model, newmods = X_new)

    data.frame(
      pred = preds$pred,
      ci.lb = preds$ci.lb,
      ci.ub = preds$ci.ub,
      pi.lb.conformal = preds$pred - q_hat,
      pi.ub.conformal = preds$pred + q_hat
    )
  }

  return(list(
    model = model,
    q_hat = q_hat,
    nonconf_scores = nonconf_scores,
    n_calibration = n_cal,
    predict = predict_conformal
  ))
}

# Usage
cp_result <- conformal_meta(
  data = primary.data.str,
  formula = ~ log1p(avg.rir) + load.set + set.rep.equated + weeks + train.status,
  alpha = 0.05
)

# Predict for new data
new_predictions <- cp_result$predict(new_data)
```

#### 12.4.2 Conformalized Quantile Regression (CQR)

CQR produces **adaptive** prediction intervals (wider where uncertainty is higher):

```r
#' Conformalized Quantile Regression for Meta-Analysis
#'
#' Uses quantile regression to get initial interval, then calibrates with CP

library(quantreg)

conformal_quantile_meta <- function(data, formula, alpha = 0.05, cal_prop = 0.2) {

  # Split data
  set.seed(42)
  n <- nrow(data)
  cal_idx <- sample(1:n, size = floor(cal_prop * n))
  train_data <- data[-cal_idx, ]
  cal_data <- data[cal_idx, ]

  # Fit quantile regression for lower and upper bounds
  # Using weighted quantile regression (weight = 1/vi)
  formula_full <- as.formula(paste("yi", paste(as.character(formula), collapse = "")))

  model_lo <- rq(formula_full, tau = alpha/2, weights = 1/vi, data = train_data)
  model_hi <- rq(formula_full, tau = 1 - alpha/2, weights = 1/vi, data = train_data)

  # Get initial quantile predictions on calibration set
  q_lo_cal <- predict(model_lo, newdata = cal_data)
  q_hi_cal <- predict(model_hi, newdata = cal_data)

  # Compute CQR nonconformity scores
  # E_i = max(q_lo - Y_i, Y_i - q_hi)
  nonconf_scores <- pmax(q_lo_cal - cal_data$yi, cal_data$yi - q_hi_cal)

  # Compute quantile
  n_cal <- length(nonconf_scores)
  q_level <- ceiling((n_cal + 1) * (1 - alpha)) / n_cal
  q_hat <- quantile(nonconf_scores, probs = min(q_level, 1))

  # Prediction function
  predict_cqr <- function(newdata) {
    q_lo_new <- predict(model_lo, newdata = newdata)
    q_hi_new <- predict(model_hi, newdata = newdata)

    data.frame(
      pred = (q_lo_new + q_hi_new) / 2,
      pi.lb.cqr = q_lo_new - q_hat,
      pi.ub.cqr = q_hi_new + q_hat
    )
  }

  return(list(
    model_lo = model_lo,
    model_hi = model_hi,
    q_hat = q_hat,
    predict = predict_cqr
  ))
}
```

#### 12.4.3 Hierarchical Conformal Prediction (Study-Level)

For valid coverage accounting for study clustering:

```r
#' Hierarchical Conformal Prediction via Subsampling
#'
#' Addresses non-exchangeability by subsampling one observation per study

hierarchical_conformal_meta <- function(data, formula, alpha = 0.05,
                                         n_subsample = 100) {

  studies <- unique(data$study)
  n_studies <- length(studies)

  # Repeated subsampling approach
  q_hats <- numeric(n_subsample)

  for (b in 1:n_subsample) {
    # Subsample one observation per study
    subsample_idx <- sapply(studies, function(s) {
      study_rows <- which(data$study == s)
      sample(study_rows, 1)
    })

    subsample_data <- data[subsample_idx, ]

    # Split into train/calibration (by study)
    n_cal_studies <- floor(0.2 * n_studies)
    cal_studies <- sample(studies, n_cal_studies)

    train_data <- subsample_data[!subsample_data$study %in% cal_studies, ]
    cal_data <- subsample_data[subsample_data$study %in% cal_studies, ]

    # Fit model
    model <- tryCatch({
      rma.mv(yi, vi,
             mods = formula,
             random = list(~1|study, ~1|obs),
             data = train_data,
             method = "REML")
    }, error = function(e) NULL)

    if (is.null(model)) next

    # Compute nonconformity scores
    X_cal <- model.matrix(formula, data = cal_data)[, -1]
    cal_preds <- predict(model, newmods = X_cal)
    nonconf_scores <- abs(cal_data$yi - cal_preds$pred)

    # Store quantile
    n_cal <- length(nonconf_scores)
    q_level <- ceiling((n_cal + 1) * (1 - alpha)) / n_cal
    q_hats[b] <- quantile(nonconf_scores, probs = min(q_level, 1), na.rm = TRUE)
  }

  # Average quantile across subsamples
  q_hat_avg <- mean(q_hats, na.rm = TRUE)

  return(list(
    q_hat = q_hat_avg,
    q_hat_distribution = q_hats,
    n_subsamples = n_subsample
  ))
}
```

#### 12.4.4 Using tidymodels Framework

```r
library(tidymodels)
library(probably)

# Prepare data for tidymodels
meta_recipe <- recipe(yi ~ avg.rir + load.set + set.rep.equated + weeks + train.status,
                      data = primary.data.str) %>%
  step_log(avg.rir, offset = 1) %>%
  step_dummy(all_nominal_predictors())

# Linear regression model (simplified; for full meta-analysis use custom model)
lm_spec <- linear_reg() %>%
  set_engine("lm")

meta_workflow <- workflow() %>%
  add_recipe(meta_recipe) %>%
  add_model(lm_spec)

# Fit
meta_fit <- fit(meta_workflow, data = primary.data.str)

# Split conformal prediction intervals
set.seed(42)
split <- initial_split(primary.data.str, prop = 0.8)

conformal_split <- int_conformal_split(
  meta_fit,
  cal_data = testing(split)
)

# Predict with conformal intervals
new_predictions <- predict(conformal_split,
                           new_data = primary.data.str,
                           level = 0.95)
```

### 12.5 Small Sample Corrections

With few studies (< 20), conformal prediction can still have coverage variability. The **Small Sample Beta Correction (SSBC)** addresses this:

```r
#' Small Sample Beta Correction for Conformal Prediction
#'
#' Adjusts significance level to achieve probabilistic coverage guarantee
#' Reference: https://arxiv.org/abs/2509.15349

ssbc_correction <- function(n_cal, alpha, delta = 0.1) {
  #' @param n_cal Number of calibration samples

  #' @param alpha Nominal miscoverage rate

  #' @param delta Probability of coverage failure (default 0.1)
  #'
  #' @return Adjusted alpha for conformal prediction

  # The coverage of split conformal follows Beta(n_cal + 1 - k, k)

  # where k = ceiling((n_cal + 1) * alpha)
  # SSBC finds adjusted alpha such that P(coverage >= 1 - alpha) >= 1 - delta

  # Simple approximation for small samples
  # Use more conservative quantile
  adjusted_alpha <- alpha * (1 - 2 / sqrt(n_cal))
  adjusted_alpha <- max(adjusted_alpha, alpha / 2)  # Don't over-correct

  return(adjusted_alpha)
}

# Usage
n_calibration <- 20  # Small calibration set
nominal_alpha <- 0.05
adjusted_alpha <- ssbc_correction(n_calibration, nominal_alpha)

cat("Nominal alpha:", nominal_alpha, "\n")
cat("Adjusted alpha:", round(adjusted_alpha, 4), "\n")
# Use adjusted_alpha in conformal prediction for better coverage
```

### 12.6 Comparison: Parametric vs Conformal PIs

```r
#' Compare Parametric and Conformal Prediction Intervals
#'
#' Visualization and coverage comparison

compare_prediction_intervals <- function(model, cp_result, data) {

  # Parametric PI from metafor
  preds_parametric <- predict(model, addx = TRUE)
  tau2_total <- sum(model$sigma2)
  t_crit <- qt(0.975, df = model$k - model$p)
  pi_width_param <- 2 * t_crit * sqrt(preds_parametric$se^2 + tau2_total)

  # Conformal PI
  pi_width_conformal <- 2 * cp_result$q_hat

  # Comparison
  comparison <- data.frame(
    Method = c("Parametric", "Conformal"),
    Mean_PI_Width = c(mean(pi_width_param), pi_width_conformal),
    Coverage_Guarantee = c("Asymptotic (assumes normality)",
                           "Finite-sample (distribution-free)")
  )

  # Visualization
  plot_data <- data.frame(
    yi = data$yi,
    pred = fitted(model),
    pi_lb_param = fitted(model) - pi_width_param/2,
    pi_ub_param = fitted(model) + pi_width_param/2,
    pi_lb_conf = fitted(model) - cp_result$q_hat,
    pi_ub_conf = fitted(model) + cp_result$q_hat
  ) %>%
    arrange(pred)

  plot_data$index <- 1:nrow(plot_data)

  p <- ggplot(plot_data, aes(x = index)) +
    geom_ribbon(aes(ymin = pi_lb_param, ymax = pi_ub_param),
                fill = "blue", alpha = 0.2) +
    geom_ribbon(aes(ymin = pi_lb_conf, ymax = pi_ub_conf),
                fill = "red", alpha = 0.2) +
    geom_line(aes(y = pred), color = "black") +
    geom_point(aes(y = yi), size = 1, alpha = 0.5) +
    labs(title = "Prediction Intervals: Parametric vs Conformal",
         subtitle = "Blue = Parametric, Red = Conformal",
         x = "Observation (sorted by prediction)",
         y = "Effect Size") +
    theme_minimal()

  return(list(comparison = comparison, plot = p))
}
```

### 12.7 When to Use Conformal Methods

#### 12.7.1 Strong Candidates for Conformal Prediction

| Scenario | Why Conformal Helps |
|----------|---------------------|
| **Few studies (k < 20)** | Finite-sample guarantees vs. asymptotic |
| **Non-normal random effects** | Distribution-free |
| **Outliers in effect sizes** | Robust to distributional violations |
| **Dose-response prediction** | Honest uncertainty in new dose regions |
| **Clinical decision-making** | Rigorous coverage guarantees |

#### 12.7.2 When Parametric May Be Preferred

| Scenario | Why Parametric Might Suffice |
|----------|------------------------------|
| **Many studies (k > 50)** | Asymptotic guarantees adequate |
| **Well-behaved data** | Normality holds |
| **Interpretability priority** | Parametric PIs more familiar |
| **Computational constraints** | Split conformal still requires calibration set |

### 12.8 Reporting Conformal Results

```markdown
## Prediction Intervals

### Method Comparison

| Method | Formula | Width | Coverage Guarantee |
|--------|---------|-------|-------------------|
| Parametric (HTS) | est ± t√(SE² + τ²) | X.XX | Asymptotic (normality assumed) |
| Split Conformal | est ± q̂ | X.XX | Finite-sample (≥ 95%, distribution-free) |
| Conformalized QR | [q̂_lo - E, q̂_hi + E] | X.XX (adaptive) | Finite-sample (≥ 95%, distribution-free) |

### Conformal Prediction Details
- Calibration set size: n = X
- Nonconformity score: Absolute residual
- Empirical quantile (q̂): X.XX
- Small sample correction applied: [Yes/No]

### Interpretation
The conformal 95% prediction interval [X.XX, X.XX] indicates that with at least 95%
probability, a new study drawn from the same population would have an effect size
within this range. Unlike parametric intervals, this guarantee holds in finite
samples without distributional assumptions.
```

### 12.9 Limitations and Considerations

| Limitation | Description | Mitigation |
|------------|-------------|------------|
| **Exchangeability required** | Observations must be exchangeable | Use hierarchical CP methods |
| **Marginal coverage only** | No conditional coverage guarantee | CQR provides adaptive widths |
| **Calibration set size** | Need sufficient calibration data | SSBC for small samples |
| **Computational cost** | Full CP very expensive | Use split CP or Jackknife+ |
| **Novel in meta-analysis** | Less established than parametric | Report both for comparison |
| **Interpretation** | Less familiar to reviewers | Provide clear explanation |

### 12.10 R Packages for Conformal Prediction

| Package | Description | Use Case |
|---------|-------------|----------|
| `probably` | tidymodels integration | General regression |
| `conformal` | Full and split conformal | General regression |
| `cfcausal` | Conformal for causal inference | Treatment effects |
| `MAPIE` (Python) | Model Agnostic PI Estimation | If using Python |

```r
# Install packages
install.packages("probably")
# devtools::install_github("ryantibs/conformal")  # If available
```

---

## 13. Implementation Priority & Roadmap

### 13.1 Priority Matrix

| Priority | Enhancement | Effort | Impact | Timeline |
|----------|-------------|--------|--------|----------|
| **1** | Publication bias assessment | Low | High | Week 1 |
| **2** | Assumption testing | Low | Medium | Week 1 |
| **3** | Sensitivity analyses | Low | High | Week 2 |
| **4** | Prediction intervals (enhance) | Low | Medium | Week 2 |
| **5** | Cross-validation | Medium | Medium | Week 3 |
| **6** | DAG development | Low | Medium | Week 3 |
| **7** | Bayesian reanalysis | Medium | High | Week 4-5 |
| **8** | Fractional polynomials | Low | Medium | Week 5 |
| **9** | Robust SE enhancements | Low | Low | Week 5 |
| **10** | **Conformal prediction** | Medium | High | Week 5-6 |
| **11** | IPD meta-analysis | High | Very High | Long-term |

### 13.2 Quick Wins (Implement Immediately)

1. **Add funnel plot + Egger's test** (< 1 hour)
2. **Add Q-Q plot of residuals** (< 30 min)
3. **Add leave-one-out analysis** (< 1 hour)
4. **Create sensitivity analysis table** (2-3 hours)
5. **Draw DAG using DAGitty** (1-2 hours)

### 13.3 Medium-Term Enhancements

1. **Bayesian meta-regression with brms** (1-2 days)
2. **5-fold cross-validation** (half day)
3. **Selection model for publication bias** (half day)
4. **Fractional polynomial comparison** (half day)

### 13.4 Long-Term Projects

1. **IPD meta-analysis** - requires data sharing agreements
2. **Meta-analytic SEM** - requires additional correlation data
3. **Prospective meta-analysis** - future study design

---

## 14. Additional R Packages Required

```r
# Publication bias
install.packages("RoBMA")

# Bayesian meta-analysis
install.packages("brms")
install.packages("cmdstanr")  # For Stan backend

# Causal inference
install.packages("ggdag")
install.packages("dagitty")

# Robust estimation
install.packages("clubSandwich")
install.packages("robumeta")

# Fractional polynomials
install.packages("mfp")

# Meta-analytic SEM
install.packages("metaSEM")

# Additional diagnostics
install.packages("car")

# Conformal prediction
install.packages("probably")
install.packages("quantreg")
# devtools::install_github("ryantibs/conformal")  # If available

# Load all
pacman::p_load(
  RoBMA, brms, ggdag, dagitty, clubSandwich, robumeta,
  mfp, metaSEM, car, probably, quantreg
)
```

---

## 15. References

### Publication Bias

- Bartoš, F., et al. (2022). Adjusting for Publication Bias in JASP and R: Selection Models, PET-PEESE, and Robust Bayesian Meta-Analysis. *Advances in Methods and Practices in Psychological Science*, 5(3). https://doi.org/10.1177/25152459221109259

- Mathur, M. B., & VanderWeele, T. J. (2024). P-hacking in meta-analyses: A formalization and new meta-analytic methods. *Research Synthesis Methods*. https://doi.org/10.1002/jrsm.1701

### Assumption Testing

- Wang, C. C., & Lee, W. C. (2023). The normality assumption on between-study random effects was questionable in a considerable number of Cochrane meta-analyses. *BMC Medicine*, 21, 112. https://doi.org/10.1186/s12916-023-02823-9

- Meng, X., et al. (2024). Sensitivity analysis with iterative outlier detection for systematic reviews and meta-analyses. *Statistics in Medicine*. https://doi.org/10.1002/sim.10008

### Bayesian Methods

- Williams, D. R., Rast, P., & Bürkner, P. C. (2018). Bayesian meta-analysis with weakly informative prior distributions. *PsyArXiv*. https://doi.org/10.31234/osf.io/7tbrm

- Bürkner, P. C. (2017). brms: An R Package for Bayesian Multilevel Models Using Stan. *Journal of Statistical Software*, 80(1). https://doi.org/10.18637/jss.v080.i01

### Causal Inference

- Dekkers, O. M., et al. (2024). From complexity to clarity: how directed acyclic graphs enhance the study design of systematic reviews and meta-analyses. *Clinical Epidemiology*, 16, 33-40. https://doi.org/10.2147/CLEP.S433587

- Tennant, P. W., et al. (2021). Use of directed acyclic graphs (DAGs) to identify confounders in applied health research: review and recommendations. *International Journal of Epidemiology*, 50(2), 620-632. https://doi.org/10.1093/ije/dyaa213

### Cross-Validation

- Willis, B. H., & Riley, R. D. (2017). Measuring the statistical validity of summary meta-analysis and meta-regression results for use in clinical practice. *Statistics in Medicine*, 36(21), 3283-3301. https://doi.org/10.1002/sim.7372

### Dose-Response

- Crippa, A., & Orsini, N. (2016). Dose-response meta-analysis of differences in means. *BMC Medical Research Methodology*, 16, 91. https://doi.org/10.1186/s12874-016-0189-0

- Bagnardi, V., et al. (2004). Flexible Meta-Regression Functions for Modeling Aggregate Dose-Response Data. *American Journal of Epidemiology*, 159(11), 1077-1086. https://doi.org/10.1093/aje/kwh142

### IPD Meta-Analysis

- Riley, R. D., et al. (2021). Individual Participant Data Meta-Analysis: A Handbook for Healthcare Research. *Wiley*. ISBN: 978-1119333722

- Burke, D. L., et al. (2017). Meta-analysis using individual participant data: one-stage and two-stage approaches. *Research Synthesis Methods*, 8(2), 204-220. https://doi.org/10.1002/jrsm.1259

### Prediction Intervals

- IntHout, J., et al. (2016). Plea for routinely presenting prediction intervals in meta-analysis. *BMJ Open*, 6(7), e010247. https://doi.org/10.1136/bmjopen-2015-010247

### Heterogeneity

- Harrer, M., et al. (2021). Doing Meta-Analysis with R: A Hands-On Guide. *Chapman & Hall/CRC Press*. https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/

### Robust Estimation

- Pustejovsky, J. E., & Tipton, E. (2022). Meta-analysis with Robust Variance Estimation: Expanding the Range of Working Models. *Prevention Science*, 23, 425-438. https://doi.org/10.1007/s11121-021-01246-3

### Conformal Prediction

- Angelopoulos, A. N., & Bates, S. (2022). A Gentle Introduction to Conformal Prediction and Distribution-Free Uncertainty Quantification. *arXiv preprint*. https://arxiv.org/abs/2107.07511

- Lei, J., G'Sell, M., Rinaldo, A., Tibshirani, R. J., & Wasserman, L. (2018). Distribution-Free Predictive Inference for Regression. *Journal of the American Statistical Association*, 113(523), 1094-1111. https://doi.org/10.1080/01621459.2017.1307116

- Dunn, R., Wasserman, L., & Ramdas, A. (2022). Distribution-Free Prediction Sets for Two-Layer Hierarchical Models. *Journal of the American Statistical Association*, 118(544), 2651-2662. https://doi.org/10.1080/01621459.2022.2060112

- Kaiser, T., & Herzog, P. (2025). A Tutorial on Distribution-Free Uncertainty Quantification Using Conformal Prediction. *Advances in Methods and Practices in Psychological Science*. https://doi.org/10.1177/25152459251380452

- Gasparin, A., et al. (2025). Conformal Prediction for Hierarchical Data. *arXiv preprint*. https://arxiv.org/abs/2411.13479

- Small Sample Beta Correction (SSBC). (2024). Probabilistic Conformal Coverage Guarantees in Small-Data Settings. *arXiv preprint*. https://arxiv.org/abs/2509.15349

---

## Appendix: Complete Enhanced Analysis Script Template

```r
# =============================================================================
# ENHANCED PTF META-ANALYSIS SCRIPT
# =============================================================================

# Load additional packages
pacman::p_load(
  # Original packages
  readxl, tidyverse, metafor, emmeans, reshape2, janitor, gt, writexl,
  splines, performance, ggdist, ggridges, orchaRd, parameters,
  ggplot2, patchwork, bayestestR,


  # New packages for enhancements
  RoBMA,        # Robust Bayesian meta-analysis
  brms,         # Bayesian meta-regression
  ggdag,        # DAG visualization
  dagitty,      # DAG analysis
  clubSandwich, # Robust variance estimation
  car,          # Diagnostics
  mfp,          # Fractional polynomials
  probably,     # Conformal prediction (tidymodels)
  quantreg      # Quantile regression for CQR
)

# =============================================================================
# After fitting primary model, add these sections:
# =============================================================================

# --- PUBLICATION BIAS ASSESSMENT ---
source("scripts/01_publication_bias.R")

# --- ASSUMPTION TESTING ---
source("scripts/02_assumption_testing.R")

# --- SENSITIVITY ANALYSES ---
source("scripts/03_sensitivity_analyses.R")

# --- CROSS-VALIDATION ---
source("scripts/04_cross_validation.R")

# --- BAYESIAN REANALYSIS (optional) ---
source("scripts/05_bayesian_analysis.R")

# --- DAG ANALYSIS ---
source("scripts/06_dag_analysis.R")

# --- CONFORMAL PREDICTION ---
source("scripts/07_conformal_prediction.R")

# =============================================================================
```

This modular approach keeps the main script clean while adding comprehensive enhancements.
