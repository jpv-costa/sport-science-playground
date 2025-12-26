# Meta-Analysis Experiment Plan: Proximity to Failure and Resistance Training Outcomes

> **Document Version:** 1.0
> **Based on:** `V2.Analysis.Script.R` (Pelland et al.)
> **Last Updated:** December 2025

---

## 1. Overview

### 1.1 Research Question

**Does training closer to muscular failure (lower Repetitions in Reserve - RIR) produce greater strength and hypertrophy adaptations compared to stopping further from failure?**

### 1.2 Study Type

Systematic Review and Multi-Level Meta-Regression Analysis

### 1.3 Outcomes of Interest

| Outcome | Description |
|---------|-------------|
| **Maximal Strength** | Changes in 1RM or maximal voluntary contraction |
| **Muscle Hypertrophy** | Changes in muscle size (CSA, thickness, lean mass) |

---

## 2. Rationale (Why)

### 2.1 Scientific Gap

Training to muscular failure is a common practice in resistance training, but its necessity for optimal adaptations remains debated. Understanding the dose-response relationship between proximity to failure (measured as RIR) and training outcomes can:

1. **Optimize training prescription** for athletes and general populations
2. **Reduce unnecessary fatigue accumulation** that may impair recovery
3. **Inform evidence-based programming decisions** for coaches and practitioners
4. **Identify potential threshold effects** (e.g., is there a point of diminishing returns?)

### 2.2 Hypotheses

1. Training closer to failure (lower RIR values) will be associated with greater effect sizes for both strength and hypertrophy
2. The relationship may follow a non-linear pattern (logarithmic or spline), suggesting diminishing returns as RIR approaches zero
3. The effect may be moderated by training status, load, volume equating method, and other study characteristics

---

## 3. Methods (How)

### 3.1 Data Collection

#### Inclusion Criteria (Implied from Data Structure)

- Randomized controlled trials or within-subject designs
- Resistance training interventions with defined or estimable proximity to failure
- Pre-post measurements of strength and/or hypertrophy outcomes
- Sufficient statistical information to calculate effect sizes (means, SDs, SEs, CIs, or p-values)

#### Data Source

Primary data extracted into: `V2.PTF.Data.xlsx`

Supplementary file for RIR estimation: `Supplementary File 1.RIR Estimation Sheet.xlsx`

### 3.2 Effect Size Calculation

Two effect size metrics are calculated for robustness:

| Metric | Code | Interpretation |
|--------|------|----------------|
| **SMCR** (Standardized Mean Change with Raw score standardization) | `escalc(measure = "SMCR")` | Hedge's g - standardized units |
| **ROMC** (Ratio of Means with Correlation) | `escalc(measure = "ROMC")` | Log response ratio - converts to % change via `(exp(x)-1)*100` |

#### 3.2.1 Handling Missing Statistics

The script implements a cascade of imputations:

```
1. SE → SD:          SD = SE × sqrt(n)
2. p-value → t:      t = qt(p/2, df=n-1, lower.tail=FALSE)
3. t → SE:           SE = delta.mean / t
4. CI → SE:          SE = (CI.upper - CI.lower) / 3.92
5. SE → SD:          SD = SE × sqrt(n)
```

#### 3.2.2 Pre-Post Correlation Imputation

When pre-post correlations (`ri`) are not reported:

1. **Calculate from available SDs:**
   ```
   ri = (pre.sd² + post.sd² - delta.sd²) / (2 × pre.sd × post.sd)
   ```

2. **Validate range:** Remove values outside [-1, +1]

3. **Meta-analytic imputation:**
   - Convert r to Fisher's z
   - Fit random-effects model to estimate pooled z
   - Convert back to r using `psych::fisherz2r()`
   - Impute separately for Strength and Hypertrophy outcomes

4. **Estimate change score SD:**
   ```
   delta.sd = sqrt(pre.sd² + post.sd² - 2×ri×pre.sd×post.sd)
   ```

### 3.3 Statistical Models

#### 3.3.1 Multi-Level Meta-Regression

All models use `rma.mv()` from the `metafor` package with:

**Random Effects Structure:**
```r
random = list(~1|study, ~1|group, ~1|obs)
```

| Level | Description |
|-------|-------------|
| **Study** | Accounts for clustering of effects within studies |
| **Group** | Accounts for multiple groups within studies |
| **Observation** | Accounts for multiple outcomes per group (residual) |

**Estimation Method:** REML (Restricted Maximum Likelihood)

**Inference:** t-distribution with containment degrees of freedom (`dfs = "contain"`)

#### 3.3.2 Functional Forms Tested

Six candidate models for the RIR-outcome relationship:

| Model | Specification | Rationale |
|-------|---------------|-----------|
| **Linear** | `avg.rir` | Simple linear dose-response |
| **Log-linear** | `log1p(avg.rir)` | Diminishing returns at higher RIR |
| **Quadratic** | `poly(avg.rir, 2, raw=TRUE)` | U-shaped or inverted-U relationship |
| **Cubic** | `poly(avg.rir, 3, raw=TRUE)` | Complex non-linear patterns |
| **Linear Spline** | `bs(spline.rir, degree=1, knots=0)` | Threshold effect at failure (knot at 0) |
| **Restricted Cubic Spline** | `rcs(avg.rir, 4)` | Flexible non-linear fit with 4 knots |

**Spline Variable:** `spline.rir` codes failure as -1, otherwise equals `avg.rir`

#### 3.3.3 Random Effects Variants

Each functional form is tested with three random structures:

| Variant | Specification | Description |
|---------|---------------|-------------|
| `ri` | `~1\|study, ~1\|group, ~1\|obs` | Random intercepts only |
| `rs` | `~avg.rir\|study, ~1\|group, ~1\|obs` | Random slopes at study level |
| `rs2` | `~avg.rir\|study, ~avg.rir\|group, ~1\|obs` | Random slopes at study and group levels |

**Total models per outcome:** 6 forms × 3 structures = **18 models**

#### 3.3.4 Model Selection

**Bayes Factor Comparison:**
- Use `bf_models()` function (from `bayestestR` package)
- Compare all 18 models pairwise
- Visualize as heatmap matrix

**Kass & Raftery (1995) Interpretation Scale:**

| 2×log(BF) | Evidence Strength |
|-----------|-------------------|
| < 0 | Favors denominator |
| 0 - 2 | Weak |
| 2 - 6 | Positive |
| 6 - 10 | Strong |
| > 10 | Very Strong |

#### 3.3.5 Final Model Estimation

1. **Refit best model** with robust variance estimation:
   ```r
   robust.rma.mv(model, cluster = study)
   ```

2. **Calculate heterogeneity statistics:**
   - I² via `i2_ml()`
   - R² via `r2_ml()`

3. **Compute prediction intervals:**
   ```r
   t_crit × sqrt(SE² + sum(tau²))
   ```

#### 3.3.6 Covariates in All Models

| Covariate | Variable | Description |
|-----------|----------|-------------|
| Training load | `load.set` | Average load per set (% of 1RM) |
| Volume equating | `set.rep.equated` | How volume was matched (set/rep/both) |
| Duration | `weeks` | Intervention length |
| Training status | `train.status` | Trained vs. untrained |

### 3.4 Moderator Analyses

Interaction terms are added to the best-fit model:
```r
update(best_model, ~. + log1p(avg.rir):moderator)
```

#### 3.4.1 Training Variables

| Moderator | Levels/Values | Analysis Type |
|-----------|---------------|---------------|
| Load (% 1RM) | 30%, 60%, 90% | Continuous interaction |
| Volume equating method | Set, Repetition, Both | Categorical interaction |
| Adjusted sets per week | Continuous | Continuous interaction |
| Training frequency | Sessions/week/muscle | Continuous interaction |
| Intervention duration | Weeks | Continuous interaction |

#### 3.4.2 Participant Variables

| Moderator | Levels/Values | Analysis Type |
|-----------|---------------|---------------|
| Training status | Trained, Untrained | Categorical interaction |
| Age | Continuous (years) | Continuous interaction |
| Sex (% male) | 0-100% | Continuous interaction |

#### 3.4.3 Design Variables

| Moderator | Levels | Analysis Type |
|-----------|--------|---------------|
| Study design | Within, Between | Categorical interaction |
| Body region | Upper, Lower | Categorical interaction |
| Exercise type | Multi-joint, Single-joint | Categorical interaction |

#### 3.4.4 Intervention Characteristics

| Moderator | Levels | Analysis Type |
|-----------|--------|---------------|
| Concurrent training | Yes, No | Categorical interaction |
| Progressive overload | Yes, No | Categorical interaction |
| Failure definition provided | Yes, No | Categorical interaction |
| Alternative set structures | Yes, No | Categorical interaction |
| Maximal intended velocity | Yes, No | Categorical interaction |

#### 3.4.5 Strength-Specific Moderators

| Moderator | Levels | Analysis Type |
|-----------|--------|---------------|
| Strength test type | Isometric, Isokinetic, Isotonic | Categorical interaction |
| Test type | 1RM, Non-1RM | Categorical interaction |

### 3.5 Marginal Effects

For each model, compute:

1. **Marginal Means:** Predicted effect sizes at RIR = 0, 1, 2, ..., 23
   ```r
   emmprep() %>% emmeans(~avg.rir, weights = "prop")
   ```

2. **Marginal Slopes:** Effect of 1-unit change in RIR at mean values
   ```r
   emmeans() %>% contrast("consec")
   ```

3. **Prediction Intervals:** Added via custom `pred_interval_esmeans()` function

### 3.6 Sensitivity Analysis

**Secondary Analysis:** Restricted to studies with **RIR ≤ 10**

- Excludes extreme values that may represent different training paradigms
- Tests robustness of primary findings
- Focuses on practically relevant training ranges

---

## 4. Expected Data Schema

### 4.1 Required Input Variables

#### Identification Variables

| Variable | Type | Description | Example |
|----------|------|-------------|---------|
| `author.year` | Character | Study citation | "Smith 2020" |
| `study` | Factor | Unique study ID | 1, 2, 3... |
| `group` | Factor | Group within study | 1, 2, 3... |
| `obs` | Factor | Observation ID | 1, 2, 3... |
| `outcome` | Factor | Outcome type | "Strength", "Hypertrophy" |

#### Descriptive Statistics (Pre-Post)

| Variable | Type | Required | Description |
|----------|------|----------|-------------|
| `n` | Integer | Yes | Sample size |
| `pre.mean` | Numeric | Yes | Pre-intervention mean |
| `post.mean` | Numeric | Yes | Post-intervention mean |
| `pre.sd` | Numeric | Yes* | Pre-intervention SD |
| `post.sd` | Numeric | Yes* | Post-intervention SD |
| `pre.se` | Numeric | No | Pre-intervention SE |
| `post.se` | Numeric | No | Post-intervention SE |

*Can be computed from SE if not available

#### Change Score Statistics

| Variable | Type | Required | Description |
|----------|------|----------|-------------|
| `delta.mean` | Numeric | No | Mean change score |
| `delta.sd` | Numeric | No | SD of change |
| `delta.se` | Numeric | No | SE of change |
| `delta.CI.lower` | Numeric | No | Lower 95% CI |
| `delta.CI.upper` | Numeric | No | Upper 95% CI |
| `t.value` | Numeric | No | t-statistic |
| `p.value` | Numeric | No | p-value |

#### Primary Predictor

| Variable | Type | Required | Description |
|----------|------|----------|-------------|
| `avg.rir` | Numeric | Yes | Average estimated RIR across sets |
| `rir.bucket` | Factor | Yes | Categorical RIR (includes "Failure") |

#### Training Characteristics

| Variable | Type | Required | Description |
|----------|------|----------|-------------|
| `load.set` | Numeric | Yes | Load per set (% 1RM) |
| `adj.sets.week` | Numeric | Yes | Adjusted sets/week/muscle |
| `frequency.per.muscle` | Numeric | Yes | Weekly frequency |
| `weeks` | Numeric | Yes | Intervention duration |
| `train.status` | Factor | Yes | "Trained" or "Untrained" |

#### Study Design Variables

| Variable | Type | Required | Description |
|----------|------|----------|-------------|
| `set.rep.equated` | Factor | Yes | "set", "rep", or "both" |
| `within.between.design` | Factor | Yes | Study design type |
| `upper.lower.other` | Factor | Yes | Body region |
| `train.exercise` | Factor | Yes | Exercise type |
| `test.exercise` | Factor | No | Testing exercise |

#### Participant Characteristics

| Variable | Type | Required | Description |
|----------|------|----------|-------------|
| `age` | Numeric | Yes | Mean age (years) |
| `sex.percent.male` | Numeric | Yes | Proportion male (0-100) |

#### Intervention Details

| Variable | Type | Required | Description |
|----------|------|----------|-------------|
| `nutrition.controlled` | Factor | Yes | Yes/No |
| `formal.cardio.intervention` | Factor | Yes | Yes/No |
| `progression` | Factor | Yes | Progressive overload |
| `failure.definition` | Factor | Yes | Failure defined |
| `alternative.set.structure` | Factor | Yes | Non-traditional sets |
| `train.intent` | Factor | Yes | Maximal velocity intent |
| `cluster.sets` | Factor | No | Cluster set use |
| `rest.redistribution` | Factor | No | Rest redistribution |
| `rest.pause` | Factor | No | Rest-pause sets |

#### Strength-Specific Variables

| Variable | Type | Required | Description |
|----------|------|----------|-------------|
| `isometric.isokinetic.isotonic` | Factor | For strength | Test modality |
| `RM.max.submax` | Factor | For strength | 1RM or submaximal |
| `str.avg.rir` | Numeric | For strength | Strength-specific RIR |
| `str.specific.load` | Numeric | For strength | Strength-specific load |
| `adj.str.specific.sets.week` | Numeric | For strength | Strength sets/week |
| `frequency.per.strength.me` | Numeric | For strength | Strength frequency |

### 4.2 Computed Variables

| Variable | Formula/Function | Description |
|----------|------------------|-------------|
| `ri` | Calculated or meta-analytically imputed | Pre-post correlation |
| `yi` | `escalc(measure="SMCR")` | Hedge's g effect size |
| `vi` | `escalc(measure="SMCR")` | Variance of Hedge's g |
| `yi.rom` | `escalc(measure="ROMC")` | Log response ratio |
| `vi.rom` | `escalc(measure="ROMC")` | Variance of log RR |
| `weights` | `1/sqrt(vi)` | Precision weights (g) |
| `weights.rom` | `1/sqrt(vi.rom)` | Precision weights (ROM) |
| `spline.rir` | `ifelse(rir.bucket=="Failure", -1, avg.rir)` | For spline models |
| `yi.rom.exp` | `(exp(yi.rom)-1)*100` | % change |

---

## 5. Expected Outputs

### 5.1 Model Comparison Outputs

| Output | Format | Description |
|--------|--------|-------------|
| `primary.model.compare.str.g` | PDF | Bayes Factor heatmap (18×18 models) for strength (g) |
| `primary.model.compare.str.rom` | PDF | Bayes Factor heatmap for strength (ROM) |
| `primary.model.compare.hyp.g` | PDF | Bayes Factor heatmap for hypertrophy (g) |
| `primary.model.compare.hyp.rom` | PDF | Bayes Factor heatmap for hypertrophy (ROM) |

### 5.2 Model Coefficient Tables

| Output | Format | Contents |
|--------|--------|----------|
| `primary.str.g.best.output.table` | PNG | Coefficients, 95% CI, 95% PI |
| `primary.str.rom.best.output.table` | PNG | Coefficients, 95% CI, 95% PI |
| `primary.hyp.g.best.output.table` | PNG | Coefficients, 95% CI, 95% PI |
| `primary.hyp.rom.best.output.table` | PNG | Coefficients, 95% CI, 95% PI |

### 5.3 Heterogeneity Statistics

| Statistic | Variable | Description |
|-----------|----------|-------------|
| I² | `*.best.model.i2` | Proportion of variance from heterogeneity |
| R² | `*.best.model.r2` | Variance explained by moderators |
| τ² | From model output | Between-study variance components |

### 5.4 Marginal Effects Plots

For each outcome × effect size combination:

| Plot Type | Description | File Pattern |
|-----------|-------------|--------------|
| Main effect | RIR vs. effect size with CI and PI bands | `*.g.marginal.plot.pdf` |
| By moderator | Stratified dose-response curves | `*.g.moderator.*.plot.pdf` |

### 5.5 Moderator Analysis Outputs

For each moderator:
1. **Interaction contrast** with 95% CI (displayed in plot caption)
2. **Stratified marginal means** at key RIR values
3. **Faceted visualization** showing effect modification

### 5.6 Summary Tables

| Table | Contents | Format |
|-------|----------|--------|
| `primary.char.table.combined.final` | Distribution of categorical study characteristics | PNG |
| `primary.density.plot.final` | Density distributions of continuous variables | PDF |

### 5.7 Data Files

| File | Contents |
|------|----------|
| `primary.data.RData` | Cleaned primary analysis dataset |
| `secondary.data.RData` | Cleaned secondary analysis dataset (RIR ≤ 10) |
| `*.best.model.i2.RData` | I² statistics |
| `*.best.model.r2.RData` | R² statistics |

---

## 6. Analysis Workflow

```
┌─────────────────────────────────────────────────────────────────┐
│                     1. DATA PREPARATION                         │
├─────────────────────────────────────────────────────────────────┤
│  ├── Load raw data (V2.PTF.Data.xlsx)                          │
│  ├── Calculate missing SDs from SEs                            │
│  ├── Convert p-values → t-values → SEs                         │
│  ├── Calculate/impute pre-post correlations (ri)               │
│  │   ├── Direct calculation where possible                     │
│  │   ├── Fisher's z meta-analysis for imputation               │
│  │   └── Separate estimates for Strength/Hypertrophy           │
│  ├── Calculate effect sizes (SMCR, ROMC)                       │
│  ├── Create spline variable (failure = -1)                     │
│  ├── Subset by outcome (Strength, Hypertrophy)                 │
│  └── Set factor levels and reference categories                │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│                  2. PRIMARY ANALYSIS (Per Outcome)              │
├─────────────────────────────────────────────────────────────────┤
│  ├── Fit 18 candidate models                                   │
│  │   ├── 6 functional forms × 3 random structures              │
│  │   └── For both g and ROM effect sizes                       │
│  ├── Compare models using Bayes Factors                        │
│  ├── Select best-fitting model                                 │
│  ├── Refit with robust variance estimation                     │
│  ├── Calculate I² and R²                                       │
│  ├── Extract coefficients with CI and PI                       │
│  └── Generate model output tables                              │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│                    3. MARGINAL EFFECTS                          │
├─────────────────────────────────────────────────────────────────┤
│  ├── Calculate marginal means at RIR = 0:23                    │
│  ├── Calculate marginal slopes at mean RIR                     │
│  ├── Add prediction intervals                                  │
│  └── Generate dose-response plots                              │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│                   4. MODERATOR ANALYSES                         │
├─────────────────────────────────────────────────────────────────┤
│  ├── Add interaction terms to best model                       │
│  ├── Calculate interaction contrasts                           │
│  ├── Generate stratified marginal means                        │
│  └── Create moderator plots                                    │
│                                                                 │
│  Moderators analyzed:                                           │
│  ├── Load, Volume method, Training status, Weeks               │
│  ├── Adjusted sets, Frequency, Age, Sex                        │
│  ├── Study design, Body region, Exercise type                  │
│  ├── Concurrent training, Progressive overload                 │
│  ├── Failure definition, Set structure, Velocity intent        │
│  └── (Strength only) Test type, 1RM vs non-1RM                 │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│              5. SECONDARY ANALYSIS (RIR ≤ 10)                   │
├─────────────────────────────────────────────────────────────────┤
│  └── Repeat steps 2-4 with restricted dataset                  │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│                   6. OUTPUT GENERATION                          │
├─────────────────────────────────────────────────────────────────┤
│  ├── Save summary tables (PNG via gt)                          │
│  ├── Save plots (PDF via ggplot2/patchwork)                    │
│  └── Save model objects (RData)                                │
└─────────────────────────────────────────────────────────────────┘
```

---

## 7. Software Requirements

### 7.1 R Packages

| Package | Version | Purpose |
|---------|---------|---------|
| `metafor` | ≥4.0 | Meta-analysis models (`rma.mv`, `escalc`, `robust`) |
| `emmeans` | ≥1.8 | Marginal means and contrasts |
| `tidyverse` | ≥2.0 | Data manipulation and piping |
| `ggplot2` | ≥3.4 | Visualization |
| `patchwork` | ≥1.1 | Plot composition |
| `gt` | ≥0.9 | Publication-quality tables |
| `splines` | base | B-spline basis functions |
| `rms` | ≥6.0 | Restricted cubic splines |
| `orchaRd` | ≥2.0 | Meta-analysis visualization |
| `bayestestR` | ≥0.13 | Bayes Factor calculation |
| `ggridges` | ≥0.5 | Density ridge plots |
| `ggdist` | ≥3.0 | Distribution visualization |
| `readxl` | ≥1.4 | Excel file reading |
| `writexl` | ≥1.4 | Excel file writing |
| `janitor` | ≥2.2 | Data cleaning |
| `reshape2` | ≥1.4 | Data reshaping |
| `parameters` | ≥0.21 | Model parameters extraction |
| `performance` | ≥0.10 | Model diagnostics |
| `psych` | ≥2.3 | Fisher's z transformations |

### 7.2 Custom Functions

| Function | Purpose | Location |
|----------|---------|----------|
| `expo(x)` | Convert log RR to % change: `(exp(x)-1)*100` | In script |
| `pred_interval_esmeans()` | Add prediction intervals to emmeans | External/custom |
| `bf_models()` | Compute Bayes Factor model comparisons | From `bayestestR` |

---

## 8. Interpretation Guidelines

### 8.1 Effect Size Benchmarks (Hedge's g)

| Magnitude | Range | Interpretation |
|-----------|-------|----------------|
| Trivial | < 0.20 | Negligible practical significance |
| Small | 0.20 - 0.49 | Small but potentially meaningful |
| Medium | 0.50 - 0.79 | Moderate practical significance |
| Large | ≥ 0.80 | Large practical significance |

### 8.2 Response Ratio Interpretation

| % Change | Interpretation |
|----------|----------------|
| < 5% | Minimal change |
| 5-10% | Small change |
| 10-20% | Moderate change |
| > 20% | Large change |

### 8.3 Heterogeneity Interpretation (I²)

| I² Value | Interpretation |
|----------|----------------|
| 0-25% | Low heterogeneity |
| 25-50% | Moderate heterogeneity |
| 50-75% | Substantial heterogeneity |
| 75-100% | Considerable heterogeneity |

### 8.4 Key Comparisons

| Comparison | Practical Relevance |
|------------|---------------------|
| RIR = 0 vs RIR = 3 | Failure vs. practical "near failure" threshold |
| Log-linear slope | Marginal effect of increasing RIR by 1 unit |
| Spline knot at 0 | Tests for distinct effect below/above failure |
| Trained vs. Untrained | Population-specific recommendations |

---

## 9. Limitations of Current Approach

### 9.1 Methodological Limitations

1. **RIR Estimation Heterogeneity:** Studies use different methods to estimate/report RIR (self-report, velocity-based, repetitions to failure test)

2. **Pre-Post Correlation Imputation:** Meta-analytic imputation assumes correlations are similar across studies, which may not hold

3. **Ecological Fallacy:** Study-level moderators may not reflect individual-level effects (e.g., mean age doesn't capture age-related heterogeneity within studies)

4. **Publication Bias:** Not explicitly assessed in current protocol

5. **Covariate Selection:** Limited to available study-level variables; unmeasured confounders may exist

### 9.2 Statistical Limitations

1. **Model Selection Uncertainty:** Bayes Factor comparison doesn't account for model uncertainty in final inference

2. **Random Effects Normality:** Assumes normally distributed random effects, which may not hold

3. **Robust SE Limitations:** Cluster-robust SEs require sufficient number of clusters (studies) for valid inference

4. **Multiple Comparisons:** Many moderator analyses increase false positive risk

### 9.3 Practical Limitations

1. **Generalizability:** Findings limited to populations and training protocols represented in included studies

2. **Effect Size Precision:** Wide prediction intervals may limit practical recommendations

3. **Temporal Dynamics:** Cannot assess how RIR-outcome relationship changes over training duration

---

## 10. File Structure

```
deadlift-study/
├── pelland/
│   ├── V2.Analysis.Script.R          # Main analysis script
│   ├── V2.PTF.Data.xlsx              # Primary data file
│   └── Supplementary File 1.RIR Estimation Sheet.xlsx
├── docs/
│   ├── 01_CURRENT_EXPERIMENT_PLAN.md  # This document
│   └── 02_METHODOLOGICAL_IMPROVEMENTS.md
├── Plots/                             # Output directory for figures
├── Tables/                            # Output directory for tables
└── Models/                            # Output directory for model outputs
```

---

## 11. References

### Key Methodological References

1. Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. *Journal of Statistical Software*, 36(3), 1-48.

2. Hedges, L. V., Tipton, E., & Johnson, M. C. (2010). Robust variance estimation in meta-regression with dependent effect size estimates. *Research Synthesis Methods*, 1(1), 39-65.

3. Kass, R. E., & Raftery, A. E. (1995). Bayes factors. *Journal of the American Statistical Association*, 90(430), 773-795.

4. IntHout, J., Ioannidis, J. P., Rovers, M. M., & Goeman, J. J. (2016). Plea for routinely presenting prediction intervals in meta-analysis. *BMJ Open*, 6(7), e010247.

### Domain-Specific References

*(To be added based on included studies)*
