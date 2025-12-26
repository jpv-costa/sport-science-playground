# scripts/replicate_pelland.R
# Replication of Pelland et al. (2024) Meta-Regression Analysis
# "Exploring the Dose-Response Relationship Between Estimated Resistance
#  Training Proximity to Failure, Strength Gain, and Muscle Hypertrophy"
#
# Original OSF: https://osf.io/7knsj/

# ==============================================================================
# Setup
# ==============================================================================

# Load packages
library(readxl)
library(tidyverse)
library(metafor)
library(data.table)

cat("=== Pelland et al. Meta-Regression Replication ===\n\n")

# ==============================================================================
# Data Loading and Preprocessing
# ==============================================================================

# Read data from OSF download
data_path <- "data/external/pelland_meta_regression/Data__Code__and_Estimation_Materials/V2.PTF.Data.xlsx"
df <- read_xlsx(data_path)
cat("Loaded", nrow(df), "rows from", data_path, "\n")

# Calculate pre-post SDs from SEs where needed
df$pre.sd <- ifelse(is.na(df$pre.se), df$pre.sd, df$pre.se * sqrt(df$n))
df$post.sd <- ifelse(is.na(df$post.se), df$post.sd, df$post.se * sqrt(df$n))

# Convert p to t (change scores)
df$t.value <- metafor::replmiss(df$t.value, with(df, qt(p.value/2, df=n-1, lower.tail=FALSE)))

# Convert t to se (change scores)
df$delta.se <- metafor::replmiss(df$delta.se, with(df, ifelse(is.na(delta.mean),
                                                    (post.mean - pre.mean)/t.value,
                                                    delta.mean/t.value)))
# Make positive
df$delta.se <- ifelse(df$delta.se < 0, df$delta.se * -1, df$delta.se)

# Convert CI to SE (change scores)
df$delta.se <- metafor::replmiss(df$delta.se, with(df, (delta.CI.upper - delta.CI.lower)/3.92))

# Convert SE to SD (change scores)
df$delta.sd <- metafor::replmiss(df$delta.sd, with(df, delta.se * sqrt(n)))

# Calculate pre-post correlation coefficient
df$ri <- (df$pre.sd^2 + df$post.sd^2 - df$delta.sd^2) / (2 * df$pre.sd * df$post.sd)

# Remove values outside the range of -1 to +1
df$ri <- ifelse(dplyr::between(df$ri, -1, 1) == FALSE, NA, df$ri)

cat("Pre-post correlations calculated\n")
cat("  Valid correlations:", sum(!is.na(df$ri)), "\n")
cat("  Missing correlations:", sum(is.na(df$ri)), "\n\n")

# ==============================================================================
# Correlation Imputation via Meta-Analysis
# ==============================================================================

cat("=== Meta-Analysis of Pre-Post Correlations ===\n")

# Convert to Fisher's z for meta-analysis
df.for.ri <- escalc(measure = "ZCOR", ri = ri, ni = n, data = df)

# Separate meta-analyses for Strength and Hypertrophy
Meta_ri_str <- rma.mv(yi, vi,
                      data = df.for.ri[df.for.ri$outcome == "Strength", ],
                      slab = paste(author.year),
                      random = list(~ 1 | study/group/obs),
                      method = "REML", test = "t")

Meta_ri_hyp <- rma.mv(yi, vi,
                      data = df.for.ri[df.for.ri$outcome == "Hypertrophy", ],
                      slab = paste(author.year),
                      random = list(~ 1 | study/group/obs),
                      method = "REML", test = "t")

# Robust variance estimates
RobuEstMeta_ri_str <- robust(Meta_ri_str, df.for.ri[df.for.ri$outcome == "Strength", ]$study)
RobuEstMeta_ri_hyp <- robust(Meta_ri_hyp, df.for.ri[df.for.ri$outcome == "Hypertrophy", ]$study)

# Convert back to r
z2r_str <- psych::fisherz2r(RobuEstMeta_ri_str$b[1])
z2r_hyp <- psych::fisherz2r(RobuEstMeta_ri_hyp$b[1])

cat("Imputed correlation for Strength:", round(z2r_str, 4), "\n")
cat("Imputed correlation for Hypertrophy:", round(z2r_hyp, 4), "\n\n")

# Impute missing correlations
df[df$outcome == "Strength", ]$ri <- ifelse(is.na(df[df$outcome == "Strength", ]$ri),
                                             z2r_str, df[df$outcome == "Strength", ]$ri)
df[df$outcome == "Hypertrophy", ]$ri <- ifelse(is.na(df[df$outcome == "Hypertrophy", ]$ri),
                                                z2r_hyp, df[df$outcome == "Hypertrophy", ]$ri)

# ==============================================================================
# Effect Size Calculation
# ==============================================================================

cat("=== Effect Size Calculation ===\n")

# Estimate change score SD where only pre-post data available
df$delta.sd <- metafor::replmiss(df$delta.sd,
                                  with(df, sqrt(pre.sd^2 + post.sd^2 - (2*ri*pre.sd*post.sd))))

# SMCR (Standardized Mean Change using Raw score standardization)
data <- escalc(measure = "SMCR",
               m1i = post.mean, m2i = pre.mean,
               sd1i = pre.sd, ni = n, ri = ri,
               data = df)

# ROMC (Response Ratio of Means, log scale)
data2 <- escalc(measure = "ROMC",
                m1i = post.mean, m2i = pre.mean,
                sd1i = post.sd, sd2i = pre.sd,
                ni = n, ri = ri,
                data = df)

# Rename ROM columns
data2 <- data2 %>%
  rename(yi.rom = yi, vi.rom = vi, ri.rom = ri)

# Combine
data <- cbind(data, subset(data2, select = c(ri.rom, yi.rom, vi.rom)))
data$weights <- (1 / sqrt(data$vi))
data$weights.rom <- (1 / sqrt(data$vi.rom))

# RIR for spline models
data$spline.rir <- ifelse(data$rir.bucket == "Failure", -1, data$avg.rir)

# Exponentiate ROM for % change
data$yi.rom.exp <- (exp(data$yi.rom) - 1) * 100

cat("SMCR effect sizes calculated:", sum(!is.na(data$yi)), "\n")
cat("ROMC effect sizes calculated:", sum(!is.na(data$yi.rom)), "\n\n")

# ==============================================================================
# Data Subsetting
# ==============================================================================

# Primary data (all RIR values)
primary.data.str <- data %>%
  filter(outcome == "Strength") %>%
  mutate(across(c(author.year, study, group, obs, outcome, train.status,
                  set.rep.equated), factor))

primary.data.hyp <- data %>%
  filter(outcome == "Hypertrophy") %>%
  mutate(across(c(author.year, study, group, obs, outcome, train.status,
                  set.rep.equated), factor))

cat("=== Data Summary ===\n")
cat("Strength outcomes:", nrow(primary.data.str), "effect sizes\n")
cat("Hypertrophy outcomes:", nrow(primary.data.hyp), "effect sizes\n\n")

# ==============================================================================
# Primary Meta-Regression Models (Strength)
# ==============================================================================

cat("=== Primary Meta-Regression: Strength (Hedge's g) ===\n\n")

# Linear model with random intercepts
primary.str.g.linear.model.ri <- rma.mv(
  yi, vi,
  data = primary.data.str,
  mods = ~ avg.rir + load.set + set.rep.equated + weeks + train.status,
  random = list(~1|study, ~1|group, ~1|obs),
  method = "REML", test = "t", dfs = "contain"
)

cat("Linear Model Results:\n")
print(summary(primary.str.g.linear.model.ri))

cat("\n\nRIR coefficient interpretation:\n")
cat("  For each 1 RIR increase (further from failure):\n")
cat("  Change in Hedge's g:", round(coef(primary.str.g.linear.model.ri)["avg.rir"], 4), "\n")
str_pval <- primary.str.g.linear.model.ri$pval[which(names(coef(primary.str.g.linear.model.ri)) == "avg.rir")]
cat("  p-value:", format.pval(str_pval, digits = 4), "\n")

# ==============================================================================
# Primary Meta-Regression Models (Hypertrophy)
# ==============================================================================

cat("\n\n=== Primary Meta-Regression: Hypertrophy (Hedge's g) ===\n\n")

primary.hyp.g.linear.model.ri <- rma.mv(
  yi, vi,
  data = primary.data.hyp,
  mods = ~ avg.rir + load.set + set.rep.equated + weeks + train.status,
  random = list(~1|study, ~1|group, ~1|obs),
  method = "REML", test = "t", dfs = "contain"
)

cat("Linear Model Results:\n")
print(summary(primary.hyp.g.linear.model.ri))

cat("\n\nRIR coefficient interpretation:\n")
cat("  For each 1 RIR increase (further from failure):\n")
cat("  Change in Hedge's g:", round(coef(primary.hyp.g.linear.model.ri)["avg.rir"], 4), "\n")
hyp_pval <- primary.hyp.g.linear.model.ri$pval[which(names(coef(primary.hyp.g.linear.model.ri)) == "avg.rir")]
cat("  p-value:", format.pval(hyp_pval, digits = 4), "\n")

# ==============================================================================
# Summary
# ==============================================================================

cat("\n\n=== REPLICATION SUMMARY ===\n")
cat("Successfully replicated Pelland et al. meta-regression analysis\n\n")

cat("Key findings:\n")
cat("1. Strength outcomes (n =", nrow(primary.data.str), "):\n")
cat("   RIR effect on Hedge's g:", round(coef(primary.str.g.linear.model.ri)["avg.rir"], 4), "\n")
cat("   p-value:", format.pval(primary.str.g.linear.model.ri$pval[2], digits = 4), "\n")

cat("\n2. Hypertrophy outcomes (n =", nrow(primary.data.hyp), "):\n")
cat("   RIR effect on Hedge's g:", round(coef(primary.hyp.g.linear.model.ri)["avg.rir"], 4), "\n")
cat("   p-value:", format.pval(primary.hyp.g.linear.model.ri$pval[2], digits = 4), "\n")

# Save results
results <- list(
  strength_model = primary.str.g.linear.model.ri,
  hypertrophy_model = primary.hyp.g.linear.model.ri,
  data_strength = primary.data.str,
  data_hypertrophy = primary.data.hyp,
  correlation_strength = z2r_str,
  correlation_hypertrophy = z2r_hyp
)

saveRDS(results, "data/processed/pelland_replication_results.rds")
cat("\nResults saved to data/processed/pelland_replication_results.rds\n")
