# ============================================================
# SIMULATION STUDY (UK-PARAMETERIZED):
# Cup-to-Disc Ratio as an Effect Modifier of GLP-1-Associated
# NAION Risk — United Kingdom Parameters
#
# Author: Holli Groff, OD, MPH
# Date:   November 2026
#
# Research Question:
#   Does pre-existing optic disc anatomy (cup-to-disc ratio)
#   modify the risk of NAION in patients taking a GLP-1
#   receptor agonist, using UK-specific incidence and
#   prescribing parameters?
#
# Design:
#   Simulated retrospective cohort study.
#   N = 50,000 UK adults with type 2 diabetes or obesity.
#   Outcome: NAION diagnosis within 3-year follow-up.
#
# Biological rationale:
#   - NAION is associated with "disc at risk" anatomy
#     (crowded optic disc, small cup-to-disc ratio)
#   - GLP-1 agonists may reduce optic nerve perfusion
#   - Hypothesis: CDR x GLP-1 interaction — small CDR
#     amplifies GLP-1-associated NAION risk
#
# UK-specific parameters:
#   - Background NAION incidence ~2.3/100,000/year
#     (Hattenhauer et al. adjusted for UK case ascertainment;
#      UK estimate lower than US due to stricter clinical criteria
#      and under-reporting in primary care)
#   - GLP-1 main effect HR ~2.17
#     (PLOS Medicine meta-analysis, 2025 — same evidence base)
#   - Ozempic (semaglutide) approved UK: January 2019 (MHRA)
#   - Wegovy (high-dose semaglutide) approved UK: June 2023 (MHRA)
#   - GLP-1 exposure ~20% by 2026 in eligible UK T2DM population
#     (lower than US due to NHS prescribing restrictions and
#      NICE guidance thresholds)
#   - UK T2DM prevalence: ~8% adults (Diabetes UK, 2024)
#   - Hypertension prevalence: ~30% UK adults (NHS Digital, 2023)
#   - CDR distribution: mean 0.4, SD 0.12 (Jonas et al. —
#     consistent across European populations)
#   - Small CDR (<=0.3) prevalence ~15%
#   - Interaction term: hypothesized OR ~3.0 for small CDR
#     x GLP-1 (novel — this is what we are estimating)
#
# Proposed real-world data source for confirmatory study:
#   Clinical Practice Research Datalink (CPRD) — UK primary
#   care EHR database covering ~20% of the UK population,
#   linked to HES (Hospital Episode Statistics) for secondary
#   care diagnoses including ophthalmology.
#   Alternative: UK Biobank (N~500,000; optic disc imaging
#   available in ~67,000 participants).
# ============================================================

# ── 0. Packages ──────────────────────────────────────────────
packages <- c("dplyr", "ggplot2", "patchwork", "broom", "tidyr")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

set.seed(42)

# ── Colour palette ───────────────────────────────────────────
COL_RED   <- "#c0392b"
COL_BLUE  <- "#2c7bb6"
COL_GREEN <- "#27ae60"
COL_ORNG  <- "#e67e22"
COL_DARK  <- "#1a1a2e"
COL_MID   <- "#4a5568"
COL_GRAY  <- "#f4f4f4"

# ============================================================
# SECTION 1: SIMULATE PATIENT POPULATION (UK parameters)
# ============================================================
N <- 500000

# logistic helper
logistic <- function(x) 1 / (1 + exp(-x))

# --- Demographics (UK T2DM/obesity population) ---
# Mean age slightly older than US cohort (NHS T2DM registry data)
age  <- pmin(pmax(rnorm(N, 60, 12), 30), 85)
male <- rbinom(N, 1, 0.55)   # T2DM slightly more prevalent in men in UK

# --- GLP-1 exposure ---
# Lower than US: NICE restricts semaglutide to BMI>=35 (or >=30 with
# comorbidities) and HbA1c>=58mmol/mol. Wegovy approved June 2023.
# Estimated ~20% of eligible UK T2DM patients on GLP-1 by 2026.
glp1_prob <- logistic(-2.0 + 0.015 * age + 0.2 * male)
glp1      <- rbinom(N, 1, glp1_prob)

# --- Cup-to-disc ratio (truncated normal) ---
# Same distribution: consistent across European populations (Jonas et al.)
cdr <- pmin(pmax(rnorm(N, 0.40, 0.12), 0.05), 0.95)

# CDR categories (clinical thresholds — same as UK optometry practice)
cdr_cat <- case_when(
  cdr <= 0.30 ~ "Small (<=0.3)",
  cdr <= 0.60 ~ "Normal (0.31-0.60)",
  TRUE        ~ "Large (>0.60)"
)
cdr_cat <- factor(cdr_cat,
                  levels = c("Normal (0.31-0.60)",   # reference
                             "Small (<=0.3)",
                             "Large (>0.60)"))

cdr_small <- as.integer(cdr <= 0.30)
cdr_large <- as.integer(cdr >  0.60)

# --- Covariates (UK-specific prevalences) ---
htn      <- rbinom(N, 1, 0.30)   # hypertension: ~30% UK adults (NHS Digital 2023)
smoker   <- rbinom(N, 1, 0.13)   # smoking: ~13% UK adults (ONS 2023)
diabetes <- rbinom(N, 1, 0.75)   # T2D (rest: obesity indication for GLP-1)

# ============================================================
# SECTION 2: SIMULATE NAION OUTCOME (UK incidence)
#
# UK NAION incidence: ~2.3/100,000/year
# 3-year cumulative incidence: ~6.9/100,000 ~ 0.000069
# logit(0.000069) ~ -9.57
#
# Same biological parameters as US model:
# GLP-1 HR 2.17, interaction OR 3.0 (hypothesized)
# ============================================================

intercept        <- -9.57          # calibrated to UK incidence 2.3/100k/yr
beta_glp1        <-  log(2.17)    # main effect of GLP-1 (same evidence base)
beta_cdr_small   <-  log(2.50)    # small CDR baseline risk (anatomical)
beta_cdr_large   <-  log(0.70)    # large CDR (lower baseline)
beta_interaction <-  log(3.00)    # KEY: interaction (novel hypothesis)
beta_age         <-  log(1.35)    # per 10 years older
beta_htn         <-  log(1.60)    # hypertension
beta_smoke       <-  log(1.25)    # smoking

log_odds <- intercept +
  beta_glp1        * glp1 +
  beta_cdr_small   * cdr_small +
  beta_cdr_large   * cdr_large +
  beta_interaction * glp1 * cdr_small +
  beta_age         * (age - 60) / 10 +
  beta_htn         * htn +
  beta_smoke       * smoker

prob_naion <- logistic(log_odds)
naion      <- rbinom(N, 1, prob_naion)

# ── Build analysis data frame ────────────────────────────────
df <- data.frame(
  naion, glp1, cdr, cdr_cat, cdr_small, cdr_large,
  age, male, htn, smoker, diabetes,
  interaction = glp1 * cdr_small,
  age_c       = (age - 60) / 10
)

# ============================================================
# SECTION 3: DESCRIPTIVE STATISTICS
# ============================================================
cat("\n", strrep("=", 60), "\n")
cat("SIMULATION STUDY (UK): CDR x GLP-1 Effect Modification\n")
cat("on NAION Risk — UK-Parameterized\n")
cat(strrep("=", 60), "\n\n")
cat(sprintf("Total N:       %d\n", N))
cat(sprintf("GLP-1 users:   %d (%.1f%%)\n",
            sum(glp1), mean(glp1) * 100))
cat(sprintf("NAION events:  %d (%.1f per 100,000)\n\n",
            sum(naion), mean(naion) * 100000))
cat("CDR distribution:\n")
cdr_table <- table(df$cdr_cat)
for (nm in names(cdr_table)) {
  cat(sprintf("  %s: %d (%.1f%%)\n",
              nm, cdr_table[nm], cdr_table[nm] / N * 100))
}
cat("\nNAION by CDR category:\n")
naion_by_cdr <- df %>%
  group_by(cdr_cat) %>%
  summarise(events = sum(naion),
            rate   = mean(naion) * 100000, .groups = "drop")
for (i in seq_len(nrow(naion_by_cdr))) {
  cat(sprintf("  %s: %d events (%.1f per 100,000)\n",
              naion_by_cdr$cdr_cat[i],
              naion_by_cdr$events[i],
              naion_by_cdr$rate[i]))
}

# ============================================================
# SECTION 4: LOGISTIC REGRESSION WITH INTERACTION TERM
# ============================================================
model <- glm(
  naion ~ glp1 + cdr_small + cdr_large +
    I(glp1 * cdr_small) +
    age_c + htn + smoker,
  data   = df,
  family = binomial(link = "logit")
)

cat("\n", strrep("=", 60), "\n")
cat("LOGISTIC REGRESSION RESULTS\n")
cat("Outcome: NAION (3-year cumulative incidence)\n")
cat(strrep("=", 60), "\n\n")

tidy_model <- tidy(model, conf.int = TRUE, exponentiate = TRUE)
tidy_model <- tidy_model %>%
  mutate(
    term = recode(term,
                  "(Intercept)"          = "Intercept",
                  "glp1"                 = "GLP-1",
                  "cdr_small"            = "CDR small",
                  "cdr_large"            = "CDR large",
                  "I(glp1 * cdr_small)"  = "GLP-1 x CDR small",
                  "age_c"                = "Age (+10yr)",
                  "htn"                  = "Hypertension",
                  "smoker"               = "Smoking"
    ),
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )

cat(sprintf("%-22s %6s  %18s  %9s\n", "Variable", "OR", "95% CI", "p-value"))
cat(strrep("-", 60), "\n")
for (i in seq_len(nrow(tidy_model))) {
  r <- tidy_model[i, ]
  if (r$term == "Intercept") {
    cat(sprintf("%-22s %6s  %18s  %9.4f\n", r$term, "--", "--", r$p.value))
  } else {
    ci_str <- sprintf("(%.2f-%.2f)", r$conf.low, r$conf.high)
    cat(sprintf("%-22s %6.2f  %18s  %9.4f %s\n",
                r$term, r$estimate, ci_str, r$p.value, r$sig))
  }
}

# ============================================================
# SECTION 5: STRATIFIED ANALYSIS
# ============================================================
cat("\n", strrep("=", 60), "\n")
cat("STRATIFIED ANALYSIS: GLP-1 NAION Risk by CDR Category\n")
cat(strrep("=", 60), "\n")

strata_list <- list(
  list(label = "Small CDR (<=0.3)",       mask = df$cdr_small == 1),
  list(label = "Normal CDR (0.31-0.60)",  mask = (df$cdr_small == 0) & (df$cdr_large == 0)),
  list(label = "Large CDR (>0.60)",       mask = df$cdr_large == 1)
)

strata_results <- list()
for (s in strata_list) {
  sub <- df[s$mask, ]
  a  <- sum(sub$glp1 == 1 & sub$naion == 1)
  b  <- sum(sub$glp1 == 1 & sub$naion == 0)
  cc <- sum(sub$glp1 == 0 & sub$naion == 1)
  d  <- sum(sub$glp1 == 0 & sub$naion == 0)
  if (a > 0 && b > 0 && cc > 0 && d > 0) {
    or_strat <- (a * d) / (b * cc)
    se_log   <- sqrt(1/a + 1/b + 1/cc + 1/d)
    lo       <- exp(log(or_strat) - 1.96 * se_log)
    hi       <- exp(log(or_strat) + 1.96 * se_log)
    rate_g   <- mean(sub$naion[sub$glp1 == 1]) * 100000
    rate_c   <- mean(sub$naion[sub$glp1 == 0]) * 100000
    cat(sprintf("\n  %s (n=%d, events=%d)\n", s$label, nrow(sub), sum(sub$naion)))
    cat(sprintf("    GLP-1 users:  %.1f per 100,000\n", rate_g))
    cat(sprintf("    Non-GLP-1:    %.1f per 100,000\n", rate_c))
    cat(sprintf("    Crude OR: %.2f (95%% CI %.2f-%.2f)\n", or_strat, lo, hi))
    strata_results[[length(strata_results) + 1]] <- list(
      label     = s$label,
      or        = or_strat,
      lo        = lo,
      hi        = hi,
      n         = nrow(sub),
      events    = sum(sub$naion),
      rate_glp1 = rate_g,
      rate_ctrl = rate_c
    )
  }
}

strata_df <- bind_rows(lapply(strata_results, as.data.frame))

# ============================================================
# SECTION 6: SAMPLE SIZE CALCULATION (UK context)
# ============================================================
cat("\n", strrep("=", 60), "\n")
cat("SAMPLE SIZE CALCULATION\n")
cat("Target: 80% power to detect interaction OR >= 3.0\n")
cat("alpha = 0.05, two-sided | P(small CDR) ~ 15%\n")
cat("UK NAION incidence: 2.3/100,000/year\n")
cat(strrep("=", 60), "\n\n")

z_alpha <- qnorm(0.975)
z_beta  <- qnorm(0.80)
p_glp1  <- 0.20            # lower GLP-1 penetration in UK vs US
p_small <- 0.15
p_naion <- 0.0000230 * 3   # UK incidence x 3 years
or_int  <- 3.0

numerator   <- (z_alpha + z_beta)^2
denominator <- p_glp1 * (1 - p_glp1) *
  p_small * (1 - p_small) *
  log(or_int)^2 *
  p_naion * (1 - p_naion) * 4

n_needed <- ceiling(numerator / denominator)
cat(sprintf("  Estimated N needed: %s\n", formatC(n_needed, format="d", big.mark=",")))
cat(sprintf("  (Interaction OR = %.1f, alpha = 0.05, power = 80%%)\n", or_int))
cat(sprintf("  Expected NAION events: ~%d\n", ceiling(n_needed * p_naion * 3)))
cat("\n  -> Feasible with CPRD (N~20M) linked to HES ophthalmology data\n")
cat("  -> UK Biobank retinal imaging subcohort (~67,000) for CDR phenotyping\n")

# ============================================================
# SECTION 7: VISUALIZATIONS
# ============================================================

# ── Plot 1: Stratified forest plot ───────────────────────────
all_strata_cols <- c(
  "Small CDR (<=0.3)"      = COL_RED,
  "Normal CDR (0.31-0.60)" = COL_BLUE,
  "Large CDR (>0.60)"      = COL_GREEN
)

strata_df$label <- factor(strata_df$label, levels = rev(strata_df$label))
strata_df$color <- all_strata_cols[as.character(strata_df$label)]

p1 <- ggplot(strata_df, aes(x = or, y = label, colour = label)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             colour = "gray50", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = lo, xmax = hi),
                 height = 0.2, linewidth = 0.9) +
  geom_point(size = 4, shape = 18) +
  scale_colour_manual(values = all_strata_cols, guide = "none") +
  geom_text(aes(x = hi, y = label,
                label = sprintf("OR %.2f (%.2f-%.2f)\nN=%s | Events=%d",
                                or, lo, hi,
                                formatC(n, format = "d", big.mark = ","),
                                events)),
            hjust = -0.05, size = 3.2, colour = COL_DARK,
            inherit.aes = FALSE, data = strata_df) +
  scale_x_log10(limits = c(0.3, 60),
                breaks  = c(0.5, 1, 2, 4, 8, 16, 32),
                labels  = c("0.5","1","2","4","8","16","32")) +
  labs(
    title    = "GLP-1 Associated NAION Risk Stratified by Cup-to-Disc Ratio",
    subtitle = "UK-Parameterized Simulation | N = 50,000 | Incidence: 2.3/100,000/yr (MHRA context)",
    x        = "Odds Ratio for GLP-1 on NAION (log scale)",
    y        = NULL,
    caption  = paste0(
      "UK parameters: NAION incidence 2.3/100,000/yr | GLP-1 HR 2.17 (PLOS Med meta-analysis 2025)\n",
      "Ozempic MHRA approval: Jan 2019 | Wegovy MHRA approval: Jun 2023 | ",
      "Proposed confirmatory data: CPRD + HES"
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13, colour = COL_DARK),
    plot.subtitle    = element_text(size = 10, colour = COL_MID),
    plot.caption     = element_text(size = 7.5, colour = "gray50"),
    panel.grid.minor = element_blank()
  )

# ── Plot 2: Interaction plot ──────────────────────────────────
interaction_long <- strata_df %>%
  select(label, rate_glp1, rate_ctrl) %>%
  tidyr::pivot_longer(
    cols      = c(rate_glp1, rate_ctrl),
    names_to  = "group",
    values_to = "rate"
  ) %>%
  mutate(group = recode(group,
                        "rate_glp1" = "GLP-1 users",
                        "rate_ctrl" = "Non-GLP-1 (diabetes meds)"
  ))

p2 <- ggplot(interaction_long,
             aes(x = label, y = rate, colour = group, group = group)) +
  geom_line(linewidth = 2) +
  geom_point(aes(shape = group), size = 4) +
  geom_text(aes(label = sprintf("%.0f", rate)),
            nudge_y = 2, size = 3.5, fontface = "bold") +
  scale_colour_manual(values = c("GLP-1 users" = COL_RED,
                                 "Non-GLP-1 (diabetes meds)" = COL_BLUE)) +
  scale_shape_manual(values  = c("GLP-1 users" = 16,
                                 "Non-GLP-1 (diabetes meds)" = 15)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.25))) +
  labs(
    title    = "Interaction Plot: NAION Incidence by GLP-1 Use and Cup-to-Disc Ratio",
    subtitle = "Diverging lines indicate effect modification (CDR as moderator) — UK parameters",
    x        = "Cup-to-Disc Ratio Category",
    y        = "NAION Incidence (per 100,000 person-years)",
    colour   = NULL,
    shape    = NULL,
    caption  = "UK-parameterized simulation | N = 50,000 | 3-year follow-up | Background incidence: 2.3/100,000/yr"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13, colour = COL_DARK),
    plot.subtitle    = element_text(size = 10, colour = COL_MID),
    plot.caption     = element_text(size = 8, colour = "gray50"),
    legend.position  = "top",
    panel.grid.minor = element_blank()
  )

# ── Plot 3: Regression coefficient OR plot ───────────────────
reg_df <- tidy_model %>%
  filter(term != "Intercept") %>%
  mutate(
    term = factor(term, levels = rev(c(
      "GLP-1", "CDR small", "CDR large",
      "GLP-1 x CDR small", "Age (+10yr)",
      "Hypertension", "Smoking"
    ))),
    var_type = case_when(
      grepl("x", term)   ~ "Interaction",
      grepl("CDR", term) ~ "CDR covariate",
      TRUE               ~ "Other covariate"
    ),
    label = sprintf("OR %.2f (%.2f-%.2f) %s",
                    estimate, conf.low, conf.high, sig)
  )

p3 <- ggplot(reg_df, aes(x = estimate, y = term, colour = var_type)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             colour = "gray50", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.25, linewidth = 0.9) +
  geom_point(size = 4, shape = 16) +
  geom_text(aes(x = conf.high, label = label),
            hjust = -0.05, size = 3.0, colour = COL_DARK) +
  scale_colour_manual(values = c(
    "Interaction"     = COL_RED,
    "CDR covariate"   = COL_ORNG,
    "Other covariate" = COL_BLUE
  )) +
  scale_x_log10(limits = c(0.3, 40),
                breaks  = c(0.5, 1, 2, 4, 8, 16),
                labels  = c("0.5","1","2","4","8","16")) +
  labs(
    title    = "Multivariable Logistic Regression: Predictors of NAION",
    subtitle = "Including GLP-1 x CDR Interaction Term — UK-Parameterized Simulation",
    x        = "Odds Ratio (95% CI, log scale)",
    y        = NULL,
    colour   = NULL,
    caption  = paste0(
      "N = 500,000 simulated UK patients | *** p<0.001, ** p<0.01, * p<0.05, ns = not significant\n",
      "UK parameters: NAION 2.3/100k/yr, GLP-1 exposure ~20%, HTN ~30% (NHS Digital 2023)"
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13, colour = COL_DARK),
    plot.subtitle    = element_text(size = 10, colour = COL_MID),
    plot.caption     = element_text(size = 7.5, colour = "gray50"),
    legend.position  = "top",
    panel.grid.minor = element_blank()
  )

# ── Save plots ────────────────────────────────────────────────
ggsave("naion_cdr_stratified_forest_UK.png", plot = p1,
       width = 11, height = 5, dpi = 150)
ggsave("naion_cdr_interaction_plot_UK.png",  plot = p2,
       width = 9,  height = 6, dpi = 150)
ggsave("naion_cdr_regression_UK.png",        plot = p3,
       width = 10, height = 6, dpi = 150)

cat("\n===== Plots saved =====\n")
cat("1. naion_cdr_stratified_forest_UK.png\n")
cat("2. naion_cdr_interaction_plot_UK.png\n")
cat("3. naion_cdr_regression_UK.png\n")
cat("\nDone!\n")
