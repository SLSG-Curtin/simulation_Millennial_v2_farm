#!/usr/bin/env Rscript
# ==============================================================================
# OPTION 1: Per-site z-score ALL SEM variables (replace originals in-place)
# ==============================================================================
# This version expands zscore_cols to include all SEM predictors, so the
# group_by(site) %>% mutate(across(..., safe_z)) applies to all of them.
# The later global scale() in sem_df is REMOVED (no double-scaling).
# ==============================================================================

library(tidyverse)
library(lavaan)

set.seed(20241127)

# Helper functions
safe_z <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  sigma <- sd(x, na.rm = TRUE)
  if (is.na(sigma) || sigma == 0) return(rep(0, length(x)))
  (x - mu) / sigma
}

assign_npp_regime <- function(npp_ratio) {
  targets <- c(0.5, 1.0, 1.5)
  labels <- c("npp_0.5", "npp_1.0", "npp_1.5")
  diffs <- abs(outer(npp_ratio, targets, "-"))
  idx <- apply(diffs, 1, function(row) {
    if (all(is.na(row))) return(NA_integer_)
    row[is.na(row)] <- Inf
    which.min(row)
  })
  factor(labels[idx], levels = labels)
}

# Data
data_path <- "d04_output/All_interactions_1118.csv"
raw_df <- read_csv(data_path, show_col_types = FALSE)

maom_base <- raw_df %>%
  mutate(site = as.factor(site), MAOM_Qmax = MAOM - Qmax)

# ---------------------------------------------------------------------------- #
# OPTION 1 KEY CHANGE: expand zscore_cols to include all SEM variables
# ---------------------------------------------------------------------------- #
#zscore_cols <- c("MAOM", "MAOM_Qmax", "pH", "ClaySilt", "CUE", "ratio_constraints")
zscore_cols <- c("MAOM", "MAOM_Qmax")

normalized_df <- maom_base %>%
  group_by(site) %>%
  mutate(across(all_of(zscore_cols), safe_z)) %>%
  ungroup()

maom_df <- normalized_df %>%
  #transmute(MAOM, pH, ClaySilt, Npp_ratio, CUE, MAOM_Qmax, ratio_constraints, site) %>%
  #filter(!is.na(MAOM)) %>%
  drop_na(pH, ClaySilt, Npp_ratio, CUE, ratio_constraints) %>%
  mutate(Npp_regime = assign_npp_regime(Npp_ratio)) %>%
  filter(!is.na(Npp_regime))

# Outputs
dir.create("d05_fig/maom_influence", recursive = TRUE, showWarnings = FALSE)
dir.create("d04_output/maom_influence", recursive = TRUE, showWarnings = FALSE)
model_artifacts_dir <- "d04_output/maom_influence"
output_dir <- "d05_fig/maom_influence"

# ---------------------------------------------------------------------------- #
# Latent factor (SEM) analysis ----------------------------------------------- #
# ---------------------------------------------------------------------------- #
# NO additional global scaling here — variables are already per-site z-scored
sem_df <- normalized_df %>%
  select(MAOM, CUE, MAOM_Qmax, pH, ClaySilt, ratio_constraints, Npp_ratio) %>%
  #mutate(Npp_regime = forcats::fct_drop(Npp_regime))
  mutate(Npp_regime = assign_npp_regime(Npp_ratio)) %>%
  filter(!is.na(Npp_regime))

sem_df <- sem_df[,c(1:6,8)]

# Candidate models based on your specification:
model_specs <- list(
  
  # Option A1
  optionA1 = '
    Efficiency =~ CUE + ratio_constraints
    Sorption =~ MAOM_Qmax + ClaySilt + pH
    MAOM ~ Efficiency + Sorption
  ',
  
  # Option A2
  optionA2 = '
    Efficiency =~ CUE + ratio_constraints
    Sorption =~ MAOM_Qmax + ClaySilt + pH
    MAOM ~ Sorption + Efficiency
    Sorption ~~ Efficiency
  ',
  
  # Option A3
  optionA3 = '
    Sorption =~ MAOM_Qmax + ClaySilt + pH
    CUE ~ ratio_constraints
    MAOM ~ Sorption + CUE
    Sorption ~~ CUE
  ',
  
  # Option A4
  optionA4 = '
    Efficiency =~ CUE + ratio_constraints
    Sorption =~ MAOM_Qmax + ClaySilt + pH
    MAOM ~ Sorption + Efficiency
    CUE ~~ ratio_constraints
    Sorption ~~ Efficiency
  ',
  
  # Option A5
  optionA5 = '
    Sorption_1 =~ MAOM_Qmax + ClaySilt
    MAOM ~ Sorption_1 + CUE + pH
    CUE ~~ ratio_constraints
  ',
  
  # Option A6
  optionA6 = '
    Efficiency =~ CUE + ratio_constraints
    Sorption_1 =~ MAOM_Qmax + ClaySilt
    Desorption =~ 1*pH
    MAOM ~ Sorption_1 + Efficiency + Desorption
    CUE ~~ ratio_constraints
  ',
  
  # Option A7
  optionA7 = '
    Efficiency =~ CUE + ratio_constraints
    Sorption_1 =~ MAOM_Qmax + ClaySilt
    Desorption =~ 1*pH
    MAOM ~ Sorption_1 + Efficiency + Desorption
    CUE ~~ ratio_constraints
  ',
  
  # Option A8
  optionA8 = '
    Efficiency =~ CUE + ratio_constraints
    Sorption_1 =~ MAOM_Qmax + ClaySilt
    Desorption =~ 1*pH
    MAOM ~ Sorption_1 + Efficiency + Desorption
    CUE ~~ ratio_constraints
    MAOM_Qmax ~~ ClaySilt
  ',
  
  # Option B: Path model (no latents) — direct effects only
  optionB = '
    MAOM ~ pH + ClaySilt + CUE + ratio_constraints + MAOM_Qmax
  '
)

# Fit all candidate models and collect fit statistics and parameters
model_results <- list()
for (nm in names(model_specs)) {
  message('Fitting model: ', nm)
  m <- model_specs[[nm]]
  fit <- tryCatch(
    lavaan::sem(m, data = sem_df, group = 'Npp_regime', estimator = 'MLR', missing = 'fiml'),
    error = function(e) {
      message('  Failed: ', e$message)
      NULL
    }
  )
  if (!is.null(fit)) {
    measures <- tryCatch(lavaan::fitMeasures(fit, c('cfi','tli','rmsea','srmr','aic','bic')), error = function(e) NULL)
    params <- tryCatch(lavaan::parameterEstimates(fit, standardized = TRUE), error = function(e) tibble::tibble())
    readr::write_csv(params, file.path(model_artifacts_dir, paste0('sem_parameters_option1_', nm, '.csv')))
    model_results[[nm]] <- list(fit = fit, measures = measures)
  }
}

# Create comparison table
comparison <- purrr::imap_dfr(model_results, function(x, nm) {
  ms <- x$measures
  tibble::tibble(
    model = nm,
    cfi = as.numeric(ms['cfi']),
    tli = as.numeric(ms['tli']),
    rmsea = as.numeric(ms['rmsea']),
    srmr = as.numeric(ms['srmr']),
    aic = as.numeric(ms['aic']),
    bic = as.numeric(ms['bic'])
  )
})
readr::write_csv(comparison, file.path(model_artifacts_dir, 'sem_model_comparison_option1.csv'))
print(comparison)

# Select a model to produce the plot (prefer 'optionA5' or first available)
chosen <- if ('optionA3' %in% names(model_results)) 'optionA3' else names(model_results)[1]
if (!is.null(chosen) && chosen %in% names(model_results)) {
  sem_params <- lavaan::parameterEstimates(model_results[[chosen]]$fit, standardized = TRUE)
  readr::write_csv(sem_params, file.path(model_artifacts_dir, paste0('sem_parameters_option1_', chosen, '.csv')))
  
  sem_reg_plot <- sem_params %>%
    filter(op == "~", lhs == "MAOM") %>%
    mutate(
      regime = levels(sem_df$Npp_regime)[group],
      regime = forcats::fct_recode(
        factor(regime, levels = levels(sem_df$Npp_regime)),
        "NPP 0.5" = "npp_0.5",
        "NPP 1.0" = "npp_1.0",
        "NPP 1.5" = "npp_1.5"
      )
    ) %>%
    ggplot(aes(x = regime, y = std.all, fill = rhs)) +
    geom_col(position = position_dodge(width = 0.6)) +
    geom_errorbar(aes(ymin = std.all - se, ymax = std.all + se), position = position_dodge(width = 0.6), width = 0.2) +
    labs(
      x = "NPP regime",
      y = "Std. effect on MAOM",
      fill = "Predictor",
      title = paste0("SEM effects on MAOM (Option 1 - ", chosen, ")")
    )
  
  ggsave(file.path(output_dir, "maom_sem_effects_option1.png"), sem_reg_plot, width = 6, height = 4.5, dpi = 320)
  print(sem_reg_plot)
} else {
  message('No SEM models converged successfully; check model specs and data.')
}

message("Option 1 script complete.")
message("Outputs saved to: ", model_artifacts_dir)
