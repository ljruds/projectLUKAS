#' Fit GAMM models for cloud fraction analysis
#'
#' @param data Data frame containing all variables
#' @param covariates Character vector of covariate names
#' @param drought_start Start date (YYYY-MM-DD)
#' @param drought_end End date (YYYY-MM-DD)
#' @param months Months to include (default April–October)
#' @param blc_values Boundary layer cloud flags to retain
#' @param k_smooth Default basis dimension for smooths
#' @param make_plots Logical; return gratia plots
#'
#' @return List containing models, summaries, diagnostics, and plots
#' @export

fit_cloud_gamm <- function(
    data,
    covariates,
    drought_start = "2025-08-01",
    drought_end   = "2025-09-25",
    months = 4:10,
    blc_values = c(1, 2),
    k_smooth = 8,
    make_plots = TRUE
) {

  require(mgcv)
  require(dplyr)
  require(lubridate)
  require(tidyr)

  if (make_plots) require(gratia)

  #---------------------------
  # 1. Input checks
  #---------------------------
  required_cols <- c("date", "mean_cf_2_6", "blc_flag", covariates)
  missing_cols <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }

  #---------------------------
  # 2. Preprocessing
  #---------------------------
  eps <- 0.001

  df <- data %>%
    mutate(
      date = as.Date(date),
      month = lubridate::month(date),
      cloud_beta = pmin(pmax(mean_cf_2_6, eps), 1 - eps),
      drought = ifelse(
        date >= as.Date(drought_start) &
          date <= as.Date(drought_end),
        "Drought", "Normal"
      ),
      drought = as.factor(drought)
    ) %>%
    filter(month %in% months) %>%
    filter(blc_flag %in% blc_values) %>%
    drop_na(all_of(c("cloud_beta", covariates, "month", "drought")))

  #---------------------------
  # 3. Formulas
  #---------------------------
  formula_base <- cloud_beta ~
    s(median_PBLH, k = 12) +
    s(LE_f_daytime_sum, k = k_smooth) +
    s(SW_IN_daytime_sum, k = k_smooth) +
    s(evening_cbh, k = k_smooth) +
    s(LapseRate_850_500_CperKm, k = k_smooth) +
    PW_mm +
    s(RH_1_7_1_daytime_mean, k = k_smooth) +
    s(G_2_1_1_daytime_sum, k = k_smooth) +
    s(TA_1_2_1_daytime_mean, k = k_smooth) +
    PA_1_2_1_daytime_mean +
    H_f_daytime_sum +
    mean_LCL_4_7 +
    LTS_K +
    SWC_3_1_1_daytime_mean +
    s(month, bs = "re")

  formula_interaction <- cloud_beta ~
    drought +
    s(LE_f_daytime_sum, by = drought, k = k_smooth) +
    s(G_2_1_1_daytime_sum, by = drought, k = k_smooth) +
    s(H_f_daytime_sum, by = drought, k = 12) +
    SWC_3_1_1_daytime_mean * drought +
    s(median_PBLH, by = drought, k = 12) +
    s(SW_IN_daytime_sum, k = k_smooth) +
    s(evening_cbh, k = k_smooth) +
    s(LapseRate_850_500_CperKm, k = k_smooth) +
    PW_mm +
    s(RH_1_7_1_daytime_mean, by = drought, k = k_smooth) +
    s(TA_1_2_1_daytime_mean, by = drought, k = k_smooth) +
    PA_1_2_1_daytime_mean +
    mean_LCL_4_7 +
    LTS_K +
    s(month, bs = "re")

  #---------------------------
  # 4. Fit models
  #---------------------------
  mod_base <- gam(formula_base, data = df,
                  family = betar(link = "logit"), method = "REML")

  mod_interaction <- gam(formula_interaction, data = df,
                         family = betar(link = "logit"), method = "REML")

  #---------------------------
  # 5. Summaries
  #---------------------------
  capture_summary <- function(model) {
    paste(capture.output(summary(model)), collapse = "\n")
  }

  summaries <- list(
    base = capture_summary(mod_base),
    interaction = capture_summary(mod_interaction)
  )

  #---------------------------
  # 6. Diagnostics
  #---------------------------
  diagnostics <- list(
    base = list(
      gam_check = capture.output(gam.check(mod_base)),
      concurvity = mgcv::concurvity(mod_base, full = TRUE)
    ),
    interaction = list(
      gam_check = capture.output(gam.check(mod_interaction)),
      concurvity = mgcv::concurvity(mod_interaction, full = TRUE)
    )
  )

  #---------------------------
  # 7. Model comparison (AIC)
  #---------------------------
  aic_vals <- AIC(mod_base, mod_interaction)

  delta_aic <- aic_vals$AIC - min(aic_vals$AIC)
  weights <- exp(-0.5 * delta_aic) / sum(exp(-0.5 * delta_aic))

  model_comparison <- data.frame(
    Model = c("Base", "Interaction"),
    AIC = aic_vals$AIC,
    Delta_AIC = delta_aic,
    Weight = weights
  )

  #---------------------------
  # 8. Results table (AUTO)
  #---------------------------
  extract_results <- function(model, model_name) {

    s <- summary(model)

    # Parametric
    param <- as.data.frame(s$p.table)
    param$term <- rownames(param)
    param$type <- "parametric"

    # Smooth
    smooth <- as.data.frame(s$s.table)
    smooth$term <- rownames(smooth)
    smooth$type <- "smooth"

    # Harmonize column names
    names(param) <- c("estimate", "std_error", "stat", "p_value", "term", "type")
    names(smooth) <- c("edf", "ref_df", "stat", "p_value", "term", "type")

    param$model <- model_name
    smooth$model <- model_name

    bind_rows(param, smooth)
  }

  results_table <- bind_rows(
    extract_results(mod_base, "Base"),
    extract_results(mod_interaction, "Interaction")
  ) %>%
    mutate(
      significance = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE ~ ""
      )
    )

  #---------------------------
  # 9. Plots
  #---------------------------
  plots <- NULL

  if (make_plots) {
    plots <- list(
      base = gratia::draw(mod_base),
      interaction = gratia::draw(mod_interaction),
      interaction_focus = gratia::draw(
        mod_interaction,
        select = "s(H_f_daytime_sum)",
        partial_match = TRUE
      )
    )
  }

  #---------------------------
  # 10. Return (S3 object)
  #---------------------------
  out <- list(
    data = df,
    models = list(base = mod_base, interaction = mod_interaction),
    summaries = summaries,
    diagnostics = diagnostics,
    model_comparison = model_comparison,
    results_table = results_table,
    plots = plots
  )

  class(out) <- "cloud_gamm"
  return(out)
}
