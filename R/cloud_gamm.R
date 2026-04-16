#' Fit GAMM for cloud fraction with structured outputs
#'
#' @param data Data frame
#' @param formula_mixed Optional full GAM formula
#' @param smooth_vars Character vector of smooth terms
#' @param linear_vars Character vector of linear terms
#' @param response Response variable (default: "mean_cf_2_6")
#' @param date_var Date column name
#' @param blc_filter Logical; filter BLC days
#' @param blc_values Values for BLC filter
#' @param months Months to include
#' @param k_list Named list of k values for smooths
#' @param output_plots Logical
#' @param plot_file File path for PDF output
#'
#' @return List with model, summary, diagnostics, plots
#' @export
fit_cloud_gamm <- function(
    data=arf_daily,
    formula_mixed = NULL,
    smooth_vars = c("median_PBLH", "LE_f_daytime_sum", "SW_IN_daytime_sum",
                    "evening_cbh", "LapseRate_850_500_CperKm",
                    "RH_1_7_1_daytime_mean", "G_2_1_1_daytime_sum",
                    "TA_1_2_1_daytime_mean"),
    linear_vars = c("PA_1_2_1_daytime_mean", "H_f_daytime_sum",
                    "mean_LCL_4_7", "LTS_K", "PW_mm",
                    "SWC_3_1_1_daytime_mean"),
    response = "mean_cf_2_6",
    date_var = "date",
    blc_filter = TRUE,
    blc_values = c(1, 2),
    months = 4:10,
    k_list = list(default = 8, median_PBLH = 12),
    output_plots = TRUE,
    plot_file = "gamm_output.pdf"
) {

  require(mgcv)
  require(dplyr)

  df <- data

  #---------------------------
  # 1. Preprocessing
  #---------------------------
  eps <- 0.001

  df[[date_var]] <- as.Date(df[[date_var]])

  df <- df %>%
    mutate(
      month = lubridate::month(.data[[date_var]]),
      cloud_beta = pmin(pmax(.data[[response]], eps), 1 - eps)
    ) %>%
    filter(month %in% months)

  if (blc_filter) {
    df <- df %>% filter(blc_flag %in% blc_values)
  }

  # Drop NA
  vars_needed <- unique(c("cloud_beta", smooth_vars, linear_vars, "month"))
  df <- df %>% drop_na(all_of(vars_needed))

  #---------------------------
  # 2. Formula construction
  #---------------------------
  if (is.null(formula_mixed)) {

    build_smooth <- function(var) {
      k_val <- ifelse(!is.null(k_list[[var]]), k_list[[var]], k_list$default)
      paste0("s(", var, ", k = ", k_val, ")")
    }

    smooth_terms <- sapply(smooth_vars, build_smooth)

    formula_str <- paste(
      "cloud_beta ~",
      paste(c(
        smooth_terms,
        linear_vars,
        "s(month, bs = 're')"
      ), collapse = " + ")
    )

    formula_mixed <- as.formula(formula_str)
  }

  #---------------------------
  # 3. Model fitting
  #---------------------------
  mod <- mgcv::gam(
    formula_mixed,
    data = df,
    family = mgcv::betar(link = "logit"),
    method = "REML"
  )

  #---------------------------
  # 4. Summary outputs
  #---------------------------
  summary_obj <- summary(mod)

  summary_str <- capture.output(summary_obj) %>%
    paste(collapse = "\n")

  gam_check <- capture.output(mgcv::gam.check(mod)) %>%
    paste(collapse = "\n")

  conc <- mgcv::concurvity(mod)

  #---------------------------
  # 5. Plotting
  #---------------------------
  if (output_plots) {

    grDevices::pdf(plot_file, width = 10, height = 8)

    # --- Smooth terms (4 per page) ---
    plot(mod, pages = ceiling(length(mod$smooth) / 4), shade = TRUE)

    # --- Diagnostics ---
    par(mfrow = c(2, 2))
    gam.check(mod)

    grDevices::dev.off()
  }

  #---------------------------
  # 6. Return object
  #---------------------------
  return(list(
    model = mod,
    formula = formula_mixed,
    summary = summary_str,
    gam_check = gam_check,
    concurvity = conc,
    data_used = df
  ))
}
