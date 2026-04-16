#' Fit GAMM for Boundary Layer Cloud Fraction (Beta Regression)
#'
#' @param data Data frame containing predictors and response
#' @param response Character, response variable (default: "mean_cf_2_6")
#' @param covariates Character vector of covariate names
#' @param months Integer vector of months to retain (default: 4:10)
#' @param blc_values Values of blc_flag to retain (default: c(1,2))
#' @param k_list Named list of k values for smooth terms
#' @param produce_plots Logical, whether to generate plots
#' @param eps Small value for beta regression bounds
#'
#' @return List containing model, summary text, diagnostics, plots
#' @export

fit_cloud_gamm <- function(
    data = arf_daily,
    response = "mean_cf_2_6",
    covariates,
    months = 4:10,
    blc_values = c(1, 2),
    k_list = list(
      median_PBLH = 12,
      LE_f_daytime_sum = 8,
      SW_IN_daytime_sum = 8,
      evening_cbh = 8,
      LapseRate_850_500_CperKm = 8,
      RH_1_7_1_daytime_mean = 8,
      G_2_1_1_daytime_sum = 8,
      TA_1_2_1_daytime_mean = 8
    ),
    produce_plots = TRUE,
    eps = 0.001
) {

  requireNamespace("mgcv")
  requireNamespace("dplyr")
  requireNamespace("ggplot2")

  library(mgcv)
  library(dplyr)
  library(ggplot2)

  #---------------------------
  # 1. Data preprocessing
  #---------------------------
  df <- data

  if (!("date" %in% names(df))) {
    stop("Data must contain a 'date' column.")
  }

  df$date <- as.Date(df$date)

  df <- df %>%
    mutate(
      month = lubridate::month(date),
      cloud_beta = pmin(pmax(.data[[response]], eps), 1 - eps)
    ) %>%
    filter(month %in% months) %>%
    filter(blc_flag %in% blc_values) %>%
    tidyr::drop_na(all_of(c("cloud_beta", covariates, "month")))

  #---------------------------
  # 2. Build formula dynamically
  #---------------------------
  smooth_term <- function(var) {
    k <- k_list[[var]]
    if (is.null(k)) k <- 8
    paste0("s(", var, ", k=", k, ")")
  }

  smooth_vars <- c(
    "median_PBLH",
    "LE_f_daytime_sum",
    "SW_IN_daytime_sum",
    "evening_cbh",
    "LapseRate_850_500_CperKm",
    "RH_1_7_1_daytime_mean",
    "G_2_1_1_daytime_sum",
    "TA_1_2_1_daytime_mean"
  )

  linear_vars <- c(
    "PW_mm",
    "PA_1_2_1_daytime_mean",
    "H_f_daytime_sum",
    "mean_LCL_4_7",
    "LTS_K",
    "SWC_3_1_1_daytime_mean"
  )

  formula_str <- paste(
    "cloud_beta ~",
    paste(c(
      sapply(smooth_vars, smooth_term),
      linear_vars,
      "s(month, bs='re')"
    ), collapse = " + ")
  )

  formula_mixed <- as.formula(formula_str)

  #---------------------------
  # 3. Fit model
  #---------------------------
  model <- gam(
    formula_mixed,
    data = df,
    family = betar(link = "logit"),
    method = "REML"
  )

  #---------------------------
  # 4. Summary outputs
  #---------------------------
  model_summary <- capture.output(summary(model))
  gam_check <- capture.output(gam.check(model))
  concurv <- capture.output(concurvity(model, full = TRUE))

  #---------------------------
  # 5. Plotting
  #---------------------------
  plot_list <- NULL

  if (produce_plots) {

    # mgcv base plots → captured as list
    plot_list <- list()

    pdf(NULL)  # prevent file output
    plot(model, pages = 1, residuals = TRUE)
    dev.off()

    # Optional: ggplot partial effects (cleaner)
    plot_smooth <- function(var) {
      vis.gam(
        model,
        view = var,
        plot.type = "response",
        main = paste("Effect of", var)
      )
    }

    smooth_plots <- lapply(smooth_vars, function(v) {
      try(plot_smooth(v), silent = TRUE)
    })

    plot_list$smooths <- smooth_plots
  }

  #---------------------------
  # 6. Return structured object
  #---------------------------
  return(list(
    model = model,
    formula = formula_mixed,
    summary = paste(model_summary, collapse = "\n"),
    gam_check = paste(gam_check, collapse = "\n"),
    concurvity = paste(concurv, collapse = "\n"),
    plots = plot_list,
    data_used = df
  ))
}
