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
    data=arf_daily,
    response = "mean_cf_2_6",
    date_var = "date",

    # Covariate specification
    smooth_vars = c("median_PBLH", "LE_f_daytime_sum", "SW_IN_daytime_sum",
                    "evening_cbh", "LapseRate_850_500_CperKm",
                    "RH_1_7_1_daytime_mean", "G_2_1_1_daytime_sum",
                    "TA_1_2_1_daytime_mean"),

    linear_vars = c("PA_1_2_1_daytime_mean", "H_f_daytime_sum",
                    "mean_LCL_4_7", "LTS_K", "SWC_3_1_1_daytime_mean",
                    "PW_mm"),

    k_smooth = 8,
    k_map = list(median_PBLH = 12),

    # Filters
    months = 4:10,
    blc_values = c(1, 2),

    # Optional direct formula override
    formula_mixed = NULL,

    eps = 0.001,
    family = mgcv::betar(link = "logit"),
    method = "REML",

    produce_plots = TRUE
) {

  require(mgcv)
  require(dplyr)
  require(lubridate)

  df <- data

  #---------------------------
  # 1. Preprocessing
  #---------------------------
  df[[date_var]] <- as.Date(df[[date_var]])

  df <- df %>%
    mutate(
      month = lubridate::month(.data[[date_var]]),
      cloud_beta = pmin(pmax(.data[[response]], eps), 1 - eps)
    ) %>%
    filter(month %in% months) %>%
    filter(blc_flag %in% blc_values)

  # Drop NA
  all_vars <- unique(c("cloud_beta", smooth_vars, linear_vars, "month"))
  df <- df %>% tidyr::drop_na(dplyr::all_of(all_vars))

  #---------------------------
  # 2. Build formula
  #---------------------------
  if (is.null(formula_mixed)) {

    smooth_terms <- sapply(smooth_vars, function(v) {
      k_val <- ifelse(!is.null(k_map[[v]]), k_map[[v]], k_smooth)
      paste0("s(", v, ", k = ", k_val, ")")
    })

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
  # 3. Fit model
  #---------------------------
  model <- mgcv::gam(
    formula_mixed,
    data = df,
    family = family,
    method = method
  )

  #---------------------------
  # 4. Outputs
  #---------------------------
  model_summary <- capture.output(summary(model))
  model_check   <- capture.output(gam.check(model))
  model_concurv <- capture.output(concurvity(model))

  #---------------------------
  # 5. Plotting
  #---------------------------
  plot_list <- NULL

  if (produce_plots) {

    n_terms <- length(model$smooth)
    n_pages <- ceiling(n_terms / 4)

    plot_list <- vector("list", n_pages)

    for (i in seq_len(n_pages)) {
      plot_list[[i]] <- recordPlot({
        plot(
          model,
          pages = 1,
          select = ((i - 1) * 4 + 1):min(i * 4, n_terms)
        )
      })
    }
  }

  #---------------------------
  # 6. Return object
  #---------------------------
  return(list(
    model = model,
    formula = formula_mixed,
    summary = model_summary,
    gam_check = model_check,
    concurvity = model_concurv,
    plots = plot_list,
    data_used = df
  ))
}
