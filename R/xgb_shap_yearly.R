

run_xgb_yearly_analysis <- function(data,
                                    response,
                                    covariates,
                                    years = c(2023, 2024, 2025),
                                    start_md = "08-01",
                                    end_md   = "09-30",
                                    test_frac = 0.2,
                                    seed = 42) {

  require(dplyr)
  require(lubridate)
  require(xgboost)
  require(ggplot2)
  require(SHAPforxgboost)

  # ---------------------------
  # 1. Data filtering
  # ---------------------------
  data$date <- as.Date(data$date)

  # Extract month-day string for filtering
  data$md <- format(data$date, "%m-%d")

  # Handle ranges that DO NOT cross year boundary (typical case)
  if (start_md <= end_md) {
    date_filter <- data$md >= start_md & data$md <= end_md
  } else {
    # Handles wrap-around ranges (e.g., Nov–Feb)
    date_filter <- data$md >= start_md | data$md <= end_md
  }

  dat <- data %>%
    dplyr::filter(
      lubridate::year(date) %in% years,
      date_filter
    ) %>%
    dplyr::select(date, dplyr::all_of(response), dplyr::all_of(covariates)) %>%
    na.omit()

  dat$year <- lubridate::year(dat$date)

  # ---------------------------
  # 2. Train / Test split
  # ---------------------------
  set.seed(seed)
  n <- nrow(dat)
  test_idx <- sample(seq_len(n), size = floor(test_frac * n))

  train_data <- dat[-test_idx, ]
  test_data  <- dat[test_idx, ]

  X_train <- as.matrix(train_data[, covariates])
  y_train <- train_data[[response]]

  X_test <- as.matrix(test_data[, covariates])
  y_test <- test_data[[response]]

  dtrain <- xgb.DMatrix(data = X_train, label = y_train)

  # ---------------------------
  # 3. Train XGBoost model
  # ---------------------------
  model <- xgb.train(
    params = list(
      objective = "reg:squarederror",
      eval_metric = "rmse",
      eta = 0.05,
      max_depth = 6,
      subsample = 0.8,
      colsample_bytree = 0.8
    ),
    data = dtrain,
    nrounds = 1000,
    verbose = 0
  )

  # ---------------------------
  # 4. Predictions + metrics
  # ---------------------------
  preds_test <- predict(model, X_test)

  test_rmse <- sqrt(mean((preds_test - y_test)^2))
  test_mae  <- mean(abs(preds_test - y_test))
  test_r2   <- 1 - sum((preds_test - y_test)^2) / sum((y_test - mean(y_test))^2)
  test_bias <- mean(preds_test - y_test)

  perf_df <- data.frame(
    observed = y_test,
    predicted = preds_test
  )

  # ---------------------------
  # 5. Performance plot
  # ---------------------------
  perf_plot <- ggplot(perf_df, aes(observed, predicted)) +
    geom_point(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    theme_minimal() +
    labs(
      title = "XGBoost Performance",
      x = "Observed",
      y = "Predicted"
    )

  # ---------------------------
  # 6. SHAP values
  # ---------------------------
  shap_values <- shap.values(
    xgb_model = model,
    X_train = X_train
  )

  shap_long <- shap.prep(
    shap_contrib = shap_values$shap_score,
    X_train = X_train
  )

  shap_long$year <- train_data$year[shap_long$ID]

  # ---------------------------
  # 7. SHAP aggregation
  # ---------------------------
  shap_year <- shap_long %>%
    mutate(year_group = ifelse(year %in% c(2023, 2024), "2023–2024", "2025")) %>%
    group_by(year_group, variable) %>%
    summarise(
      mean_abs_shap = mean(abs(value)),
      .groups = "drop"
    )

  shap_bar_plot <- ggplot(shap_year,
                          aes(reorder(variable, mean_abs_shap),
                              mean_abs_shap,
                              fill = year_group)) +
    geom_col(position = "dodge") +
    coord_flip() +
    theme_minimal() +
    labs(
      x = "Feature",
      y = "Mean |SHAP value|",
      fill = "Period",
      title = paste0("Drivers (", start_md, " to ", end_md, ")")
    )

  # ---------------------------
  # 8. SHAP evolution plot
  # ---------------------------
  shap_time <- shap_long %>%
    group_by(year, variable) %>%
    summarise(mean_abs_shap = mean(abs(value)), .groups = "drop")

  shap_line_plot <- ggplot(shap_time,
                           aes(year, mean_abs_shap,
                               color = variable,
                               group = variable)) +
    geom_line() +
    geom_point(size = 2) +
    theme_minimal() +
    labs(
      y = "Mean |SHAP value|",
      title = "Evolution of Feature Importance by Year"
    )

  # ---------------------------
  # 9. SHAP summary
  # ---------------------------
  shap_summary_plot <- shap.plot.summary(shap_long)

  # ---------------------------
  # 10. SHAP dependence
  # ---------------------------
  create_shap_dependence <- function(var) {
    shap.plot.dependence(
      data_long = shap_long,
      x = var,
      smooth = TRUE
    )
  }

  shap_dependence_plots <- lapply(covariates, create_shap_dependence)
  names(shap_dependence_plots) <- covariates

  # ---------------------------
  # 11. Output
  # ---------------------------
  return(list(
    model = model,
    metrics = list(
      RMSE = test_rmse,
      MAE  = test_mae,
      R2   = test_r2,
      Bias = test_bias
    ),
    data = list(
      train = train_data,
      test  = test_data
    ),
    plots = list(
      performance = perf_plot,
      shap_bar = shap_bar_plot,
      shap_line = shap_line_plot,
      shap_summary = shap_summary_plot,
      shap_dependence = shap_dependence_plots
    ),
    shap = list(
      shap_long = shap_long,
      shap_year = shap_year
    )
  ))
}
