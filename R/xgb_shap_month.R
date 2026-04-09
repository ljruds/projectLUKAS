xgb_shap_monthly <- function(
    response,
    covariates,
    filter_expr = NULL,   # e.g. quote(blc_flag %in% c(1,2))
    test_size = 0.2,
    seed = 42,
    xgb_params = list(
      objective = "reg:squarederror",
      eval_metric = "rmse",
      eta = 0.01,
      max_depth = 6,
      subsample = 0.8,
      colsample_bytree = 0.8
    ),
    nrounds = 1000
) {

  library(dplyr)
  library(lubridate)
  library(xgboost)
  library(caret)
  library(ggplot2)
  library(SHAPforxgboost)

  #---------------------------
  # 1. Load + preprocess
  #---------------------------
  data <- arf_daily
  data$date <- as.Date(data$date)

  if (!is.null(filter_expr)) {
    data <- data %>% filter(!!filter_expr)
  }

  dat <- data %>%
    select(date, all_of(response), all_of(covariates)) %>%
    na.omit() %>%
    mutate(month = lubridate::month(date, label = TRUE))

  #---------------------------
  # 2. Matrix + split
  #---------------------------
  X <- as.matrix(dat[, covariates])
  y <- dat[[response]]

  set.seed(seed)
  train_index <- createDataPartition(y, p = 1 - test_size, list = FALSE)

  X_train <- X[train_index, ]
  X_test  <- X[-train_index, ]

  y_train <- y[train_index]
  y_test  <- y[-train_index]

  dtrain <- xgb.DMatrix(X_train, label = y_train)
  dtest  <- xgb.DMatrix(X_test, label = y_test)

  #---------------------------
  # 3. Train model
  #---------------------------
  model <- xgb.train(
    params = xgb_params,
    data = dtrain,
    nrounds = nrounds,
    watchlist = list(train = dtrain, test = dtest),
    verbose = 0
  )

  #---------------------------
  # 4. Performance
  #---------------------------
  preds <- predict(model, X_test)

  rmse <- sqrt(mean((preds - y_test)^2))
  r2   <- cor(preds, y_test)^2

  perf_df <- data.frame(observed = y_test, predicted = preds)

  perf_plot <- ggplot(perf_df, aes(observed, predicted)) +
    geom_point(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    theme_minimal() +
    labs(
      title = "XGBoost Performance",
      x = "Observed",
      y = "Predicted"
    )

  #---------------------------
  # 5. SHAP values
  #---------------------------
  shap <- shap.values(
    xgb_model = model,
    X_train = X
  )

  shap_long <- shap.prep(
    shap_contrib = shap$shap_score,
    X_train = X
  )

  # attach month correctly
  shap_long$month <- dat$month[shap_long$ID]

  # aggregate monthly SHAP importance
  shap_month <- shap_long %>%
    group_by(month, variable) %>%
    summarise(mean_abs_shap = mean(abs(value)), .groups = "drop")

  #---------------------------
  # 6. SHAP plots
  #---------------------------
  shap_bar <- ggplot(shap_month,
                     aes(month, mean_abs_shap, fill = variable)) +
    geom_col() +
    theme_minimal() +
    labs(
      title = "Monthly Drivers (SHAP)",
      y = "Mean |SHAP|",
      x = "Month"
    )

  shap_line <- ggplot(shap_month,
                      aes(month, mean_abs_shap,
                          color = variable,
                          group = variable)) +
    geom_line(size = 1.2) +
    geom_point(size = 2) +
    theme_minimal() +
    labs(
      title = "Seasonal SHAP Evolution",
      y = "Mean |SHAP|"
    )

  shap_summary <- shap.plot.summary(shap_long)

  #---------------------------
  # 7. Return everything
  #---------------------------
  return(list(
    model = model,
    metrics = list(RMSE = rmse, R2 = r2),
    performance_plot = perf_plot,
    shap = shap,
    shap_long = shap_long,
    shap_month = shap_month,
    shap_bar = shap_bar,
    shap_line = shap_line,
    shap_summary = shap_summary
  ))
}
