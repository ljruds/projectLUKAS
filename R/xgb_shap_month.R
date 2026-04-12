#' Train XGBoost Model and Compute Monthly SHAP Contributions
#'
#' This function trains an XGBoost regression model using a cloud fraction as a
#' response variable ('cf_2_6_mean') and set of covariates, evaluates model performance, and computes
#' SHAP (SHapley Additive exPlanations) values by month to visualize feature behaviour.
#'  SHAP values are aggregated by month to assess seasonal variability in feature importance.
#'
#' @param response Character string. Name of the response (target) variable,
#' should be cloud fraction ('cf_2_6_mean') but user has freedom to chose any response from set of covariates in data.
#' @param covariates Character vector. Names of predictor variables used in the model.
#' @param filter_expr Optional quoted expression used to filter the dataset
#'   (e.g., \code{quote(blc_flag \%in\% c(1,2))}). Default is \code{NULL}.
#' @param test_size Numeric. Proportion of data to use for testing (default = 0.2).
#' @param seed Integer. Random seed for reproducibility (default = 42).
#' @param xgb_params List. Parameters passed to \code{xgboost::xgb.train}.
#'   Defaults include:
#'   \itemize{
#'     \item \code{objective = "reg:squarederror"}
#'     \item \code{eval_metric = "rmse"}
#'     \item \code{eta = 0.01}
#'     \item \code{max_depth = 6}
#'     \item \code{subsample = 0.8}
#'     \item \code{colsample_bytree = 0.8}
#'   }
#' @param nrounds Integer. Number of boosting iterations (default = 1000).
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Loads and preprocesses the dataset \code{arf_daily}, including optional filtering. Preprocessing includes partitioning into months, and filtering for BLC days.
#'   \item Splits the data into training and testing sets.
#'   \item Trains an XGBoost regression model.
#'   \item Evaluates performance using RMSE and R-squared.
#'   \item Computes SHAP values using \code{SHAPforxgboost}.
#'   \item Aggregates SHAP values by month to quantify seasonal feature importance.
#'   \item Produces visualization outputs for performance and SHAP diagnostics.
#' }
#'
#' @return A list containing:
#' \describe{
#'   \item{model}{Trained XGBoost model object.}
#'   \item{metrics}{List with RMSE and R-squared values.}
#'   \item{performance_plot}{\code{ggplot2} scatter plot of observed vs predicted values.}
#'   \item{shap}{Raw SHAP output from \code{shap.values}.}
#'   \item{shap_long}{Long-format SHAP data used for plotting.}
#'   \item{shap_month}{Data frame of monthly aggregated SHAP importance.}
#'   \item{shap_bar}{Bar plot of mean absolute SHAP values by month.}
#'   \item{shap_line}{Line plot showing seasonal evolution of SHAP importance.}
#'   \item{shap_summary}{SHAP summary plot.}
#' }
#'
#' @examples
#' \dontrun{
#' result <- xgb_shap_monthly(
#'   response = "cloud_fraction",
#'   covariates = c("temperature", "humidity", "wind_speed"),
#'   filter_expr = quote(month(date) %in% 6:8)
#' )
#'
#' result$metrics
#' result$performance_plot
#' result$shap_bar
#' }
#'
#' @import dplyr
#' @import lubridate
#' @import xgboost
#' @import caret
#' @import ggplot2
#' @import SHAPforxgboost
#'
#' @export
#'
#' fish love structure


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
