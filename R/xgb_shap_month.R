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
#' @param use_blc_filter Logical. Whether to filter data to specific
#'   boundary layer cloud (BLC) classes (default = TRUE).
#' @param blc_values Numeric vector. Values of \code{blc_flag} to retain
#'   when \code{use_blc_filter = TRUE} (default = c(1, 2)).
#' @param months Integer vector. Months to include in the analysis
#'   (default = 4:8, i.e., April–August).
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
#'


xgb_shap_monthly <- function(
    response="mean_cf_2_6",
    covariates,
    use_blc_filter = TRUE,
    blc_values = c(1, 2),
    months = 4:10,
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
    nrounds = 1000,
    produce_plots = FALSE,
    top_n = 8
) {

  library(dplyr)
  library(lubridate)
  library(xgboost)
  library(caret)
  library(ggplot2)
  library(SHAPforxgboost)
  library(tidyr)

  #---------------------------
  # 1. Load + preprocess
  #---------------------------
  data <- arf_daily
  data$date <- as.Date(data$date)

  # --- Apply filters FIRST (like your QMD) ---

  # Month filter (default April–August)
  data <- data %>%
    dplyr::filter(lubridate::month(date) %in% months)

  # Optional BLC filter
  if (use_blc_filter) {
    data <- data %>%
      dplyr::filter(blc_flag %in% blc_values)
  }

  # --- THEN build modeling dataset (clean pipeline) ---
  dat <- data %>%
    dplyr::select(date, dplyr::all_of(response), dplyr::all_of(covariates)) %>%
    stats::na.omit() %>%
    dplyr::mutate(
      month = lubridate::month(date, label = TRUE),
      ID = dplyr::row_number()
    )

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
    theme_minimal()

  #---------------------------
  # 5. SHAP values (single computation)
  #---------------------------
  shap_values <- shap.values(
    xgb_model = model,
    X_train = X
  )

  shap_long <- shap.prep(
    shap_contrib = shap_values$shap_score,
    X_train = X
  )

  # Join instead of indexing (this is the fix)
  shap_long <- shap_long %>%
    dplyr::left_join(
      dat[, c("ID", "month")],
      by = "ID"
    )

  dat$ID <- seq_len(nrow(dat))

  #---------------------------
  # 6. Monthly SHAP (original method)
  #---------------------------
  shap_month <- shap_long %>%
    group_by(month, variable) %>%
    summarise(mean_abs_shap = mean(abs(value)), .groups = "drop")

  shap_bar <- ggplot(shap_month,
                     aes(month, mean_abs_shap, fill = variable)) +
    geom_col() +
    theme_minimal()

  shap_line <- ggplot(shap_month,
                      aes(month, mean_abs_shap,
                          color = variable,
                          group = variable)) +
    geom_line() +
    geom_point() +
    theme_minimal()

  #---------------------------
  # 7. NEW: Top-variable monthly SHAP workflow
  #---------------------------
  shap_df <- as.data.frame(shap_values$shap_score)

  shap_df$date  <- dat$date
  shap_df$month <- dat$month
  shap_df$doy   <- dat$doy

  top_vars <- names(sort(shap_values$mean_shap_score,
                         decreasing = TRUE))[1:top_n]

  shap_top <- shap_df[, c("month", top_vars)]

  monthly_shap <- shap_top %>%
    group_by(month) %>%
    summarise(across(all_of(top_vars), mean, na.rm = TRUE),
              .groups = "drop")

  monthly_long <- monthly_shap %>%
    pivot_longer(
      cols = -month,
      names_to = "feature",
      values_to = "shap_value"
    )

  monthly_long$month <- factor(
    monthly_long$month,
    levels = 1:12,
    labels = month.abb
  )

  shap_top_line <- ggplot(monthly_long,
                          aes(month, shap_value,
                              colour = feature,
                              group = feature)) +
    geom_line() +
    geom_point() +
    theme_minimal()

  shap_top_bar <- ggplot(monthly_long,
                         aes(month, shap_value, fill = feature)) +
    geom_col(position = "stack") +
    theme_minimal()

  #---------------------------
  # 8. OPTIONAL: dependence plots
  #---------------------------
  shap_dependence <- NULL

  if (produce_plots) {

    create_shap_dependence <- function(var) {
      shap.plot.dependence(
        data_long = shap_long,
        x = var,
        smooth = TRUE
      )
    }

    shap_dependence <- lapply(covariates, create_shap_dependence)
    names(shap_dependence) <- covariates
  }

  shap_summary <- shap.plot.summary(shap_long)

  #---------------------------
  # 9. Return
  #---------------------------
  return(list(
    model = model,
    metrics = list(RMSE = rmse, R2 = r2),
    performance_plot = perf_plot,

    shap_values = shap_values,
    shap_long = shap_long,

    shap_month = shap_month,
    shap_bar = shap_bar,
    shap_line = shap_line,

    top_variables = top_vars,
    shap_top_line = shap_top_line,
    shap_top_bar = shap_top_bar,

    shap_dependence = shap_dependence,
    shap_summary = shap_summary
  ))
}
