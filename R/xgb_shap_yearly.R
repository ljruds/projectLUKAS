#' Train XGBoost Model and Compute Yearly SHAP Contributions
#'
#' This function trains an XGBoost regression model using a cloud fraction as a
#' response variable ('cf_2_6_mean') and set of covariates, evaluates model performance, and computes
#' SHAP (SHapley Additive exPlanations) values by year to visualize feature behaviour.
#'  SHAP values are aggregated by year to assess seasonal variability in feature importance.
#' We sort yearly data by month such that we eliminate most of the seasonal variability that would
#' otherwise exist for an entire year's worth of data, at the expense of data amount. We maintain
#' flexibility to change sampling period however needed.
#'
#' @param response Character string. Name of the response (target) variable,
#' default is cloud fraction ('cf_2_6_mean') but user has freedom to chose any response from set of covariates in data.
#' @param covariates Character vector. Names of predictor variables used in the model.
#' @param use_blc_filter Logical. Whether to filter data to specific
#'   boundary layer cloud (BLC) classes (default = FALSE).
#' @param blc_values Numeric vector. Values of \code{blc_flag} to retain
#'   when \code{use_blc_filter = TRUE} (default = c(1, 2)).
#' @param start_md String. Month and days of start of sampling period in the form
#'  "MM-DD". Default is "08-01", to encompass drought period.
#' @param end_md String. Month and days of end of sampling period in the form
#'  "MM-DD". Default is "09-30", to encompass drought period.
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
#'   \item Loads and preprocesses the dataset \code{arf_daily}, including optional filtering. Preprocessing includes partitioning into
#'   years based on start_md, end_md, and filtering for BLC days.
#'   \item Splits the data into training and testing sets.
#'   \item Trains an XGBoost regression model.
#'   \item Evaluates performance using RMSE and R-squared.
#'   \item Computes SHAP values using \code{SHAPforxgboost}.
#'   \item Aggregates SHAP values by year to quantify feature importance differences
#'   between years.
#'   \item Produces visualization outputs for performance and SHAP feature importance.
#' }
#'
#' @return A list containing:
#' \describe{
#'   \item{model}{Trained XGBoost model object.}
#'   \item{metrics}{List with RMSE, MAE, R-squared, and bias.}
#'   \item{data}{List containing train and test datasets.}
#'   \item{plots}{List of ggplot objects:
#'     \describe{
#'       \item{performance}{Observed vs predicted scatter plot.}
#'       \item{shap_bar}{Bar plot of mean absolute SHAP values by year group.}
#'       \item{shap_line}{Line plot showing evolution of SHAP importance by year.}
#'       \item{shap_summary}{SHAP summary plot.}
#'       \item{shap_dependence}{List of SHAP dependence plots for each feature.}
#'     }
#'   }
#'   \item{shap}{List containing:
#'     \describe{
#'       \item{shap_values}{Raw SHAP output from \code{shap.values}.}
#'       \item{shap_long}{Long-format SHAP data.}
#'       \item{shap_year}{Aggregated SHAP importance by year group.}
#'     }
#'   }
#' }
#'
#' @examples
#' \dontrun{
#' result <- xgb_shap_yearly(
#'   data = my_data,
#'   response = "cloud_fraction",
#'   covariates = c("temperature", "humidity", "wind_speed"),
#'   years = c(2023, 2024, 2025),
#'   start_md = "07-15",
#'   end_md = "09-25",
#'   use_blc_filter = TRUE,
#'   xgb_params = list(
#'     objective = "reg:squarederror",
#'     eval_metric = "rmse",
#'     eta = 0.01,
#'     max_depth = 6
#'   ),
#'   nrounds = 500
#' )
#'
#' result$metrics
#' result$plots$performance
#' result$plots$shap_bar
#' }
#'
#' @import dplyr
#' @import lubridate
#' @import xgboost
#' @import ggplot2
#' @import SHAPforxgboost
#'
#'
#'
#' @export
#'
#'

xgb_shap_yearly <- function(data = arf_daily,
                                    response="mean_cf_2_6",
                                    covariates,
                                    years = c(2023, 2024, 2025),
                                    start_md = "08-01",
                                    end_md   = "09-30",
                                    use_blc_filter = FALSE,
                                    blc_values = c(1, 2),
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
                                    nrounds = 1000) {

  require(dplyr)
  require(lubridate)
  require(xgboost)
  require(ggplot2)
  require(SHAPforxgboost)

  # ---------------------------
  # 1. Preprocessing + filtering
  # ---------------------------

  data$date <- as.Date(data$date)

  # Optional BLC filter
  if (use_blc_filter) {
    data <- data %>%
      dplyr::filter(blc_flag %in% blc_values)
  }

  data$md <- format(data$date, "%m-%d")

  if (start_md <= end_md) {
    date_filter <- data$md >= start_md & data$md <= end_md
  } else {
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
  # 2. Train/test split
  # ---------------------------
  set.seed(seed)
  n <- nrow(dat)
  test_idx <- sample(seq_len(n), size = floor(test_size * n))

  train_data <- dat[-test_idx, ]
  test_data  <- dat[test_idx, ]

  X_train <- as.matrix(train_data[, covariates])
  y_train <- train_data[[response]]

  X_test <- as.matrix(test_data[, covariates])
  y_test <- test_data[[response]]

  dtrain <- xgb.DMatrix(data = X_train, label = y_train)
  dtest  <- xgb.DMatrix(data = X_test, label = y_test)

  # ---------------------------
  # 3. Train model
  # ---------------------------
  model <- xgb.train(
    params = xgb_params,
    data = dtrain,
    nrounds = nrounds,
    watchlist = list(train = dtrain, test = dtest),
    verbose = 0
  )

  # ---------------------------
  # 4. Metrics
  # ---------------------------
  preds_test <- predict(model, X_test)

  test_rmse <- sqrt(mean((preds_test - y_test)^2))
  test_mae  <- mean(abs(preds_test - y_test))
  test_r2   <- 1 - sum((preds_test - y_test)^2) / sum((y_test - mean(y_test))^2)
  test_bias <- mean(preds_test - y_test)

  perf_df <- data.frame(observed = y_test, predicted = preds_test)

  perf_plot <- ggplot(perf_df, aes(observed, predicted)) +
    geom_point(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    theme_minimal() +
    labs(title = "XGBoost Performance", x = "Observed", y = "Predicted")

  # ---------------------------
  # 5. SHAP
  # ---------------------------
  shap_values <- shap.values(model, X_train)

  shap_long <- shap.prep(
    shap_contrib = shap_values$shap_score,
    X_train = X_train
  )

  shap_long$year <- train_data$year[shap_long$ID]

  shap_year <- shap_long %>%
    mutate(year_group = ifelse(year %in% c(2023, 2024), "2023–2024", "2025")) %>%
    group_by(year_group, variable) %>%
    summarise(mean_abs_shap = mean(abs(value)), .groups = "drop")

  shap_bar_plot <- ggplot(shap_year,
                          aes(reorder(variable, mean_abs_shap),
                              mean_abs_shap,
                              fill = year_group)) +
    geom_col(position = "dodge") +
    coord_flip() +
    theme_minimal()

  shap_time <- shap_long %>%
    group_by(year, variable) %>%
    summarise(mean_abs_shap = mean(abs(value)), .groups = "drop")

  shap_line_plot <- ggplot(shap_time,
                           aes(year, mean_abs_shap,
                               color = variable,
                               group = variable)) +
    geom_line() +
    geom_point() +
    theme_minimal()

  shap_summary_plot <- shap.plot.summary(shap_long)

  create_shap_dependence <- function(var) {
    shap.plot.dependence(shap_long, x = var, smooth = TRUE)
  }

  shap_dependence_plots <- lapply(covariates, create_shap_dependence)
  names(shap_dependence_plots) <- covariates

  # ---------------------------
  # 6. Output
  # ---------------------------
  return(list(
    # --- model + data ---
    model = model,

    data = list(
      train = train_data,
      test  = test_data
    ),

    # --- metrics ---
    metrics = list(
      RMSE = test_rmse,
      MAE  = test_mae,
      R2   = test_r2,
      Bias = test_bias
    ),

    # --- plots ---
    performance_plot = perf_plot,
    shap_bar = shap_bar_plot,
    shap_line = shap_line_plot,
    shap_summary = shap_summary_plot,
    shap_dependence = shap_dependence_plots,

    # --- shap outputs ---
    shap_values = shap_values,
    shap_long = shap_long,
    shap_year = shap_year
  ))
}
