
#' Drought vs Non-Drought Distribution Analysis
#'
#' Performs comparative distribution analysis between drought and non-drought periods
#' for a set of covariates. The function generates faceted density and boxplots,
#' and computes nonparametric statistical tests and effect sizes.
#'
#' @description
#' This function partitions a dataset into drought and non-drought periods based on
#' user-defined date thresholds. For each specified covariate, it evaluates whether
#' the distributions differ between periods using:
#'
#' - Wilcoxon rank-sum test (Mann–Whitney U equivalent)
#' - Kolmogorov–Smirnov test
#' - Rank-biserial correlation
#' - Cliff’s Delta
#'
#' It also produces faceted visualizations to compare empirical distributions.
#'
#' @param data A data.frame containing the dataset.
#' @param covariates Character vector of variable names to analyze.
#' @param date_col Name of the date column (default: "date").
#' @param drought_start Start date of drought period (character or POSIXct).
#' @param drought_end End date of drought period (character or POSIXct).
#' @param lower_bound Optional numeric lower bound for variable clipping.
#' @param upper_bound Optional numeric upper bound for variable clipping.
#'
#' @details
#' ## Data Partitioning
#' Observations are classified into two groups:
#'
#' \deqn{
#' \text{period} =
#' \begin{cases}
#' \text{Drought}, & t \in [t_{start}, t_{end}] \\
#' \text{Non-drought}, & \text{otherwise}
#' \end{cases}
#' }
#'
#' ## Statistical Tests
#'
#' ### 1. Wilcoxon Rank-Sum Test (Mann–Whitney U)
#' Tests for differences in central tendency between two independent samples.
#'
#' Relationship between reported statistic \eqn{W} and Mann–Whitney \eqn{U}:
#'
#' \deqn{
#' U = W - \frac{n_1(n_1 + 1)}{2}
#' }
#'
#' where:
#' - \eqn{n_1}: number of drought observations
#' - \eqn{n_2}: number of non-drought observations
#'
#' ### 2. Kolmogorov–Smirnov Test
#' Measures the maximum difference between empirical cumulative distribution functions:
#'
#' \deqn{
#' D = \sup_x |F_1(x) - F_2(x)|
#' }
#'
#' where:
#' - \eqn{F_1(x)}, \eqn{F_2(x)} are empirical CDFs
#'
#' ## Effect Sizes
#'
#' ### 1. Rank-Biserial Correlation
#'
#' \deqn{
#' r = 1 - \frac{2U}{n_1 n_2}
#' }
#'
#' Interpretation:
#' - \eqn{r \approx 0}: little separation
#' - \eqn{r > 0}: drought values tend to be larger
#' - \eqn{r < 0}: non-drought values tend to be larger
#'
#' ### 2. Cliff’s Delta
#'
#' \deqn{
#' \delta = P(X > Y) - P(X < Y)
#' }
#'
#' where:
#' - \eqn{X}: drought sample
#' - \eqn{Y}: non-drought sample
#'
#' Interpretation:
#' - \eqn{\delta = 0}: no difference
#' - \eqn{\delta = \pm1}: complete separation
#'
#' ## Visualization
#'
#' - Kernel density estimates (KDE) are used to approximate probability densities.
#' - Boxplots summarize median, interquartile range, and dispersion.
#' - Faceting allows comparison across multiple variables simultaneously.
#'
#' ## Notes
#'
#' - KDE estimation is sensitive to boundary effects when variables are bounded.
#' - Optional clipping ensures values lie within physically meaningful limits.
#'
#' @return A list with:
#' \describe{
#'   \item{plots}{
#'     \itemize{
#'       \item density: Faceted density plot
#'       \item box: Faceted boxplot
#'     }
#'   }
#'   \item{stats}{
#'     Data frame containing:
#'     \itemize{
#'       \item variable name
#'       \item sample sizes
#'       \item test statistics
#'       \item p-values
#'       \item effect sizes
#'     }
#'   }
#' }
#'
#' @examples
#' \dontrun{
#' result <- drought_density_analysis(
#'   data = df,
#'   covariates = c("mean_cf_2_6", "LCL"),
#'   drought_start = "2025-08-10",
#'   drought_end = "2025-09-25",
#'   lower_bound = 0,
#'   upper_bound = 1
#' )
#'
#' result$plots$density
#' result$stats
#' }
#'
#' @section Interpretation Guide:
#'
#' - [Add domain-specific interpretation here]
#' - [Explain expected behavior of variables during drought]
#' - [Relate statistical significance to physical mechanisms]
#'
#' @section Limitations:
#'
#' - [Discuss sampling biases]
#' - [Discuss temporal autocorrelation]
#' - [Discuss KDE boundary effects if relevant]
#'
#' @export



drought_eda <- function(data=arf_daily,
                                     covariates,
                                     date_col = "date",
                                     drought_start,
                                     drought_end,
                                     lower_bound = NULL,
                                     upper_bound = NULL) {

  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(purrr)
  library(effsize)

  # --- Ensure datetime ---
  data[[date_col]] <- as.POSIXct(data[[date_col]])
  drought_start <- as.POSIXct(drought_start)
  drought_end   <- as.POSIXct(drought_end)

  # --- Label periods ---
  data <- data %>%
    mutate(period = ifelse(
      .data[[date_col]] >= drought_start & .data[[date_col]] <= drought_end,
      "Drought", "Non-drought"
    ))

  # --- Pivot to long format for faceting ---
  df_long <- data %>%
    select(all_of(c(date_col, covariates, "period"))) %>%
    pivot_longer(cols = all_of(covariates),
                 names_to = "variable",
                 values_to = "value") %>%
    drop_na()

  # --- Apply bounds if provided ---
  if (!is.null(lower_bound) | !is.null(upper_bound)) {
    df_long <- df_long %>%
      mutate(value = pmin(pmax(value, lower_bound %||% -Inf),
                          upper_bound %||% Inf))
  }

  # --- Split helper ---
  split_data <- df_long %>%
    group_by(variable) %>%
    group_split()

  # --- Stats computation per variable ---
  stats_list <- map(split_data, function(df_var) {

    var_name <- unique(df_var$variable)

    drought_vals <- df_var %>%
      filter(period == "Drought") %>%
      pull(value)

    nondrought_vals <- df_var %>%
      filter(period == "Non-drought") %>%
      pull(value)

    n1 <- length(drought_vals)
    n2 <- length(nondrought_vals)

    # --- Wilcoxon ---
    w_test <- wilcox.test(drought_vals, nondrought_vals, exact = FALSE)
    W <- as.numeric(w_test$statistic)
    U <- W - n1 * (n1 + 1) / 2

    # --- Rank-biserial ---
    r_rb <- 1 - (2 * U) / (n1 * n2)

    # --- KS test ---
    ks <- ks.test(drought_vals, nondrought_vals)

    # --- Cliff's delta ---
    cliff <- effsize::cliff.delta(drought_vals, nondrought_vals)

    data.frame(
      variable = var_name,
      n_drought = n1,
      n_nondrought = n2,
      U = U,
      p_value_wilcox = w_test$p.value,
      p_value_ks = ks$p.value,
      ks_statistic = as.numeric(ks$statistic),
      rank_biserial_r = r_rb,
      cliffs_delta = cliff$estimate
    )
  })

  stats_table <- bind_rows(stats_list)

  # --- Faceted density plot ---
  density_plot <- ggplot(df_long, aes(x = value, fill = period)) +
    geom_density(alpha = 0.3) +
    facet_wrap(~variable, scales = "free") +
    labs(
      title = "Density: Drought vs Non-drought",
      x = "Value",
      y = "Density"
    ) +
    theme_minimal()

  if (!is.null(lower_bound) & !is.null(upper_bound)) {
    density_plot <- density_plot +
      xlim(lower_bound, upper_bound) +
      geom_vline(xintercept = c(lower_bound, upper_bound),
                 linetype = "dashed")
  }

  # --- Faceted boxplot ---
  box_plot <- ggplot(df_long, aes(x = period, y = value)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.25) +
    facet_wrap(~variable, scales = "free") +
    labs(
      title = "Boxplots: Drought vs Non-drought",
      x = "",
      y = "Value"
    ) +
    theme_minimal()

  if (!is.null(lower_bound) & !is.null(upper_bound)) {
    box_plot <- box_plot + ylim(lower_bound, upper_bound)
  }

  return(list(
    plots = list(
      density = density_plot,
      box = box_plot
    ),
    stats = stats_table
  ))
}
