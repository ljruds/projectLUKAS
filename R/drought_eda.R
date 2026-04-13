




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
