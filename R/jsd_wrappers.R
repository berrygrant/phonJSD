#' Estimate Jensen–Shannon divergence or distance between two categories
#'
#' @param data Data frame with at least `category_col` and `features`.
#' @param features Character vector of feature column names (e.g., c("F1","F2")).
#' @param category_col Name of the column giving the two-way category factor.
#' @param group_col Optional grouping column (not yet implemented in this simple version).
#' @param do_boot Logical; if TRUE, run nonparametric bootstrap.
#' @param n_boot Number of bootstrap resamples.
#' @param min_tokens Minimum tokens per category (currently not enforced here).
#' @param est_distance Logical; if TRUE, return Jensen–Shannon *distance* (sqrt of divergence).
#' @param conf_level Confidence level for bootstrap interval.
#'
#' @return A tibble with one row and columns:
#'   scope, n_tokens, n_boot, jsd_point, jsd_mean, jsd_sd, jsd_low, jsd_high.
#'
#' @export
estimate_jsd <- function(data,
                         features,
                         category_col,
                         group_col    = NULL,
                         do_boot      = FALSE,
                         n_boot       = 1000,
                         min_tokens   = 5,
                         est_distance = FALSE,
                         conf_level   = 0.95) {
  
  # For now, only global (no group_col support in this implementation)
  if (!is.null(group_col)) {
    stop("estimate_jsd(): `group_col` support not yet implemented in this version.")
  }
  
  # Basic checks
  if (!all(features %in% names(data))) {
    stop("estimate_jsd(): All `features` must be columns in `data`.")
  }
  if (!category_col %in% names(data)) {
    stop("estimate_jsd(): `category_col` must be a column in `data`.")
  }
  
  # Subset and drop NAs
  df <- data %>%
    dplyr::select(dplyr::all_of(c(category_col, features))) %>%
    tidyr::drop_na()
  
  # Ensure exactly two categories
  g <- droplevels(factor(df[[category_col]]))
  if (nlevels(g) != 2L) {
    stop("estimate_jsd(): `category_col` must have exactly 2 levels in the filtered data.")
  }
  
  # ---- Point estimate via KDE ----
  jsd_div_point <- jsd_kde_nd(
    data     = df,
    features = features,
    group    = category_col  # pass column name, not vector
  )
  
  jsd_point <- if (est_distance) sqrt(jsd_div_point) else jsd_div_point
  
  # If no bootstrap, just return the point estimate row
  if (!do_boot) {
    return(tibble::tibble(
      scope     = "global",
      n_tokens  = nrow(df),
      n_boot    = 0L,
      jsd_point = jsd_point,
      jsd_mean  = NA_real_,
      jsd_sd    = NA_real_,
      jsd_low   = NA_real_,
      jsd_high  = NA_real_
    ))
  }
  
  # ---- Bootstrap case ----
  alpha <- 1 - conf_level
  
  boot_vals <- replicate(n_boot, {
    idx <- sample.int(nrow(df), replace = TRUE)
    df_boot <- df[idx, , drop = FALSE]
    
    jsd_div_boot <- jsd_kde_nd(
      data     = df_boot,
      features = features,
      group    = category_col  # ← THIS is the key fix
    )
    
    if (est_distance) sqrt(jsd_div_boot) else jsd_div_boot
  })
  
  jsd_mean <- mean(boot_vals)
  jsd_sd   <- stats::sd(boot_vals)
  qs       <- stats::quantile(
    boot_vals,
    probs = c(alpha / 2, 1 - alpha / 2),
    names = FALSE
  )
  
  tibble::tibble(
    scope     = "global",
    n_tokens  = nrow(df),
    n_boot    = as.integer(n_boot),
    jsd_point = jsd_point,
    jsd_mean  = jsd_mean,
    jsd_sd    = jsd_sd,
    jsd_low   = qs[1],
    jsd_high  = qs[2]
  )
}