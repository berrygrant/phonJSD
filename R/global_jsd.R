#' Global JSD with bootstrap confidence interval
#'
#' Computes a single Jensen–Shannon divergence (JSD) value for two categories
#' in an n-dimensional acoustic space, together with bootstrap-based
#' confidence intervals obtained by resampling tokens with replacement.
#'
#' This is the "group-wise" version of JSD: it ignores speakers and treats
#' all tokens as coming from a single population for each category.
#'
#' @param data Data frame containing at least the category column and
#'   the feature columns.
#' @param features Character vector of column names giving the acoustic
#'   dimensions (e.g., c("f1", "f2") or paste0("mfcc", 1:13)).
#' @param category_col String; name of the column giving the two categories
#'   to compare (e.g., "vowel"). Must have exactly two unique values.
#' @param n_boot Integer; number of bootstrap resamples.
#' @param min_tokens Minimum total number of non-missing tokens required.
#' @param ... Additional arguments passed to \code{jsd_kde_nd()}.
#'
#' @return A one-row data frame with columns:
#'   \itemize{
#'     \item \code{n_tokens} – total number of tokens used
#'     \item \code{n_boot} – number of successful bootstrap samples
#'     \item \code{jsd_point} – JSD on the full dataset
#'     \item \code{jsd_mean} – mean JSD across bootstrap samples
#'     \item \code{jsd_sd} – standard deviation of bootstrap JSD
#'     \item \code{jsd_low}, \code{jsd_high} – 95% bootstrap CI (2.5%, 97.5%)
#'   }
#' @export
#' @importFrom stats sd quantile
global_boot_jsd <- function(data,
                            features,
                            category_col,
                            n_boot     = 300,
                            min_tokens = 20,
                            ...) {

  if (!category_col %in% names(data)) {
    stop("`category_col` must be a column in `data`.")
  }
  if (!all(features %in% names(data))) {
    stop("All `features` must be columns in `data`.")
  }

  # Drop rows with missing category or features
  keep_cols <- c(category_col, features)
  df <- data[stats::complete.cases(data[, keep_cols, drop = FALSE]), keep_cols, drop = FALSE]

  n <- nrow(df)
  if (n < min_tokens) {
    stop("Not enough tokens after removing missing values. Got ", n, ", need at least ", min_tokens, ".")
  }

  cats <- unique(df[[category_col]])
  if (length(cats) != 2L) {
    stop("`category_col` must have exactly two categories in `data`.")
  }

  # Point estimate on full data
  jsd_point <- jsd_kde_nd(
    data     = df,
    features = features,
    group    = category_col,
    ...
  )

  # Bootstrap
  jsd_vals <- numeric(n_boot)
  n_kept   <- 0L

  for (b in seq_len(n_boot)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    boot_df <- df[idx, , drop = FALSE]

    # Guard: some bootstrap samples may lose one category if data are tiny
    if (length(unique(boot_df[[category_col]])) < 2L) {
      jsd_vals[b] <- NA_real_
    } else {
      jsd_vals[b] <- jsd_kde_nd(
        data     = boot_df,
        features = features,
        group    = category_col,
        ...
      )
      n_kept <- n_kept + 1L
    }
  }

  jsd_vals <- jsd_vals[!is.na(jsd_vals)]

  if (!length(jsd_vals)) {
    stop("All bootstrap samples lost one of the categories. Try increasing n or reducing n_boot.")
  }

  data.frame(
    n_tokens = n,
    n_boot   = length(jsd_vals),
    jsd_point = jsd_point,
    jsd_mean  = mean(jsd_vals),
    jsd_sd    = stats::sd(jsd_vals),
    jsd_low   = stats::quantile(jsd_vals, 0.025),
    jsd_high  = stats::quantile(jsd_vals, 0.975),
    stringsAsFactors = FALSE
  )
}
