#' Prepare JSD estimates for beta regression / GAMs
#'
#' JSD lives in \eqn{[0,1]}. This helper adds a `jsd_beta` column bounded in
#' \eqn{(0,1)} so it can be used with the Beta family (e.g., in \pkg{mgcv}).
#'
#' @param jsd_df Data frame containing a JSD column.
#' @param jsd_col String: name of the JSD column (default "jsd_mean").
#' @param eps Small constant used to bound JSD away from 0 and 1.
#'
#' @return A modified data frame with an added `jsd_beta` column.
#' @export
prepare_jsd_beta <- function(jsd_df,
                             jsd_col = "jsd_mean",
                             eps     = 1e-6) {

  if (!jsd_col %in% names(jsd_df)) {
    stop("`jsd_col` must be a column in `jsd_df`.")
  }
  jsd <- jsd_df[[jsd_col]]
  jsd_beta <- pmin(pmax(jsd, eps), 1 - eps)
  jsd_df$jsd_beta <- jsd_beta
  jsd_df
}
