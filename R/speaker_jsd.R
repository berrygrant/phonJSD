#' Group-level JSD point estimates
#'
#' Computes JSD for each group (e.g., speaker) comparing two categories
#' (e.g., vowels) in an n-dimensional acoustic space.
#'
#' @param data Data frame containing acoustic measurements.
#' @param group_col String: name of column giving the grouping unit
#'   (e.g., "speaker").
#' @param category_col String: name of column giving the category
#'   to compare (e.g., "vowel"). Each group must have exactly two categories.
#' @param features Character vector of column names giving the acoustic space.
#' @param min_tokens Minimum number of tokens per group required to compute
#'   JSD. Groups with fewer tokens are dropped.
#' @param ... Additional arguments passed to \code{jsd_kde_nd()}.
#'
#' @return A tibble with one row per group and columns:
#'   \code{group}, \code{n_tokens}, and \code{jsd}.
#' @export
#' @importFrom dplyr group_by summarize n n_distinct filter ungroup rename
#' @importFrom tibble tibble
#' @importFrom rlang .data
speaker_jsd <- function(data,
                        group_col,
                        category_col,
                        features,
                        min_tokens = 20,
                        ...) {

  if (!group_col %in% names(data)) {
    stop("`group_col` must be a column in `data`.")
  }
  if (!category_col %in% names(data)) {
    stop("`category_col` must be a column in `data`.")
  }

  dplyr::group_by(data, .data[[group_col]]) |>
    dplyr::filter(
      dplyr::n_distinct(.data[[category_col]]) == 2L,
      dplyr::n() >= min_tokens
    ) |>
    dplyr::summarize(
      n_tokens = dplyr::n(),
      jsd      = jsd_kde_nd(
        data     = dplyr::cur_data_all(),
        features = features,
        group    = category_col,
        ...
      ),
      .groups = "drop"
    ) |>
    dplyr::rename(group = 1)
}

#' Bootstrap JSD for each group
#'
#' Computes bootstrap mean, SD, and confidence interval for JSD within each
#' group (e.g., speaker), using resampling with replacement.
#'
#' @inheritParams speaker_jsd
#' @param n_boot Number of bootstrap resamples per group.
#'
#' @return A tibble with one row per group and columns:
#'   \code{group}, \code{n_tokens}, \code{jsd_mean}, \code{jsd_sd},
#'   \code{jsd_low}, and \code{jsd_high}.
#' @export
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @importFrom stats quantile sd
boot_jsd <- function(data,
                     group_col,
                     category_col,
                     features,
                     n_boot     = 300,
                     min_tokens = 30,
                     ...) {

  if (!group_col %in% names(data)) {
    stop("`group_col` must be a column in `data`.")
  }
  if (!category_col %in% names(data)) {
    stop("`category_col` must be a column in `data`.")
  }

  groups <- split(data, data[[group_col]])

  res_list <- purrr::map(groups, function(df_g) {

    if (nrow(df_g) < min_tokens ||
        dplyr::n_distinct(df_g[[category_col]]) < 2L) {
      return(tibble::tibble(
        group    = df_g[[group_col]][1],
        n_tokens = nrow(df_g),
        jsd_mean = NA_real_,
        jsd_sd   = NA_real_,
        jsd_low  = NA_real_,
        jsd_high = NA_real_
      ))
    }

    jsd_vals <- replicate(
      n_boot,
      {
        samp <- df_g[sample.int(nrow(df_g), size = nrow(df_g), replace = TRUE), ]
        jsd_kde_nd(
          data     = samp,
          features = features,
          group    = category_col,
          ...
        )
      }
    )

    tibble::tibble(
      group    = df_g[[group_col]][1],
      n_tokens = nrow(df_g),
      jsd_mean = mean(jsd_vals),
      jsd_sd   = stats::sd(jsd_vals),
      jsd_low  = stats::quantile(jsd_vals, 0.025),
      jsd_high = stats::quantile(jsd_vals, 0.975)
    )
  })

  dplyr::bind_rows(res_list)
}
