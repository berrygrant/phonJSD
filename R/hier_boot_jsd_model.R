#' Hierarchical bootstrap for JSD-based models
#'
#' Performs a hierarchical bootstrap: resample groups with replacement,
#' resample tokens within each sampled group, compute JSD per group,
#' fit a model to the bootstrap JSD values, and repeat.
#'
#' This lets you propagate measurement uncertainty in JSD into model
#' parameters (e.g., GAM/LMM coefficients).
#'
#' @param data Data frame with at least: group_col, category_col, features,
#'   and any predictors used in the model.
#' @param group_col String: grouping variable (e.g., "speaker").
#' @param category_col String: category variable with 2 levels (e.g., "vowel").
#' @param features Character vector of acoustic feature columns.
#' @param formula Model formula to pass to `fit_fun` (e.g.,
#'   `jsd_beta ~ s(age) + s(region, bs = "re")`).
#' @param fit_fun A function that takes `(formula, data, ...)` and returns a
#'   fitted model. Defaults to `mgcv::gam` if available, otherwise `stats::lm`.
#' @param n_outer Number of hierarchical bootstrap replicates.
#' @param min_tokens Minimum within-group tokens required.
#' @param eps Small epsilon for bounding JSD in (0, 1) if using Beta family.
#' @param progress Logical; if TRUE, prints progress every 10 replicates.
#' @param ... Additional arguments passed to `fit_fun`.
#'
#' @return A tibble with columns:
#'   \itemize{
#'     \item \code{boot_id} – bootstrap replicate index
#'     \item \code{term} – model term
#'     \item \code{estimate} – estimate for that term in that replicate
#'   }
#' @export
#' @importFrom dplyr group_by ungroup summarize n first bind_rows left_join across all_of
#' @importFrom purrr map
#' @importFrom tibble tibble
#' @importFrom stats coef
#' @importFrom rlang .data
hier_boot_jsd_model <- function(data,
                                group_col,
                                category_col,
                                features,
                                formula,
                                fit_fun   = NULL,
                                n_outer   = 200,
                                min_tokens = 30,
                                eps       = 1e-6,
                                progress  = TRUE,
                                ...) {

  if (is.null(fit_fun)) {
    if (requireNamespace("mgcv", quietly = TRUE)) {
      fit_fun <- mgcv::gam
    } else {
      fit_fun <- stats::lm
    }
  }

  groups <- unique(data[[group_col]])
  n_groups <- length(groups)

  boot_results <- vector("list", n_outer)

  for (b in seq_len(n_outer)) {
    if (progress && b %% 10 == 0) {
      message("Bootstrap replicate ", b, " / ", n_outer)
    }

    # Resample groups with replacement
    boot_group_ids <- sample(groups, size = n_groups, replace = TRUE)

    boot_df_list <- lapply(boot_group_ids, function(g) {
      df_g <- data[data[[group_col]] == g, , drop = FALSE]
      if (nrow(df_g) < min_tokens ||
          dplyr::n_distinct(df_g[[category_col]]) < 2L) {
        return(NULL)
      }
      # Resample tokens within group
      df_g[sample.int(nrow(df_g), size = nrow(df_g), replace = TRUE), ]
    })

    boot_df <- dplyr::bind_rows(boot_df_list)
    if (nrow(boot_df) == 0L) {
      next
    }

    # Compute JSD per group for this bootstrap sample
    jsd_sum <- jsd_summary(
      data         = boot_df,
      group_col    = group_col,
      category_col = category_col,
      features     = features,
      do_boot      = FALSE,
      min_tokens   = min_tokens
    )

    # Merge in predictors for modeling (e.g., age, region)
    preds <- data |>
      dplyr::group_by(.data[[group_col]]) |>
      dplyr::summarize(
        dplyr::across(
          .cols = !dplyr::all_of(c(category_col, features)),
          .fns  = dplyr::first
        ),
        .groups = "drop"
      )

    model_df <- jsd_sum |>
      dplyr::left_join(preds, by = c("group" = group_col))

    # Prepare jsd_beta using point JSD
    model_df <- prepare_jsd_beta(
      jsd_df = model_df,
      jsd_col = "jsd_point",
      eps = eps
    )

    # Fit model
    fit <- try(fit_fun(formula = formula, data = model_df, ...),
               silent = TRUE)
    if (inherits(fit, "try-error")) {
      next
    }

    cf <- stats::coef(fit)
    boot_results[[b]] <- tibble::tibble(
      boot_id  = b,
      term     = names(cf),
      estimate = unname(cf)
    )
  }

  dplyr::bind_rows(boot_results)
}
