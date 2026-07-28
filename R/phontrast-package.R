#' phontrast: Contrast and Separation Metrics for Phonological Categories
#'
#' A unified toolkit for quantifying the separation and overlap between
#' phonological categories (e.g., vowels, consonants) in arbitrary
#' n-dimensional acoustic spaces such as formants, MFCCs, spectral features,
#' duration, or learned embeddings. The main entry point, \code{phontrast()},
#' computes and compares multiple contrast metrics in one call: Jensen-Shannon
#' divergence and distance, the Pillai-Bartlett trace, Bhattacharyya distance
#' and affinity, Mahalanobis distance, and proportional overlap, globally or by
#' group, with optional bootstrap intervals. Percent overlap is returned as a
#' 0--1 proportion, not a 0--100 percentage.
#'
#' phontrast was formerly released as phonJSD (through version 1.2.0); see the
#' package NEWS for migration notes.
#'
#' @section Recommended workflow:
#' \enumerate{
#'   \item Start with \code{phontrast()} to compute one or more contrast metrics
#'     for a two-category contrast, globally or by group, in tidy long or wide
#'     form.
#'   \item Use \code{estimate_jsd()} when Jensen-Shannon divergence or
#'     Jensen-Shannon distance is the primary outcome and you need optional
#'     bootstrap intervals.
#'   \item Use lower-level helpers such as \code{jsd_kde_nd()},
#'     \code{percent_overlap_kde()}, \code{pillai_overlap()}, and
#'     \code{bhattacharyya_mvnorm()} when validating methods, debugging one
#'     contrast, or reproducing a specific metric.
#'   \item Use \code{plot_contrast()} for a distribution-aware, annotated view
#'     of one contrast that draws the same density model the metrics use, and
#'     \code{plot_overlap_metrics()} (also available as \code{plot()} /
#'     \code{ggplot2::autoplot()} on \code{phontrast()} results),
#'     \code{plot_category_space()}, and \code{plot_category_pca()} for
#'     ggplot2-backed diagnostics and presentation figures, all sharing the
#'     colorblind-safe \code{theme_phontrast()} visual identity.
#' }
#'
#' @section Choosing metrics:
#' JSD, Jensen-Shannon distance, Pillai trace, Bhattacharyya distance, and
#' Mahalanobis distance increase as categories become more separated. Percent
#' overlap and Bhattacharyya affinity increase as categories overlap more. The
#' long output from \code{phontrast()} includes an \code{orientation} column and
#' a separation-oriented \code{separation_value} column to make these directions
#' explicit. JSD and percent overlap estimate distributional separation/overlap
#' using KDE by default; Pillai and Mahalanobis emphasize mean separation;
#' Bhattacharyya metrics use a multivariate-normal approximation.
#'
#' @section Density backends:
#' The distributional metrics (Jensen-Shannon divergence and proportional
#' overlap) are computed from a density estimate for each category. The
#' \code{density} argument selects that estimate: \code{"kde"} (the default)
#' uses kernel density estimation, and \code{"mvnorm"} fits one multivariate
#' normal per category and estimates the metric between the two Gaussians by
#' Monte-Carlo (with \code{mc_n} samples and reproducible \code{eval_seed}).
#' The \code{"mvnorm"} backend matches the estimator behind JSD and overlap to
#' the same multivariate-normal assumptions the Pillai, Bhattacharyya, and
#' Mahalanobis columns already make, and is convenient for higher-dimensional
#' feature spaces where multivariate KDE is impractical. It is available on
#' \code{phontrast()}, \code{estimate_jsd()}, \code{estimate_overlap()},
#' \code{jsd_summary()}, \code{global_boot_jsd()}, \code{jsd_kde_nd()}, and
#' \code{percent_overlap_kde()}; the parametric metrics are unaffected by it.
#'
#' @section High-dimensional workflows:
#' Metrics can be estimated in arbitrary n-dimensional numeric feature spaces,
#' including MFCCs and learned embeddings. Use \code{plot_category_pca()} for a
#' two-dimensional PCA diagnostic, but report metric estimates from the intended
#' full feature set.
#'
#' @section Confidence intervals:
#' Confidence intervals use \code{ci_lower} and \code{ci_upper} columns.
#' Legacy JSD aliases \code{jsd_low} and \code{jsd_high} are retained for
#' compatibility.
#'
#' @keywords internal
"_PACKAGE"
