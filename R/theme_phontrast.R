# Shared visual identity for phontrast plots: an Okabe-Ito palette and a
# publication-oriented theme. All plotting functions in the package apply these
# by default; they are exported so users can restyle their own ggplots to match.

# Okabe-Ito palette, reordered so the first two colors give the
# highest-contrast, colorblind-safe pair for a two-category contrast.
.phontrast_colors <- c(
  blue       = "#0072B2",
  vermillion = "#D55E00",
  green      = "#009E73",
  orange     = "#E69F00",
  purple     = "#CC79A7",
  skyblue    = "#56B4E9",
  yellow     = "#F0E442",
  grey       = "#999999"
)

# Ink and grid colors used by theme_phontrast() and annotation layers.
.phontrast_ink <- "#1A1A1A"
.phontrast_ink_muted <- "#4D4D4D"
.phontrast_grid <- "#E6E6E6"

#' Okabe-Ito color palette used by phontrast plots
#'
#' Returns \code{n} colors from the Okabe-Ito palette, a qualitative palette
#' designed to be distinguishable under the common forms of color-vision
#' deficiency. The palette is reordered so the first two colors (blue and
#' vermillion) form the highest-contrast pair for two-category contrasts.
#' For more than eight groups the palette is interpolated, with a warning,
#' since interpolated qualitative colors lose their guarantees.
#'
#' @param n Number of colors to return. Defaults to the full palette.
#'
#' @return A character vector of \code{n} hex colors, named for \code{n <= 8}.
#' @examples
#' phontrast_palette()
#' phontrast_palette(2)
#' @export
phontrast_palette <- function(n = NULL) {
  if (is.null(n)) {
    return(.phontrast_colors)
  }
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1 ||
      n != as.integer(n)) {
    stop("`n` must be a single positive integer.", call. = FALSE)
  }
  n <- as.integer(n)
  if (n <= length(.phontrast_colors)) {
    return(.phontrast_colors[seq_len(n)])
  }
  warning(
    "phontrast_palette() has ", length(.phontrast_colors),
    " qualitative colors; interpolating to ", n,
    " loses the colorblind-safety guarantees.",
    call. = FALSE
  )
  grDevices::colorRampPalette(unname(.phontrast_colors))(n)
}

#' Okabe-Ito discrete color and fill scales
#'
#' Discrete \pkg{ggplot2} scales built on \code{phontrast_palette()}. These are
#' the default scales for all phontrast plotting functions and can be added to
#' any ggplot to match the package's visual identity.
#'
#' @param ... Passed to \code{ggplot2::discrete_scale()} (e.g. \code{name},
#'   \code{labels}, \code{guide}).
#'
#' @return A \pkg{ggplot2} scale object.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   ggplot2::ggplot(iris, ggplot2::aes(Sepal.Length, Sepal.Width,
#'                                      color = Species)) +
#'     ggplot2::geom_point() +
#'     scale_color_phontrast() +
#'     theme_phontrast()
#' }
#' @export
scale_colour_phontrast <- function(...) {
  .phontrast_discrete_scale("colour", ...)
}

#' @rdname scale_colour_phontrast
#' @export
scale_color_phontrast <- function(...) {
  .phontrast_discrete_scale("colour", ...)
}

#' @rdname scale_colour_phontrast
#' @export
scale_fill_phontrast <- function(...) {
  .phontrast_discrete_scale("fill", ...)
}

.phontrast_discrete_scale <- function(aesthetics, ...) {
  .require_ggplot2()
  args <- list(
    aesthetics = aesthetics,
    palette = function(n) unname(phontrast_palette(n)),
    ...
  )
  # ggplot2 < 3.5.0 requires `scale_name` (no default); 3.5.0+ deprecates it.
  # Supply it only where it is still mandatory. The empty-symbol check is done
  # inline: binding the missing-argument marker to a local and reading it
  # would itself raise a missing-argument error.
  if (identical(formals(ggplot2::discrete_scale)[["scale_name"]], quote(expr = ))) {
    args$scale_name <- "phontrast"
  }
  do.call(ggplot2::discrete_scale, args)
}

#' Publication theme for phontrast plots
#'
#' A minimal, publication-oriented \pkg{ggplot2} theme: quiet major grid, no
#' minor grid, bold plot-aligned title, muted subtitle, left-aligned caption
#' for estimator provenance, and a top-aligned legend. Applied by default in
#' all phontrast plotting functions.
#'
#' @param base_size Base font size in points.
#' @param base_family Base font family.
#'
#' @return A \pkg{ggplot2} theme object.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   ggplot2::ggplot(iris, ggplot2::aes(Sepal.Length, Sepal.Width)) +
#'     ggplot2::geom_point() +
#'     theme_phontrast()
#' }
#' @export
theme_phontrast <- function(base_size = 12, base_family = "") {
  .require_ggplot2()
  .check_plot_number(base_size, "base_size", lower = 1)
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        color = .phontrast_grid, linewidth = 0.4
      ),
      axis.title = ggplot2::element_text(color = .phontrast_ink),
      axis.text = ggplot2::element_text(color = .phontrast_ink_muted),
      plot.title = ggplot2::element_text(
        color = .phontrast_ink, face = "bold", size = ggplot2::rel(1.1)
      ),
      plot.subtitle = ggplot2::element_text(
        color = .phontrast_ink_muted, size = ggplot2::rel(0.9)
      ),
      plot.caption = ggplot2::element_text(
        color = .phontrast_ink_muted, size = ggplot2::rel(0.75), hjust = 0
      ),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      strip.text = ggplot2::element_text(color = .phontrast_ink, face = "bold"),
      legend.position = "top",
      legend.justification = "left",
      legend.title = ggplot2::element_text(color = .phontrast_ink),
      legend.text = ggplot2::element_text(color = .phontrast_ink_muted),
      plot.margin = ggplot2::margin(12, 16, 10, 12)
    )
}

# Apply the house style (theme + discrete scales) to a package-built plot.
# Fill scale is attached only when a layer maps the fill aesthetic, so plots
# without fills do not carry a dormant scale.
.phontrast_style <- function(p, fill = FALSE, base_size = 12) {
  p <- p + theme_phontrast(base_size = base_size) + scale_colour_phontrast()
  if (isTRUE(fill)) {
    p <- p + scale_fill_phontrast()
  }
  p
}
