# Tests for the plotting family added in phontrast 2.3.0: the Okabe-Ito
# theme/palette module, the distribution-aware plot_contrast(), and the
# plot()/autoplot() methods on phontrast() results.

.plot_test_data <- function(n = 60) {
  set.seed(2026)
  data.frame(
    vowel = rep(c("ih", "eh"), each = n),
    f1 = c(rnorm(n, 500, 55), rnorm(n, 565, 60)),
    f2 = c(rnorm(n, 1980, 150), rnorm(n, 1870, 155))
  )
}

.plot_test_grouped <- function(n = 40) {
  set.seed(7)
  mk <- function(id, sep) data.frame(
    speaker = id,
    vowel = rep(c("ih", "eh"), each = n),
    f1 = c(rnorm(n, 500, 50), rnorm(n, 500 + sep, 55)),
    f2 = c(rnorm(n, 1980, 140), rnorm(n, 1980 - 2 * sep, 150))
  )
  rbind(mk("s01", 90), mk("s02", 30))
}

.layer_geoms <- function(p) {
  vapply(p$layers, function(l) class(l$geom)[1], character(1))
}

test_that("phontrast_palette returns Okabe-Ito colors and validates n", {
  pal <- phontrast_palette()
  expect_length(pal, 8L)
  # The leading pair is the maximum-contrast CVD-safe blue/vermillion pair.
  expect_identical(unname(pal[1:2]), c("#0072B2", "#D55E00"))
  expect_identical(unname(phontrast_palette(2)), c("#0072B2", "#D55E00"))
  expect_error(phontrast_palette(0), "positive integer")
  expect_warning(big <- phontrast_palette(12), "interpolating")
  expect_length(big, 12L)
})

test_that("theme and scales construct ggplot2 objects", {
  skip_if_not_installed("ggplot2")
  expect_s3_class(theme_phontrast(), "theme")
  expect_s3_class(scale_colour_phontrast(), "Scale")
  expect_s3_class(scale_color_phontrast(), "Scale")
  expect_s3_class(scale_fill_phontrast(), "Scale")
  expect_error(theme_phontrast(base_size = -1), "base_size")
})

test_that("plot_contrast() builds an annotated 2-D KDE plot", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ks")
  d <- .plot_test_data()
  p <- plot_contrast(d, c("f2", "f1"), "vowel", grid_n = 61,
                     reverse_x = TRUE, reverse_y = TRUE)
  expect_s3_class(p, "ggplot")
  geoms <- .layer_geoms(p)
  expect_true("GeomTile" %in% geoms)    # overlap shading
  expect_true("GeomPoint" %in% geoms)   # tokens
  expect_true("GeomPath" %in% geoms)    # HDR contours
  expect_true("GeomLabel" %in% geoms)   # metric annotation
  expect_match(p$labels$caption, "density: kde")
  expect_match(p$labels$caption, "highest-density")

  ann <- attr(p, "contrast_metrics")
  expect_s3_class(ann, "tbl_df")
  expect_true(all(c("jsd", "percent_overlap") %in% names(ann)))
  expect_true(is.finite(ann$jsd[1]))
})

test_that("plot_contrast() draws the mvnorm backend's fitted model", {
  skip_if_not_installed("ggplot2")
  d <- .plot_test_data()
  p <- plot_contrast(d, c("f2", "f1"), "vowel", density = "mvnorm",
                     eval_seed = 11, grid_n = 61)
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$caption, "density: mvnorm")
  expect_match(p$labels$caption, "seed = 11")
  expect_match(p$labels$caption, "coverage")

  # The annotated values come from the mvnorm backend, so they differ from the
  # KDE annotation on the same data.
  pk <- plot_contrast(d, c("f2", "f1"), "vowel", grid_n = 61)
  jm <- attr(p, "contrast_metrics")$jsd[1]
  jk <- attr(pk, "contrast_metrics")$jsd[1]
  expect_false(isTRUE(all.equal(jm, jk)))
})

test_that("plot_contrast() 1-D keeps densities and overlap on one scale", {
  skip_if_not_installed("ggplot2")
  d <- .plot_test_data()
  p <- plot_contrast(d, "f1", "vowel")
  geoms <- .layer_geoms(p)
  expect_true("GeomArea" %in% geoms)    # overlap ribbon
  expect_true("GeomLine" %in% geoms)    # density curves
  expect_true("GeomRug" %in% geoms)     # tokens
  # The ribbon is min(p, q) on the density scale: it can never exceed the
  # curves (this failing was a real bug caught in visual QA).
  ov <- p$layers[[1]]$data
  curves <- p$layers[[3]]$data
  expect_lte(max(ov$ov), max(curves$dens) + 1e-12)
  # 1-D draws no region contours, so the caption must not claim any.
  expect_false(grepl("regions:", p$labels$caption))
})

test_that("plot_contrast() facets by group and skips small groups", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ks")
  d <- .plot_test_grouped()
  tiny <- d[d$speaker == "s02", ][1:10, ]
  tiny$speaker <- "s03"
  expect_warning(
    p <- plot_contrast(rbind(d, tiny), c("f1", "f2"), "vowel",
                       group_col = "speaker", grid_n = 41),
    "skipped 1 group"
  )
  expect_s3_class(p$facet, "FacetWrap")
  ann <- attr(p, "contrast_metrics")
  expect_setequal(ann$group, c("s01", "s02"))
})

test_that("plot_contrast() bootstrap annotation carries intervals", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ks")
  d <- .plot_test_data(40)
  p <- plot_contrast(d, c("f1", "f2"), "vowel", n_boot = 5, grid_n = 41)
  ann <- attr(p, "contrast_metrics")
  expect_true(all(c("jsd_ci_lower", "jsd_ci_upper") %in% names(ann)))
  expect_match(p$labels$caption, "bootstrap")
})

test_that("plot_contrast() annotate = FALSE omits labels and metrics", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ks")
  d <- .plot_test_data(40)
  p <- plot_contrast(d, c("f1", "f2"), "vowel", annotate = FALSE, grid_n = 41)
  expect_false("GeomLabel" %in% .layer_geoms(p))
  expect_null(attr(p, "contrast_metrics"))
})

test_that("plot_contrast() validates its inputs", {
  skip_if_not_installed("ggplot2")
  d <- .plot_test_data(40)
  expect_error(
    plot_contrast(d, c("f1", "f2", "f1"), "vowel"),
    "one or two column names"
  )
  expect_error(
    plot_contrast(transform(d, vowel = "ih"), "f1", "vowel"),
    "exactly two"
  )
  expect_error(
    plot_contrast(d, "f1", "vowel", levels = c(0.5, 1.2)),
    "levels"
  )
  expect_error(
    plot_contrast(d, "f1", "vowel", n_boot = -1),
    "n_boot"
  )
  expect_error(
    plot_contrast(d[c(1:5, 41:45), ], "f1", "vowel"),  # both categories, n = 10
    "min_tokens"
  )
})

test_that("phontrast() results plot directly via plot() and autoplot()", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ks")
  d <- .plot_test_data(40)
  m <- phontrast(d, c("f1", "f2"), "vowel", output = "long", progress = FALSE)
  expect_s3_class(m, "phontrast_contrast")

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  p <- plot(m)
  expect_s3_class(p, "ggplot")
  expect_s3_class(ggplot2::autoplot(m), "ggplot")

  # Wide output carries the class too.
  mw <- phontrast(d, c("f1", "f2"), "vowel", output = "wide", progress = FALSE)
  expect_s3_class(mw, "phontrast_contrast")
})

test_that("retrofitted plots use the phontrast theme", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ks")
  d <- .plot_test_data(40)
  m <- phontrast(d, c("f1", "f2"), "vowel", output = "long", progress = FALSE)
  for (p in list(
    plot_overlap_metrics(m),
    plot_category_space(d, c("f1", "f2"), "vowel"),
    plot_category_space(d, "f1", "vowel"),
    plot_category_pca(d, c("f1", "f2"), "vowel")
  )) {
    expect_s3_class(p, "ggplot")
    # theme_phontrast() positions the title at the plot (not panel) edge.
    expect_identical(p$theme$plot.title.position, "plot")
  }
})
