expect_relative_equal <- function(actual, expected, tolerance = 1e-9) {
  relative_error <- abs(actual - expected) /
    max(abs(expected), .Machine$double.eps)
  expect_true(
    relative_error <= tolerance,
    info = paste0(
      "relative error ", format(relative_error, digits = 16),
      " exceeds ", format(tolerance, scientific = TRUE)
    )
  )
}

test_that("default Pillai return and base-R value remain unchanged", {
  set.seed(101)
  data <- data.frame(
    category = rep(c("a", "b"), each = 30),
    f1 = c(rnorm(30), rnorm(30, 0.8)),
    f2 = c(rnorm(30), rnorm(30, 0.3))
  )

  default <- pillai_overlap(data, c("f1", "f2"), "category")
  explicit_default <- pillai_overlap(
    data,
    c("f1", "f2"),
    "category",
    proportion_standardized = FALSE
  )
  base_v <- summary(
    stats::manova(cbind(f1, f2) ~ category, data = data),
    test = "Pillai"
  )$stats[1, "Pillai"]

  expect_named(default, c("pillai", "p_value"))
  expect_identical(default, explicit_default)
  expect_equal(default$pillai, base_v, tolerance = 1e-6)

  one_feature <- pillai_overlap(data, "f1", "category")
  one_fit <- stats::lm(f1 ~ category, data = data)
  one_aov <- stats::anova(one_fit)
  expected_one <- one_aov[1, "Sum Sq"] / sum(one_aov[, "Sum Sq"])
  expect_named(one_feature, c("pillai", "p_value"))
  expect_equal(one_feature$pillai, expected_one, tolerance = 1e-12)

  one_standardized <- pillai_overlap(
    data,
    "f1",
    "category",
    proportion_standardized = TRUE
  )
  one_expected <- getFromNamespace(".pillai_standardized_fields", "phontrast")(
    one_feature$pillai, 30, 30, 1
  )
  expect_identical(one_standardized[c("pillai", "p_value")], one_feature)
  expect_equal(
    one_standardized[setdiff(names(one_standardized), c("pillai", "p_value"))],
    one_expected,
    tolerance = 1e-12
  )
})

test_that("worked estimator-chain vectors are reproduced", {
  chain <- getFromNamespace(".pillai_standardized_fields", "phontrast")
  vectors <- data.frame(
    id = paste0("TV", 1:8),
    n1 = c(50, 50, 20, 10, 5, 20, 2, 30),
    n2 = c(50, 50, 80, 90, 95, 80, 58, 30),
    p = c(2, 2, 2, 2, 2, 3, 2, 2),
    V = c(0.15, 0.02, 0.12, 0.02, 0.15, 0.12, 0.25, 0.035088),
    H = c(50, 50, 32, 18, 9.5, 32, 3.8666666666666667, 30),
    bias = c(
      0.08, 0.08, 0.125, 0.2222222222222222,
      0.42105263157894735, 0.1875, 1.0344827586206897,
      0.13333333333333333
    ),
    d2_plugin = c(
      0.6917647058823531, 0.08, 0.8352272727272727,
      0.2222222222222222, 3.640866873065016, 0.8352272727272727,
      10, 0.14060722635846584
    ),
    d2_unbiased = c(
      0.5905882352941177, -0.0024489795918367363,
      0.6846590909090908, -0.006802721088435382,
      3.108359133126936, 0.6136363636363635,
      8.448275862068964, 0.0000011054548670341724
    ),
    fallback = c(FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE),
    pillai_eq = c(
      0.12865197334700154, NA, 0.14614918132201332, NA,
      0.43728222996515687, 0.13300492610837436,
      0.6786703601108033, 2.763636403816603e-7
    ),
    d2_fallback = c(
      NA, 0.07755102040816327, NA, 0.21541950113378683,
      NA, NA, NA, NA
    ),
    fragile = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE)
  )

  expected_names <- c(
    "n1", "n2", "H", "d2_plugin", "d2_unbiased", "pillai_eq",
    "pillai_eq_fallback", "d2_fallback", "bias_2p_over_H",
    "fragile_minority"
  )

  for (i in seq_len(nrow(vectors))) {
    v <- vectors[i, ]
    out <- chain(v$V, v$n1, v$n2, v$p)

    expect_identical(names(out), expected_names, info = v$id)
    expect_identical(out$n1, as.integer(v$n1), info = v$id)
    expect_identical(out$n2, as.integer(v$n2), info = v$id)
    expect_relative_equal(out$H, v$H)
    expect_relative_equal(out$bias_2p_over_H, v$bias)
    expect_relative_equal(out$d2_plugin, v$d2_plugin)
    expect_relative_equal(out$d2_unbiased, v$d2_unbiased)
    expect_identical(out$pillai_eq_fallback, v$fallback, info = v$id)
    expect_identical(out$fragile_minority, v$fragile, info = v$id)

    if (is.na(v$pillai_eq)) {
      expect_true(is.na(out$pillai_eq), info = v$id)
    } else {
      expect_relative_equal(out$pillai_eq, v$pillai_eq)
    }
    if (is.na(v$d2_fallback)) {
      expect_true(is.na(out$d2_fallback), info = v$id)
    } else {
      expect_relative_equal(out$d2_fallback, v$d2_fallback)
    }
  }

  boundary <- chain(0.035088, 30, 30, 2)
  expect_gt(boundary$d2_unbiased, 0)
  expect_false(boundary$pillai_eq_fallback)
})

test_that("public opt-in appends fields and labels fallback results", {
  points <- data.frame(
    f1 = c(0, 1, 0, 1),
    f2 = c(0, 0, 1, 1)
  )
  data <- rbind(
    data.frame(category = "a", points),
    data.frame(category = "b", points)
  )

  raw <- pillai_overlap(data, c("f1", "f2"), "category")
  standardized <- pillai_overlap(
    data,
    c("f1", "f2"),
    "category",
    proportion_standardized = TRUE
  )

  expect_identical(standardized[c("pillai", "p_value")], raw)
  expect_named(
    standardized,
    c(
      "pillai", "p_value", "n1", "n2", "H", "d2_plugin",
      "d2_unbiased", "pillai_eq", "pillai_eq_fallback",
      "d2_fallback", "bias_2p_over_H", "fragile_minority"
    )
  )
  expect_true(standardized$pillai_eq_fallback)
  expect_true(is.na(standardized$pillai_eq))
  expect_equal(standardized$d2_fallback, 0, tolerance = 1e-15)
  expect_lt(standardized$d2_unbiased, 0)
  expect_equal(
    standardized$d2_unbiased,
    standardized$d2_fallback - standardized$bias_2p_over_H,
    tolerance = 1e-15
  )
})

test_that("category counts use model factor-level order after filtering", {
  set.seed(102)
  data <- data.frame(
    category = factor(
      c(rep("a", 20), rep("b", 80), "a"),
      levels = c("b", "a", "unused")
    ),
    f1 = rnorm(101),
    f2 = rnorm(101)
  )
  data$f1[101] <- NA_real_

  out <- pillai_overlap(
    data,
    c("f1", "f2"),
    "category",
    proportion_standardized = TRUE
  )

  expect_identical(out$n1, 80L)
  expect_identical(out$n2, 20L)
  expect_equal(out$H, 32)
  expect_false(out$fragile_minority)
})

test_that("definedness and design failures are explicit", {
  undefined <- data.frame(
    category = c("a", "a", "b", "b", "b"),
    f1 = c(0, 1, 0, 1, 2),
    f2 = c(0, 0, 1, 2, 1)
  )
  expect_warning(
    out <- pillai_overlap(
      undefined,
      c("f1", "f2"),
      "category",
      proportion_standardized = TRUE
    ),
    "undefined"
  )
  optional <- setdiff(names(out), c("pillai", "p_value"))
  expect_true(all(vapply(out[optional], function(x) length(x) == 1L && is.na(x),
                         logical(1))))
  expect_type(out$n1, "integer")
  expect_type(out$H, "double")
  expect_type(out$pillai_eq_fallback, "logical")

  too_small <- data.frame(
    category = c("a", rep("b", 5)),
    f1 = c(0, 1:5)
  )
  expect_error(
    pillai_overlap(
      too_small,
      "f1",
      "category",
      proportion_standardized = TRUE
    ),
    "at least two"
  )

  singular <- data.frame(
    category = rep(c("a", "b"), each = 5),
    f1 = c(0:4, 1:5),
    f2 = 2 * c(0:4, 1:5)
  )
  expect_error(
    pillai_overlap(
      singular,
      c("f1", "f2"),
      "category",
      proportion_standardized = TRUE
    ),
    "nonsingular"
  )

  # Near-collinear features must fail the same way. The guard decides rank by
  # tolerance-based QR (with orders of magnitude of margin here), not by
  # whether chol() throws -- the latter differs between BLAS builds on
  # (near-)singular input, which is what CRAN's tests-MKL check tripped on.
  near_singular <- singular
  near_singular$f2 <- near_singular$f2 +
    1e-9 * c(1, -1, 2, -2, 3, -3, 4, -4, 5, -5)
  expect_error(
    pillai_overlap(
      near_singular,
      c("f1", "f2"),
      "category",
      proportion_standardized = TRUE
    ),
    "nonsingular"
  )

  three_classes <- data.frame(
    category = rep(c("a", "b", "c"), each = 4),
    f1 = seq_len(12)
  )
  expect_error(
    pillai_overlap(
      three_classes,
      "f1",
      "category",
      proportion_standardized = TRUE
    ),
    "exactly two"
  )
  expect_error(
    pillai_overlap(
      undefined,
      c("f1", "f2"),
      "category",
      proportion_standardized = NA
    ),
    "TRUE or FALSE"
  )
})

test_that("internal chain guards invalid scores and sample sizes", {
  chain <- getFromNamespace(".pillai_standardized_fields", "phontrast")

  expect_error(chain(-0.01, 20, 20, 2), "\\[0, 1\\)")
  expect_error(chain(1, 20, 20, 2), "\\[0, 1\\)")
  expect_error(chain(Inf, 20, 20, 2), "\\[0, 1\\)")
  expect_error(chain(NA_real_, 20, 20, 2), "\\[0, 1\\)")
  expect_error(chain(0.2, 1, 20, 2), "at least 2")
  expect_error(chain(0.2, 20.5, 20, 2), "positive integer")
})

test_that("end-to-end values agree with base-R MANOVA and the chain", {
  set.seed(2026)
  n1 <- 20L
  n2 <- 80L
  A <- cbind(rnorm(n1, 0, 1), rnorm(n1, 0, 1))
  B <- cbind(rnorm(n2, 0.9, 1), rnorm(n2, 0.4, 1))
  data <- data.frame(
    class = rep(c("a", "b"), c(n1, n2)),
    F1 = c(A[, 1], B[, 1]),
    F2 = c(A[, 2], B[, 2])
  )

  base_v <- summary(
    stats::manova(cbind(F1, F2) ~ class, data = data),
    test = "Pillai"
  )$stats[1, "Pillai"]
  out <- pillai_overlap(
    data,
    c("F1", "F2"),
    "class",
    proportion_standardized = TRUE
  )
  expected <- getFromNamespace(".pillai_standardized_fields", "phontrast")(
    base_v, n1, n2, 2
  )

  expect_equal(out$pillai, base_v, tolerance = 1e-6)
  expect_equal(base_v, 0.205456307, tolerance = 1e-6)
  expect_relative_equal(out$d2_plugin, expected$d2_plugin, tolerance = 1e-12)
  expect_relative_equal(out$d2_unbiased, expected$d2_unbiased, tolerance = 1e-12)
  expect_identical(out$pillai_eq_fallback, FALSE)
  expect_relative_equal(out$pillai_eq, expected$pillai_eq, tolerance = 1e-12)
  expect_equal(out$d2_plugin, 1.58382716, tolerance = 1e-8)
  expect_equal(out$d2_unbiased, 1.41034265, tolerance = 1e-8)
  expect_equal(out$pillai_eq, 0.260675293, tolerance = 1e-8)
})

test_that("closed-form distribution and expectation anchors hold", {
  expect_equal(
    stats::qbeta(0.95, 1, 28.5),
    0.09977758022528305,
    tolerance = 1e-14
  )
  expect_equal(
    stats::qbeta(0.95, 1, 48.5),
    0.05989873024908154,
    tolerance = 1e-14
  )
  expect_equal(
    1.44 / (4 + 1.44),
    0.2647058823529411,
    tolerance = 1e-14
  )

  fallback_probability <- function(N, p) {
    nu_e <- N - 2
    stats::pf(
      (nu_e - p + 1) / (nu_e - p - 1),
      p,
      nu_e - p + 1
    )
  }
  expect_equal(fallback_probability(40, 2), 0.642293423, tolerance = 1e-9)
  expect_equal(fallback_probability(100, 2), 0.635946052, tolerance = 1e-9)
  expect_equal(fallback_probability(500, 2), 0.632862003, tolerance = 1e-9)
  expect_equal(1 - 1 / exp(1), 0.632120559, tolerance = 1e-9)
  expect_equal(fallback_probability(100, 3), 0.613214289, tolerance = 1e-9)

  expected_plugin_at_merger <- function(H) {
    N <- 100
    p <- 2
    nu_e <- N - 2
    2 * nu_e * p / (H * (nu_e - p - 1))
  }
  expect_equal(
    expected_plugin_at_merger(50),
    0.08252631578947368,
    tolerance = 1e-14
  )
  expect_equal(
    expected_plugin_at_merger(9.5),
    0.4343490304709141,
    tolerance = 1e-14
  )
})

test_that("merger fallback rate follows the split-independent law", {
  chain <- getFromNamespace(".pillai_standardized_fields", "phontrast")
  replicates <- 4000L
  target <- stats::pf(97 / 95, 2, 97)
  fallback_tolerance <- 3 * sqrt(target * (1 - target) / replicates)
  size_tolerance <- 3 * sqrt(0.05 * 0.95 / replicates)
  critical_v <- stats::qbeta(0.95, 1, 48.5)

  simulate_split <- function(n1, n2) {
    N <- n1 + n2
    nu_e <- N - 2
    fallback <- logical(replicates)
    d2_unbiased <- numeric(replicates)
    reject <- logical(replicates)

    for (i in seq_len(replicates)) {
      A <- matrix(rnorm(n1 * 2), ncol = 2)
      B <- matrix(rnorm(n2 * 2), ncol = 2)
      mean_A <- colMeans(A)
      mean_B <- colMeans(B)
      centered_A <- sweep(A, 2, mean_A)
      centered_B <- sweep(B, 2, mean_B)
      error_sscp <- crossprod(centered_A) + crossprod(centered_B)
      delta <- mean_A - mean_B
      t2 <- n1 * n2 / N * drop(
        crossprod(delta, solve(error_sscp / nu_e, delta))
      )
      V <- t2 / (t2 + nu_e)
      estimate <- chain(V, n1, n2, 2)

      fallback[i] <- estimate$pillai_eq_fallback
      d2_unbiased[i] <- estimate$d2_unbiased
      reject[i] <- V > critical_v
    }

    list(
      fallback_rate = mean(fallback),
      d2_mean = mean(d2_unbiased),
      d2_se = stats::sd(d2_unbiased) / sqrt(replicates),
      size = mean(reject)
    )
  }

  # Fixed R acceptance seed; do not change without recording new diagnostics.
  set.seed(2026)
  results <- list(
    balanced = simulate_split(50L, 50L),
    imbalanced = simulate_split(5L, 95L)
  )

  for (name in names(results)) {
    result <- results[[name]]
    expect_lt(
      abs(result$fallback_rate - target),
      fallback_tolerance
    )
    expect_gt(result$fallback_rate, 0.5)
    expect_lt(abs(result$d2_mean), 3 * result$d2_se)
    expect_lt(abs(result$size - 0.05), size_tolerance)
  }
})
