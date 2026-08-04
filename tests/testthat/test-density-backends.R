# Regression tests for the pluggable multivariate-normal density backend
# (density = "mvnorm") added in phontrast 2.2.0. The backend fits one Gaussian
# per category (base-R Cholesky log-density and sampler, no extra dependency)
# and estimates JSD / overlap between the fitted Gaussians by fresh-sample
# Monte-Carlo.

test_that(".mvn_logdens matches the closed-form Gaussian log-density", {
  # Diagonal covariance reduces to a sum of univariate normal log-densities,
  # which we can check against stats::dnorm without any extra package.
  set.seed(1)
  X <- matrix(rnorm(30), ncol = 3)
  mu <- c(0.2, -0.4, 1.1)
  s2 <- c(0.5, 2.0, 1.3)
  ref <- rowSums(vapply(
    seq_along(mu),
    function(j) stats::dnorm(X[, j], mu[j], sqrt(s2[j]), log = TRUE),
    numeric(nrow(X))
  ))
  got <- phontrast:::.mvn_logdens(X, mu, diag(s2))
  expect_equal(got, ref, tolerance = 1e-10)
})

test_that(".mvn_rsample reproduces the target mean and covariance", {
  set.seed(2)
  mu <- c(1, -2, 0.5)
  S <- crossprod(matrix(rnorm(9), 3)) + diag(3)   # a valid positive-definite S
  Y <- phontrast:::.mvn_rsample(2e5, mu, S)
  expect_equal(ncol(Y), 3L)
  expect_equal(colMeans(Y), mu, tolerance = 0.03)
  expect_equal(stats::cov(Y), S, tolerance = 0.05)
})

test_that("mvnorm backend recovers the JSD/overlap between the fitted Gaussians", {
  skip_if_not_installed("mvtnorm")
  # The mvnorm backend targets the JSD / overlap between the two *fitted*
  # Gaussians (one per category), a well-defined estimand that differs from the
  # population truth only by the sampling variability in the fit. We therefore
  # fit the same Gaussians the backend fits (ridge = 1e-6, its default),
  # grid-integrate the exact JSD / overlap between them, and confirm the
  # fresh-sample Monte-Carlo estimate reproduces that reference to within
  # Monte-Carlo error.
  set.seed(101)
  n <- 600
  A <- mvtnorm::rmvnorm(n, c(0, 0), diag(2))
  B <- mvtnorm::rmvnorm(n, c(1.2, 0.8), diag(2))
  d <- data.frame(cat = rep(c("A", "B"), each = n),
                  x = c(A[, 1], B[, 1]), y = c(A[, 2], B[, 2]))

  X1 <- as.matrix(d[d$cat == "A", c("x", "y")])
  X2 <- as.matrix(d[d$cat == "B", c("x", "y")])
  f_mu1 <- colMeans(X1); f_S1 <- stats::cov(X1) + diag(1e-6, 2)
  f_mu2 <- colMeans(X2); f_S2 <- stats::cov(X2) + diag(1e-6, 2)

  gx <- seq(-6, 8, length.out = 320); da <- (gx[2] - gx[1])^2
  G <- as.matrix(expand.grid(gx, gx))
  p <- mvtnorm::dmvnorm(G, f_mu1, f_S1); q <- mvtnorm::dmvnorm(G, f_mu2, f_S2)
  m <- 0.5 * (p + q)
  jsd_fit <- 0.5 * sum(ifelse(p > 0, p * log2(p / m), 0)) * da +
             0.5 * sum(ifelse(q > 0, q * log2(q / m), 0)) * da
  ovl_fit <- sum(pmin(p, q)) * da

  jsd_mvn <- jsd_kde_nd(d, c("x", "y"), "cat", density = "mvnorm",
                        mc_n = 40000, eval_seed = 7)
  ovl_mvn <- percent_overlap_kde(d, c("x", "y"), "cat", density = "mvnorm",
                                 mc_n = 40000, eval_seed = 7)
  expect_equal(jsd_mvn, jsd_fit, tolerance = 0.01)
  expect_equal(ovl_mvn, ovl_fit, tolerance = 0.01)
})

test_that("mvnorm backend runs in high dimensions and stays in range", {
  # No KDE bandwidth selection is involved, so the parametric backend handles
  # dimensionalities where multivariate KDE is impractical.
  set.seed(3)
  n <- 120; k <- 8
  feats <- paste0("m", seq_len(k))
  shift <- c(1.1, 0.9, rep(0, k - 2))
  mat <- matrix(rnorm(2 * n * k), ncol = k)
  mat[(n + 1):(2 * n), ] <- sweep(mat[(n + 1):(2 * n), ], 2, shift, "+")
  d <- data.frame(cat = rep(c("a", "b"), each = n)); d[feats] <- mat

  jsd_hd <- jsd_kde_nd(d, feats, "cat", density = "mvnorm", eval_seed = 5)
  ovl_hd <- percent_overlap_kde(d, feats, "cat", density = "mvnorm", eval_seed = 5)
  expect_true(is.finite(jsd_hd) && jsd_hd >= 0 && jsd_hd <= 1)
  expect_true(is.finite(ovl_hd) && ovl_hd >= 0 && ovl_hd <= 1)
})

test_that("mvnorm backend regularizes near-singular covariance and validates ridge", {
  # Near-collinear features give a strongly rank-deficient covariance. The
  # default ridge keeps the fitted covariance positive definite, so the backend
  # still returns a finite, in-range estimate rather than failing.
  set.seed(4)
  f1 <- rnorm(60)
  d <- data.frame(cat = rep(c("a", "b"), each = 30),
                  f1 = f1, f2 = 2 * f1 + rnorm(60, 0, 0.05))
  jsd_nc <- jsd_kde_nd(d, c("f1", "f2"), "cat", density = "mvnorm", eval_seed = 1)
  expect_true(is.finite(jsd_nc) && jsd_nc >= 0 && jsd_nc <= 1)

  # A non-positive ridge is rejected up front.
  expect_error(
    phontrast:::.mvnorm_mc_pair(d, c("f1", "f2"), "cat", ridge = 0),
    "positive"
  )
})

test_that("mvnorm JSD is invariant to category label order and ignores loo", {
  set.seed(9)
  d <- data.frame(cat = rep(c("a", "b"), each = 70),
                  f1 = c(rnorm(70, 0), rnorm(70, 1.3)),
                  f2 = c(rnorm(70, 0), rnorm(70, 0.9)))
  d2 <- d; d2$cat <- factor(d2$cat, levels = c("b", "a"))
  expect_equal(jsd_kde_nd(d, c("f1", "f2"), "cat", density = "mvnorm", eval_seed = 1),
               jsd_kde_nd(d2, c("f1", "f2"), "cat", density = "mvnorm", eval_seed = 1))
  # The Gaussian fit has no self-kernel, so `loo` has no effect under mvnorm.
  expect_equal(jsd_kde_nd(d, c("f1", "f2"), "cat", density = "mvnorm", loo = TRUE,  eval_seed = 1),
               jsd_kde_nd(d, c("f1", "f2"), "cat", density = "mvnorm", loo = FALSE, eval_seed = 1))
})

test_that("eval_seed makes the mvnorm estimate reproducible without touching the RNG stream", {
  set.seed(9)
  d <- data.frame(cat = rep(c("a", "b"), each = 60),
                  f1 = c(rnorm(60, 0), rnorm(60, 1)),
                  f2 = c(rnorm(60, 0), rnorm(60, 0.7)))
  a <- jsd_kde_nd(d, c("f1", "f2"), "cat", density = "mvnorm", eval_seed = 42)
  b <- jsd_kde_nd(d, c("f1", "f2"), "cat", density = "mvnorm", eval_seed = 42)
  expect_identical(a, b)
  cc <- jsd_kde_nd(d, c("f1", "f2"), "cat", density = "mvnorm", eval_seed = 7)
  expect_false(isTRUE(all.equal(a, cc)))

  # The private seeded generator must not read or change the caller's stream.
  set.seed(123); r1 <- runif(1)
  set.seed(123)
  invisible(jsd_kde_nd(d, c("f1", "f2"), "cat", density = "mvnorm", eval_seed = 42))
  r2 <- runif(1)
  expect_equal(r1, r2)
})

test_that("seeded evaluation-point sampling leaves the RNG stream unchanged", {
  points <- matrix(seq_len(200), ncol = 2)

  set.seed(456); expected <- runif(1)
  set.seed(456)
  first <- phontrast:::.sample_kde_eval_points(points, 12, eval_seed = 99)
  second <- phontrast:::.sample_kde_eval_points(points, 12, eval_seed = 99)
  observed <- runif(1)

  expect_identical(first, second)
  expect_equal(observed, expected)
})

test_that("density argument threads through the high-level API", {
  set.seed(3)
  d <- data.frame(
    speaker = rep(c("s1", "s2"), each = 80),
    cat = rep(rep(c("a", "b"), each = 40), 2),
    f1 = c(rnorm(40, 0), rnorm(40, 1), rnorm(40, 0.2), rnorm(40, 1.2)),
    f2 = c(rnorm(40, 0), rnorm(40, 1), rnorm(40, 0.1), rnorm(40, 1.1))
  )
  ft <- c("f1", "f2")

  # estimate_jsd / estimate_overlap accept density and differ from the KDE default
  jk <- estimate_jsd(d, ft, "cat")$jsd_point
  jm <- estimate_jsd(d, ft, "cat", density = "mvnorm", eval_seed = 1)$jsd_point
  expect_true(is.finite(jm))
  expect_false(isTRUE(all.equal(jk, jm)))

  ok <- estimate_overlap(d, ft, "cat")$overlap
  om <- estimate_overlap(d, ft, "cat", density = "mvnorm", eval_seed = 1)$overlap
  expect_true(is.finite(om) && om >= 0 && om <= 1)
  expect_false(isTRUE(all.equal(ok, om)))

  # global_boot_jsd and jsd_summary accept density
  gb <- global_boot_jsd(d, ft, "cat", n_boot = 6, density = "mvnorm",
                        mc_n = 2000, eval_seed = 1)
  expect_true(is.finite(gb$jsd_point) && is.finite(gb$jsd_mean))
  js <- jsd_summary(d, "speaker", "cat", ft, do_boot = FALSE,
                    density = "mvnorm", eval_seed = 1)
  expect_true(all(is.finite(js$jsd_point)))
})

test_that("phontrast density = 'mvnorm' changes only the distributional metrics", {
  set.seed(21)
  d <- data.frame(
    cat = rep(c("a", "b"), each = 80),
    f1 = c(rnorm(80, 0), rnorm(80, 1)),
    f2 = c(rnorm(80, 0), rnorm(80, 0.7))
  )
  ft <- c("f1", "f2")
  rk <- phontrast(d, ft, "cat", density = "kde", progress = FALSE)
  rm <- phontrast(d, ft, "cat", density = "mvnorm", eval_seed = 1, progress = FALSE)

  # JSD and overlap respond to the density backend...
  expect_false(isTRUE(all.equal(rk$jsd, rm$jsd)))
  expect_false(isTRUE(all.equal(rk$percent_overlap, rm$percent_overlap)))
  # ...the parametric metrics are unaffected by construction.
  expect_equal(rk$pillai, rm$pillai)
  expect_equal(rk$bhatt_dist, rm$bhatt_dist)
  expect_equal(rk$mahalanobis_dist, rm$mahalanobis_dist)
})

test_that("invalid density is rejected", {
  d <- data.frame(cat = rep(c("a", "b"), each = 30),
                  f1 = c(rnorm(30), rnorm(30, 1)),
                  f2 = c(rnorm(30), rnorm(30, 1)))
  expect_error(jsd_kde_nd(d, c("f1", "f2"), "cat", density = "banana"))
  expect_error(estimate_jsd(d, c("f1", "f2"), "cat", density = "banana"))
})
