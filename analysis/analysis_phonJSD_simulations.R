## ============================================================
## 0. Setup
## ============================================================

set.seed(5241988) # My birthday

# Core packages
library(phonJSD)   # This package
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)
library(mvtnorm)   # multivariate normal simulation
library(MASS)      # for LDA (MVN section only)

## ============================================================
## 0a. Helpers
## ============================================================

# Ensure feature columns are named x1, x2, ..., xd
rename_features_x <- function(df, category_col = "category") {
  feat_idx <- which(names(df) != category_col)
  names(df)[feat_idx] <- paste0("x", seq_along(feat_idx))
  df
}

# Guardrail: ensure exactly two categories
check_two_categories <- function(df, category_col = "category") {
  if (!category_col %in% names(df)) {
    stop("check_two_categories(): `", category_col, "` is not a column in `df`.")
  }
  k <- dplyr::n_distinct(df[[category_col]])
  if (k != 2L) {
    stop(
      "check_two_categories(): `", category_col, "` must have exactly 2 levels; found ",
      k, " (", paste(unique(df[[category_col]]), collapse = ", "), ")."
    )
  }
  invisible(df)
}

## ============================================================
## 0b. Metrics for MVN section (full battery)
## ============================================================

compute_all_metrics_mvn <- function(df,
                                    features,
                                    category_col = "category") {
  # Two categories only
  check_two_categories(df, category_col)
  
  # ---- JSD via KDE (divergence in bits) ----
  jsd_div <- jsd_kde_nd(
    data     = df,
    features = features,
    group    = category_col
  )
  jsd_val <- jsd_div
  
  # ---- Pillai (MANOVA) ----
  pill_row <- tryCatch(
    estimate_pillai(
      data         = df,
      features     = features,
      category_col = category_col,
      group_col    = NULL
    ),
    error = function(e) {
      tibble::tibble(
        pillai  = NA_real_,
        p_value = NA_real_
      )
    }
  )
  
  # ---- Bhattacharyya (MVN) ----
  bh_row <- tryCatch(
    estimate_bhatt(
      data         = df,
      features     = features,
      category_col = category_col,
      group_col    = NULL
    ),
    error = function(e) {
      tibble::tibble(
        bhatt_dist     = NA_real_,
        bhatt_affinity = NA_real_
      )
    }
  )
  
  tibble::tibble(
    jsd        = jsd_val,
    pillai     = pill_row$pillai[1],
    pillai_p   = pill_row$p_value[1],
    bhatt_dist = bh_row$bhatt_dist[1],
    bhatt_aff  = bh_row$bhatt_affinity[1]
  )
}

lda_accuracy_mvn <- function(df,
                             features,
                             category_col = "category",
                             test_prop    = 0.3) {
  df <- df %>% mutate(id_row = dplyr::row_number())
  n  <- nrow(df)
  n_test <- max(2, floor(test_prop * n))
  
  test_ids  <- sample(df$id_row, n_test)
  train_ids <- setdiff(df$id_row, test_ids)
  
  train <- df %>% filter(id_row %in% train_ids)
  test  <- df %>% filter(id_row %in% test_ids)
  
  form <- as.formula(
    paste(category_col, "~", paste(features, collapse = " + "))
  )
  
  lda_fit <- tryCatch(
    MASS::lda(form, data = train),
    error = function(e) NULL
  )
  if (is.null(lda_fit)) return(NA_real_)
  
  pred <- tryCatch(
    predict(lda_fit, newdata = test),
    error = function(e) NULL
  )
  if (is.null(pred)) return(NA_real_)
  
  mean(pred$class == test[[category_col]])
}

## ============================================================
## 0c. MVN simulator used in Section 1
## ============================================================

simulate_gaussian_pair <- function(n_per_cat = 100,
                                   d         = 3,
                                   mean_sep  = 1,
                                   cov_scale = 1) {
  # Means for the two categories
  mu1 <- rep(0, d)
  mu2 <- c(mean_sep, rep(0, d - 1))
  
  # Shared covariance matrix
  Sigma <- diag(cov_scale, d)
  
  # Sample from MVN
  X1 <- mvtnorm::rmvnorm(n_per_cat, mean = mu1, sigma = Sigma)
  X2 <- mvtnorm::rmvnorm(n_per_cat, mean = mu2, sigma = Sigma)
  
  # Bind into a single data frame with a category label
  df <- rbind(
    tibble::tibble(category = "A", as.data.frame(X1)),
    tibble::tibble(category = "B", as.data.frame(X2))
  )
  
  # Make sure feature columns are named x1, x2, ..., xd
  df <- rename_features_x(df, category_col = "category")
  df
}

## ============================================================
## 1. Convergent validity under MVN assumptions (sequential)
## ============================================================

# Grid of sample sizes and separations
n_grid    <- c(20, 50, 100, 200)
sep_grid  <- c(0.5, 1, 1.5, 2)  # in SD units along x1
d_dim     <- 3                  # 3-dimensional space

# *** If runtime is painful, drop this to 50 or 30; or reduce n_grid/sep_grid ***
n_rep     <- 50                 # repetitions per condition

sim_mvnorm_results <- tidyr::crossing(
  n_per_cat = n_grid,
  mean_sep  = sep_grid,
  rep       = 1:n_rep
) %>%
  dplyr::mutate(
    data = purrr::map2(
      n_per_cat,
      mean_sep,
      ~ simulate_gaussian_pair(
        n_per_cat = .x,
        d         = d_dim,
        mean_sep  = .y,
        cov_scale = 1
      )
    )
  ) %>%
  dplyr::mutate(
    metrics = purrr::map(
      data,
      ~ compute_all_metrics_mvn(
        df           = .x,
        features     = paste0("x", 1:d_dim),
        category_col = "category"
      )
    ),
    acc_lda = purrr::map_dbl(
      data,
      ~ lda_accuracy_mvn(
        df           = .x,
        features     = paste0("x", 1:d_dim),
        category_col = "category"
      )
    )
  ) %>%
  dplyr::select(-data) %>%
  tidyr::unnest(cols = metrics)

# Summaries: correlation and scaling across metrics
mvnorm_summary <- sim_mvnorm_results %>%
  dplyr::group_by(n_per_cat, mean_sep) %>%
  dplyr::summarize(
    jsd_mean        = mean(jsd),
    pillai_mean     = mean(pillai),
    bhatt_mean      = mean(bhatt_dist),
    lda_acc_mean    = mean(acc_lda, na.rm = TRUE),
    cor_jsd_pillai  = cor(jsd, pillai,     use = "pairwise.complete.obs"),
    cor_jsd_bhatt   = cor(jsd, bhatt_dist, use = "pairwise.complete.obs"),
    cor_jsd_acc     = cor(jsd, acc_lda,    use = "pairwise.complete.obs"),
    .groups         = "drop"
  )

mvnorm_summary

# Example plot: JSD vs Pillai as a function of separation
ggplot(mvnorm_summary,
       aes(x = jsd_mean, y = pillai_mean,
           color = factor(mean_sep), size = n_per_cat)) +
  geom_point() +
  labs(color = "Mean separation\n(SD units)",
       size  = "n per\ncategory",
       x     = "Mean JSD (divergence, bits)",
       y     = "Mean Pillai trace",
       title = "JSD vs. Pillai as a function of separation\n and n per category") +
  theme_minimal()+theme_tq()+
  scale_color_tq()+theme(legend.position='bottom')+
  guides(
    size  = guide_legend(order = 1, ncol = 1),
    color = guide_legend(order = 2, ncol = 1)
  )+coord_flip()

## ============================================================
## 2. Stress Test Simulations
## ============================================================

compute_jsd <- function(df) {
  jsd_kde_nd(
    data     = df,
    features = c("x1", "x2"),
    group    = "category"
  )
}

# Mixture vs Gaussian (non-normality, non-linear-ish)
simulate_mixture <- function(n = 200, d = 2, sep = 1.5) {
  muA    <- rep(0, d)
  SigmaA <- diag(1, d)
  
  muB1   <- c(sep,  rep(0, d - 1))
  muB2   <- c(-sep, rep(0, d - 1))
  SigmaB <- diag(1, d)
  
  ## Category A: unimodal
  X_A <- mvtnorm::rmvnorm(n, mean = muA, sigma = SigmaA)
  
  ## Category B: mixture, generated with a simple for-loop (no map_dfr)
  comp <- sample(c(1L, 2L), n, replace = TRUE)
  X_B  <- matrix(NA_real_, nrow = n, ncol = d)
  for (i in seq_len(n)) {
    mu_i <- if (comp[i] == 1L) muB1 else muB2
    X_B[i, ] <- mvtnorm::rmvnorm(1, mean = mu_i, sigma = SigmaB)
  }
  
  colnames(X_A) <- paste0("x", 1:d)
  colnames(X_B) <- paste0("x", 1:d)
  
  df <- rbind(
    data.frame(category = "A", X_A, check.names = FALSE),
    data.frame(category = "B", X_B, check.names = FALSE)
  )
  
  tibble::as_tibble(df)
}

# Heteroskedastic Gaussians (same means, different covariances)
simulate_heterosk <- function(n = 200, d = 2, scale = 3) {
  mu1 <- rep(0, d)
  mu2 <- rep(0, d)
  
  Sigma1 <- diag(1,    d)
  Sigma2 <- diag(scale, d)
  
  X1 <- mvtnorm::rmvnorm(n, mean = mu1, sigma = Sigma1)
  X2 <- mvtnorm::rmvnorm(n, mean = mu2, sigma = Sigma2)
  
  colnames(X1) <- paste0("x", 1:d)
  colnames(X2) <- paste0("x", 1:d)
  
  df <- rbind(
    data.frame(category = "A", X1, check.names = FALSE),
    data.frame(category = "B", X2, check.names = FALSE)
  )
  
  tibble::as_tibble(df)
}

# Nonlinear boundary: two half-moons
simulate_moons <- function(n = 200, sep = 1.5) {
  theta <- runif(n, 0, pi)
  
  A <- cbind(
    x1 = cos(theta),
    x2 = sin(theta)
  )
  B <- cbind(
    x1 = cos(theta) + sep,
    x2 = -sin(theta)
  )
  
  df <- rbind(
    data.frame(category = "A", A, check.names = FALSE),
    data.frame(category = "B", B, check.names = FALSE)
  )
  
  tibble::as_tibble(df)
}

stress_grid <- tibble::tibble(
  scenario = c("Mixture vs Gaussian", "Heteroskedastic Gaussians", "Nonlinear Boundary (Moons)")
)

stress_results <- stress_grid %>%
  dplyr::mutate(
    samples = purrr::map(
      scenario,
      function(scen) {
        purrr::map_dfr(seq_len(200), function(i) {
          df <- switch(
            scen,
            "Mixture vs Gaussian"  = simulate_mixture(),
            "Heteroskedastic Gaussians" = simulate_heterosk(),
            "Nonlinear Boundary (Moons)"    = simulate_moons(),
            stop("Unknown scenario: ", scen)
          )
          
          tibble::tibble(
            scenario = scen,
            rep      = i,
            jsd      = compute_jsd(df)
          )
        })
      }
    )
  ) %>%
  dplyr::select(-scenario) %>%   # <— drop outer scenario
  tidyr::unnest(samples)


## ============================================================
## Pillai under the same stress scenarios
## ============================================================

# Helper to get global Pillai on x1,x2 ~ category
compute_pillai_stress <- function(df) {
  res <- tryCatch(
    estimate_pillai(
      data         = df,
      features     = c("x1", "x2"),
      category_col = "category",
      group_col    = NULL
    ),
    error = function(e) {
      # If MANOVA fails (e.g., singular covariance), just return NA
      tibble::tibble(pillai = NA_real_, p_value = NA_real_)
    }
  )
  res$pillai[1]
}

# Run the same 200 reps per scenario, now collecting Pillai
stress_results_pillai <- stress_grid %>%
  dplyr::mutate(
    samples = purrr::map(
      scenario,
      function(scen) {
        purrr::map_dfr(seq_len(200), function(i) {
          df <- switch(
            scen,
            "Mixture vs Gaussian"  = simulate_mixture(),
            "Heteroskedastic Gaussians" = simulate_heterosk(),
            "Nonlinear Boundary (Moons)"    = simulate_moons(),
            stop("Unknown scenario: ", scen)
          )
          tibble::tibble(
            scenario = scen,
            rep      = i,
            pillai   = compute_pillai_stress(df)
          )
        })
      }
    )
  ) %>%
  dplyr::select(-scenario) %>%   # avoid duplicate column on unnest
  tidyr::unnest(samples)

# Quick check / summary
dplyr::group_by(stress_results_pillai, scenario) %>%
  dplyr::summarise(
    pillai_mean = mean(pillai, na.rm = TRUE),
    pillai_sd   = sd(pillai,   na.rm = TRUE),
    .groups = "drop"
  )

# Comparison Graph
## ============================================================
## Combine JSD + Pillai for a single faceted comparison
## ============================================================

# 1. Convert JSD results to long format
df_jsd <- stress_results %>%
  dplyr::select(scenario, rep, jsd) %>%
  dplyr::mutate(
    metric = "JSD",
    value  = jsd
  ) %>%
  dplyr::select(scenario, rep, metric, value)

# 2. Convert Pillai results to long format
df_pillai <- stress_results_pillai %>%
  dplyr::select(scenario, rep, pillai) %>%
  dplyr::mutate(
    metric = "Pillai",
    value  = pillai
  ) %>%
  dplyr::select(scenario, rep, metric, value)

# 3. Bind the two long dataframes
df_both <- dplyr::bind_rows(df_jsd, df_pillai)

# 4. Combined figure
ggplot(df_both,
       aes(x = scenario, y = value, fill = metric)) +
  geom_boxplot(position = position_dodge(width = 0.75),
               outlier.alpha = 0.4) +
  labs(
    x = "Scenario",
    y = "Metric value",
    fill = "Metric"  ) +
  theme_tq()+
  scale_fill_tq()+
  theme_minimal(base_size = 14)+theme(legend.position='bottom')+
  coord_flip()



ggplot(mvnorm_long,
       aes(x = value, y = acc_lda, color = metric)) +
  geom_point(alpha = 0.25, size = 1) +
  geom_smooth(method = "gam", se = FALSE, size = 1) +
  facet_grid(
    n_per_cat ~ mean_sep,
    labeller = labeller(
      n_per_cat = function(x) paste0("n = ", x, " per cat."),
      mean_sep  = function(x) paste0("Δμ = ", x, " SD")
    )
  ) +
  labs(
    x     = "Metric value",
    y     = "LDA classification accuracy",
    color = "Metric",
    title = "Perceptual proxy: LDA accuracy vs JSD and Pillai under multivariate normal assumptions"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.spacing = unit(1.2, "lines"),
    legend.position='bottom'
  )+scale_color_tq()

mvnorm_acc_corr <- sim_mvnorm_results %>%
  group_by(n_per_cat, mean_sep) %>%
  summarise(
    cor_jsd_acc    = cor(jsd,    acc_lda, use = "pairwise.complete.obs"),
    cor_pillai_acc = cor(pillai, acc_lda, use = "pairwise.complete.obs"),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols      = c(cor_jsd_acc, cor_pillai_acc),
    names_to  = "metric",
    values_to = "correlation"
  )

mvnorm_acc_corr$metric <- factor(
  mvnorm_acc_corr$metric,
  levels = c("cor_jsd_acc", "cor_pillai_acc"),
  labels = c("JSD", "Pillai")
)

ggplot(mvnorm_acc_corr,
       aes(x = factor(mean_sep),
           y = correlation,
           fill = metric)) +
  geom_col(position = position_dodge(width = 0.7)) +
  facet_wrap(~ n_per_cat, labeller = label_both) +
  labs(
    x     = "Mean separation (Δμ in SD units)",
    y     = "Correlation with LDA accuracy",
    fill  = "Metric",
    title = "Convergent 'perceptual' validity:\ncor(JSD/Pillai, LDA accuracy) by separation and sample size"
  ) +
  ylim(0, 1) +
  theme_minimal(base_size = 11)+theme_tq()+
  scale_fill_tq()

## ============================================================
## 3. Empirical sanity check on vowel data
## ============================================================

## Option A: PB52 /E/ vs /I/ (Petersen & Barney) via phonTools
## --------------------------------------------------------

if (requireNamespace("phonTools", quietly = TRUE)) {
  data(pb52, package = "phonTools")
}

if (exists("pb52")) {
  eh_ih <- pb52 %>%
    dplyr::filter(vowel %in% c("E", "I")) %>%
    droplevels()
  
  features_formant <- c("f1", "f2")
  
  # Global estimates
  jsd_pb <- estimate_jsd(
    data         = eh_ih,
    features     = features_formant,
    category_col = "vowel",
    do_boot      = TRUE,
    n_boot       = 300
  )
  
  pillai_pb <- estimate_pillai(
    data         = eh_ih,
    features     = features_formant,
    category_col = "vowel"
  )
  
  bhatt_pb <- estimate_bhatt(
    data         = eh_ih,
    features     = features_formant,
    category_col = "vowel"
  )
  
  overlap_pb <- estimate_overlap(
    data         = eh_ih,
    features     = features_formant,
    category_col = "vowel"
  )
  
  jsd_pb
  pillai_pb
  bhatt_pb
  overlap_pb
}

library(ggplot2)
library(dplyr)

if (exists("eh_ih")) {
  
  # Make F1/F2 Bark-like orientation (bigger F1 = lower)
  eh_ih %>% mutate(vowel=case_when(vowel=='E'~'/ɛ/',
                                   vowel=='I'~'/ɪ/')) %>%
  ggplot(aes(x = f2, y = f1, color = vowel)) +
    stat_ellipse(level = 0.95, linewidth = 1) +
    geom_point(alpha = 0.25) +
    scale_x_reverse() + scale_y_reverse() +
    labs(
      title = "Petersen & Barney (1952): /ɛ/ vs /ɪ/",
      subtitle = paste0(
        "JSD = ", round(jsd_pb$jsd_point, 3),
        "   |   Pillai = ", round(pillai_pb$pillai, 3),
        "   |   Bhatt = ", round(bhatt_pb$bhatt_dist, 3),
        "   |   Overlap = ", round(overlap_pb$overlap, 3)
      ),
      x = "F2",
      y = "F1",
      color = "Vowel"
    ) +
    scale_color_brewer(palette='Set1')+
    theme_minimal(base_size = 14)+theme_tq()
}

if (exists("eh_ih")) {
  eh_ih %>% mutate(vowel=case_when(vowel=='E'~'/ɛ/',
                                   vowel=='I'~'/ɪ/')) %>%
  ggplot(., aes(x = f2, y = f1, color = vowel, fill = vowel)) +
    stat_density_2d(geom = "polygon", alpha = 0.25, contour_var = "ndensity") +
    scale_x_reverse() + scale_y_reverse() +
    labs(
      title = "Density Overlap for /ɛ/ vs /ɪ/ (PB52)",
      subtitle = paste0(
        "JSD = ", round(jsd_pb$jsd_point, 3),
        " | Pillai = ", round(pillai_pb$pillai, 3),
        " | Bhatt = ", round(bhatt_pb$bhatt_dist, 3),
        " | Overlap = ", round(overlap_pb$overlap, 3)
      ),
      x = "F2",
      y = "F1"
    ) +
    theme_minimal(base_size = 14) + theme_tq()+
    scale_color_tq()+scale_fill_tq()+
    guides(fill = "none")
}

library(tibble)
library(ggplot2)

if (exists("eh_ih")) {
  
  metrics_df <- tibble(
    Metric = c("JSD", "Pillai", "Bhattacharyya", "Overlap"),
    Value  = c(
      jsd_pb$jsd_point,
      pillai_pb$pillai,
      bhatt_pb$bhatt_dist,
      overlap_pb$overlap
    )
  )
  
  ggplot(metrics_df, aes(x = Metric, y = Value)) +
    geom_col(fill = "#4E79A7") +
    geom_text(aes(label = round(Value, 3)), vjust = -0.4, size = 5) +
    labs(
      title = "Contrast Metrics for /ɛ/ vs /ɪ/ (PB52)",
      y = "Value",
      x = ""
    ) +
    theme_minimal(base_size = 14)
}

## Option B: your own vowel data with EH/IH (commented example)
## ------------------------------------------------------------

# vowels_eh_ih <- vowels_df %>%
#   dplyr::filter(vowel %in% c("EH", "IH")) %>%
#   dplyr::mutate(
#     category = if_else(vowel == "EH", "eh", "ih")
#   )
#
# features_formant <- c("F1", "F2")
#
# jsd_vowels <- estimate_jsd(
#   data         = vowels_eh_ih,
#   features     = features_formant,
#   category_col = "category",
#   group_col    = "speaker",
#   do_boot      = TRUE,
#   n_boot       = 300,
#   min_tokens   = 6
# )
#
# pillai_vowels <- estimate_pillai(
#   data         = vowels_eh_ih,
#   features     = features_formant,
#   category_col = "category",
#   group_col    = "speaker",
#   min_tokens   = 6
# )
#
# vowel_compare <- jsd_vowels %>%
#   dplyr::left_join(pillai_vowels, by = c("group", "n_tokens"))
#
# ggplot(vowel_compare,
#        aes(x = jsd_mean, y = pillai)) +
#   geom_point() +
#   geom_smooth(method = "lm", se = FALSE) +
#   labs(x = "JSD (EH~IH, mean bootstrap)",
#        y = "Pillai (F1, F2 MANOVA)") +
#   theme_minimal()