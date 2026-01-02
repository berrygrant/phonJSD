## ============================================================
## 0. Setup
## ============================================================

set.seed(2026)

# Core packages
library(phonJSD)   # your package
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)

# For multivariate normal simulation & LDA
library(mvtnorm)
library(MASS)

## Helper: generic metric calculator for two-category data
## ------------------------------------------------------------

compute_all_metrics <- function(df,
                                features,
                                category_col = "category") {
  # JSD via KDE (non-parametric)
  jsd_val <- jsd_kde_nd(
    data     = df,
    features = features,
    group    = category_col
  )
  
  # Pillai (MANOVA) on same features
  pill <- pillai_overlap(
    data        = df,
    features    = features,
    category_col = category_col
  )
  
  # Bhattacharyya under MVN
  bh <- bhattacharyya_mvnorm(
    data        = df,
    features    = features,
    category_col = category_col
  )
  
  tibble(
    jsd        = jsd_val,
    pillai     = pill$pillai,
    pillai_p   = pill$p_value,
    bhatt_dist = bh$distance,
    bhatt_aff  = bh$affinity
  )
}

## Helper: classification accuracy (LDA) as a "perceptual proxy"
## ------------------------------------------------------------

lda_accuracy <- function(df,
                         features,
                         category_col = "category",
                         test_prop    = 0.3) {
  df <- df %>% mutate(id_row = row_number())
  n  <- nrow(df)
  n_test <- max(2, floor(test_prop * n))
  
  test_ids  <- sample(df$id_row, n_test)
  train_ids <- setdiff(df$id_row, test_ids)
  
  train <- df %>% filter(id_row %in% train_ids)
  test  <- df %>% filter(id_row %in% test_ids)
  
  # LDA model
  form <- as.formula(
    paste(category_col, "~", paste(features, collapse = " + "))
  )
  lda_fit <- MASS::lda(form, data = train)
  
  pred <- predict(lda_fit, newdata = test)
  mean(pred$class == test[[category_col]])
}

## ============================================================
## 1. Convergent validity under MVN assumptions
## ============================================================

# Simulate two Gaussian categories with controlled separation & equal cov

simulate_gaussian_pair <- function(n_per_cat = 100,
                                   d         = 3,
                                   mean_sep  = 1,
                                   cov_scale = 1) {
  # Category A mean at 0, Category B mean shifted along first dimension
  mu1 <- rep(0, d)
  mu2 <- c(mean_sep, rep(0, d - 1))
  
  # Shared covariance
  Sigma <- diag(cov_scale, d)
  
  X1 <- mvtnorm::rmvnorm(n_per_cat, mean = mu1, sigma = Sigma)
  X2 <- mvtnorm::rmvnorm(n_per_cat, mean = mu2, sigma = Sigma)
  
  df <- rbind(
    tibble(category = "A", as.data.frame(X1)),
    tibble(category = "B", as.data.frame(X2))
  )
  names(df)[names(df) %in% paste0("V", 1:d)] <- paste0("x", 1:d)
  df
}

# Grid of sample sizes and separations
n_grid    <- c(20, 50, 100, 200)
sep_grid  <- c(0.5, 1, 1.5, 2)  # in SD units along x1
d         <- 3                  # 3-dimensional space
n_rep     <- 200                # repetitions per condition

sim_mvnorm_results <- crossing(
  n_per_cat = n_grid,
  mean_sep  = sep_grid,
  rep       = 1:n_rep
) %>%
  mutate(
    data = purrr::pmap(
      list(n_per_cat, mean_sep),
      ~ simulate_gaussian_pair(
        n_per_cat = ..1,
        d         = d,
        mean_sep  = ..2,
        cov_scale = 1
      )
    )
  ) %>%
  mutate(
    metrics = purrr::map(data, ~ compute_all_metrics(
      df        = .x,
      features  = paste0("x", 1:d),
      category_col = "category"
    )),
    acc_lda = purrr::map_dbl(data, ~ lda_accuracy(
      df        = .x,
      features  = paste0("x", 1:d),
      category_col = "category"
    ))
  ) %>%
  select(-data) %>%
  unnest(cols = metrics)

# Example summaries: correlation and scaling across metrics
mvnorm_summary <- sim_mvnorm_results %>%
  group_by(n_per_cat, mean_sep) %>%
  summarize(
    jsd_mean        = mean(jsd),
    pillai_mean     = mean(pillai),
    bhatt_mean      = mean(bhatt_dist),
    lda_acc_mean    = mean(acc_lda),
    cor_jsd_pillai  = cor(jsd, pillai),
    cor_jsd_bhatt   = cor(jsd, bhatt_dist),
    cor_jsd_acc     = cor(jsd, acc_lda),
    .groups         = "drop"
  )

mvnorm_summary

# Example plot: JSD vs Pillai as a function of separation
ggplot(mvnorm_summary,
       aes(x = jsd_mean, y = pillai_mean,
           color = factor(mean_sep), size = n_per_cat)) +
  geom_point() +
  labs(color = "Mean separation",
       size  = "n per category",
       x     = "Mean JSD",
       y     = "Mean Pillai") +
  theme_minimal()

## ============================================================
## 2. Stress tests: non-Gaussian, heteroskedastic, nonlinear
## ============================================================

# 2a. Non-Gaussian: mixture for one category, unimodal for other
simulate_mixture_vs_gaussian <- function(n_per_cat = 100, d = 3) {
  # Category A: unimodal at 0
  muA <- rep(0, d)
  SigmaA <- diag(1, d)
  
  # Category B: mixture of two equally weighted Gaussians
  muB1 <- c(1.5, rep(0, d - 1))
  muB2 <- c(-1.5, rep(0, d - 1))
  SigmaB <- diag(0.5, d)
  
  comp    <- sample(c(1, 2), n_per_cat, replace = TRUE)
  means_B <- rbind(muB1, muB2)
  
  X_A <- mvtnorm::rmvnorm(n_per_cat, mean = muA, sigma = SigmaA)
  X_B <- purrr::map_dfr(comp, function(k) {
    as_tibble(t(mvtnorm::rmvnorm(1, mean = means_B[k, ], sigma = SigmaB)))
  })
  
  df <- bind_rows(
    tibble(category = "A", X_A),
    tibble(category = "B", X_B)
  )
  names(df)[names(df) %in% paste0("V", 1:d)] <- paste0("x", 1:d)
  df
}

# 2b. Heteroskedastic: same means, different covariances
simulate_heteroskedastic <- function(n_per_cat = 100, d = 3) {
  mu1 <- rep(0, d)
  mu2 <- rep(0, d)
  
  Sigma1 <- diag(c(0.5, rep(0.5, d - 1)))
  Sigma2 <- diag(c(2, rep(0.5, d - 1)))
  
  X1 <- mvtnorm::rmvnorm(n_per_cat, mean = mu1, sigma = Sigma1)
  X2 <- mvtnorm::rmvnorm(n_per_cat, mean = mu2, sigma = Sigma2)
  
  df <- rbind(
    tibble(category = "A", as.data.frame(X1)),
    tibble(category = "B", as.data.frame(X2))
  )
  names(df)[names(df) %in% paste0("V", 1:d)] <- paste0("x", 1:d)
  df
}

# 2c. Nonlinear boundary: two “half-moons” in 2D, plus noise dims
simulate_nonlinear <- function(n_per_cat = 100, d = 3) {
  # base 2D half-moons
  theta_A <- runif(n_per_cat, 0, pi)
  theta_B <- runif(n_per_cat, 0, pi)
  
  r <- 1
  sep <- 1.5
  
  xA1 <- r * cos(theta_A)
  xA2 <- r * sin(theta_A)
  
  xB1 <- r * cos(theta_B) + sep
  xB2 <- -r * sin(theta_B)
  
  X_A <- cbind(xA1, xA2)
  X_B <- cbind(xB1, xB2)
  
  # Add (d-2) noise dimensions
  if (d > 2) {
    noise_A <- matrix(rnorm(n_per_cat * (d - 2), 0, 0.5), ncol = d - 2)
    noise_B <- matrix(rnorm(n_per_cat * (d - 2), 0, 0.5), ncol = d - 2)
    X_A <- cbind(X_A, noise_A)
    X_B <- cbind(X_B, noise_B)
  }
  
  df <- rbind(
    tibble(category = "A", as.data.frame(X_A)),
    tibble(category = "B", as.data.frame(X_B))
  )
  names(df)[names(df) %in% paste0("V", 1:d)] <- paste0("x", 1:d)
  df
}

# Run stress test simulations
n_per_cat_stress <- 100
d_stress <- 3
n_rep_stress <- 200

stress_design <- tibble(
  scenario = c("mixture_vs_gaussian", "heteroskedastic", "nonlinear")
)

simulate_stress_one <- function(scenario) {
  sim_fun <- switch(
    scenario,
    mixture_vs_gaussian = simulate_mixture_vs_gaussian,
    heteroskedastic     = simulate_heteroskedastic,
    nonlinear           = simulate_nonlinear
  )
  
  map_dfr(1:n_rep_stress, function(rep) {
    df <- sim_fun(n_per_cat = n_per_cat_stress, d = d_stress)
    metrics <- compute_all_metrics(
      df        = df,
      features  = paste0("x", 1:d_stress),
      category_col = "category"
    )
    acc <- lda_accuracy(
      df        = df,
      features  = paste0("x", 1:d_stress),
      category_col = "category"
    )
    
    metrics %>%
      mutate(
        scenario = scenario,
        rep      = rep,
        acc_lda  = acc
      )
  })
}

stress_results <- stress_design %>%
  mutate(
    res = purrr::map(scenario, simulate_stress_one)
  ) %>%
  select(-scenario) %>%
  unnest(res)

# Summaries: how do metrics correlate with accuracy under each stress scenario?
stress_summary <- stress_results %>%
  group_by(scenario) %>%
  summarize(
    jsd_mean       = mean(jsd),
    pillai_mean    = mean(pillai),
    bhatt_mean     = mean(bhatt_dist),
    lda_acc_mean   = mean(acc_lda),
    cor_jsd_acc    = cor(jsd, acc_lda),
    cor_pillai_acc = cor(pillai, acc_lda),
    cor_bhatt_acc  = cor(bhatt_dist, acc_lda),
    .groups        = "drop"
  )

stress_summary

# Example plot: metric vs accuracy by scenario
ggplot(stress_results,
       aes(x = jsd, y = acc_lda, color = scenario)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE) +
  labs(x = "JSD", y = "LDA classification accuracy") +
  theme_minimal()

## ============================================================
## 3. Empirical sanity check on vowel data (e.g., PB /ɛ/ vs /ɪ/)
## ============================================================

## Assuming you have a vowel dataset `vowels_df` with:
## - speaker: speaker ID
## - vowel: vowel label (includes eh/ih)
## - F1, F2: in Bark or Hz
## If you already have MFCCs, just change `features`.

# Example: filter to /ɛ/ vs /ɪ/, re-label as two categories
vowels_eh_ih <- vowels_df %>%
  filter(vowel %in% c("EH", "IH")) %>%  # adjust labels for your dataset
  mutate(
    category = if_else(vowel == "EH", "eh", "ih")
  )

features_formant <- c("F1", "F2")

# Per-speaker JSD (point + bootstrap)
jsd_vowels <- jsd_summary(
  data         = vowels_eh_ih,
  group_col    = "speaker",
  category_col = "category",
  features     = features_formant,
  do_boot      = TRUE,
  n_boot       = 300,
  min_tokens   = 6
)

jsd_vowels

# Per-speaker Pillai on same features
pillai_vowels <- speaker_pillai(
  data         = vowels_eh_ih,
  group_col    = "speaker",
  category_col = "category",
  features     = features_formant,
  min_tokens   = 6
)

pillai_vowels

# Merge metrics
vowel_compare <- jsd_vowels %>%
  left_join(pillai_vowels, by = c("group", "n_tokens"))

vowel_compare

# Simple scatter: JSD vs Pillai across speakers
ggplot(vowel_compare,
       aes(x = jsd_mean, y = pillai)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "JSD (eh ~ ih, mean bootstrap)",
       y = "Pillai (F1,F2 MANOVA)") +
  theme_minimal()

# Optional: LDA accuracy for each speaker as perceptual proxy
speaker_lda_acc <- vowels_eh_ih %>%
  group_by(speaker) %>%
  filter(n_distinct(category) == 2) %>%
  summarize(
    lda_acc = lda_accuracy(
      df        = cur_data_all(),
      features  = features_formant,
      category_col = "category"
    ),
    .groups = "drop"
  )

vowel_compare <- vowel_compare %>%
  left_join(speaker_lda_acc, by = c("group" = "speaker"))

# Correlations on real data
cor(vowel_compare$jsd_mean, vowel_compare$pillai,   use = "complete.obs")
cor(vowel_compare$jsd_mean, vowel_compare$lda_acc, use = "complete.obs")
cor(vowel_compare$pillai,    vowel_compare$lda_acc, use = "complete.obs")
