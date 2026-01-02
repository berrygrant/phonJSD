#' Kullback–Leibler divergence for discrete distributions
#'
#' Computes KL(p || q) in bits for discrete probability vectors.
#' Zeros are handled by replacing them with a small epsilon.
#'
#' @param p,q Numeric probability vectors of the same length.
#'
#' @return A single numeric value: the KL divergence in bits.
#' @export
kl_div <- function(p, q) {
  stopifnot(length(p) == length(q))
  eps <- .Machine$double.eps
  p_safe <- ifelse(p == 0, eps, p)
  q_safe <- ifelse(q == 0, eps, q)
  sum(p_safe * (log(p_safe / q_safe) / log(2)))
}

#' Jensen–Shannon divergence for discrete distributions
#'
#' Computes JSD(p, q) in bits for discrete probability vectors.
#' The value is bounded in \code{[0, 1]} for equally weighted mixtures.
#'
#' @inheritParams kl_div
#'
#' @return A single numeric value: the Jensen–Shannon divergence in bits.
#' @export
jsd <- function(p, q) {
  stopifnot(length(p) == length(q))
  p <- p / sum(p)
  q <- q / sum(q)
  m <- 0.5 * (p + q)
  0.5 * kl_div(p, m) + 0.5 * kl_div(q, m)
}
