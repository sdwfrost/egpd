#' egpd: Extended Generalized Pareto Distribution GAMs
#'
#' The egpd package fits Extended Generalized Pareto Distribution (EGPD),
#' Discrete EGPD (DEGPD), Zero-Inflated EGPD (ZIEGPD), and Zero-Inflated
#' Discrete EGPD (ZIDEGPD) models within a GAM framework, with additional
#' support for fitting via \code{gamlss} and \code{bamlss}. The fitting
#' infrastructure is adapted from the \code{evgam} package by Ben Youngman.
#' P-value testing code is adapted from the \code{mgcv} package by Simon Wood.
#'
#' @section Main function:
#' The main function is \code{\link{egpd}}.
#'
#' @references
#' Abbas, A., Ahmad, T. and Ahmad, I. (2025).
#' Modeling zero-inflated precipitation extremes.
#' \emph{arXiv preprint} arXiv:2504.11058.
#' \url{https://arxiv.org/abs/2504.11058}
#'
#' Ahmad, T. and Arshad, I. A. (2024).
#' New flexible versions of extended generalized Pareto model for count data.
#' \emph{arXiv preprint} arXiv:2409.18719.
#' \url{https://arxiv.org/abs/2409.18719}
#'
#' Ahmad, T. and Hussain, A. (2025).
#' Flexible model for varying levels of zeros and outliers in count data.
#' \emph{arXiv preprint} arXiv:2510.27365.
#' \url{https://arxiv.org/abs/2510.27365}
#'
#' @docType package
#' @name egpd-package
#'
#' @useDynLib egpd, .registration = TRUE
#' @importFrom Rcpp sourceCpp evalCpp
#' @import stats graphics grDevices utils
NULL
