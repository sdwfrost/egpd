// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "gfunc_derivs.h"

// Shared stabilized Jet + log-sum-exp evaluator for ZIDEGPD models 1..6.
arma::mat zidegpd_loglink_eval(const Rcpp::List& pars, const std::vector<arma::mat>& X,
                               const arma::vec& y, const arma::uvec& dupid, int dcate,
                               const Rcpp::List& offsets, int m, bool value_only = false);

const double xieps = 0.0;

// [[Rcpp::export]]
double zidegpd1d0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2, const arma::mat& X3, const arma::mat& X4, arma::vec yvec, const arma::uvec& dupid, int dcate, const Rcpp::List& offsets)
{
  return arma::accu(zidegpd_loglink_eval(pars, {X1, X2, X3, X4}, yvec, dupid, dcate, offsets, 1, true).col(0));
}

// [[Rcpp::export]]
arma::mat zidegpd1d12(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat X4, arma::vec yvec, const arma::uvec& dupid, int dcate, const Rcpp::List& offsets)
{
  return zidegpd_loglink_eval(pars, {X1, X2, X3, X4}, yvec, dupid, dcate, offsets, 1).cols(1, 14);
}

// [[Rcpp::export]]
double zidegpd2d0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2, const arma::mat& X3, const arma::mat& X4, const arma::mat& X5, const arma::mat& X6, arma::vec yvec, const arma::uvec& dupid, int dcate, const Rcpp::List& offsets)
{
  return arma::accu(zidegpd_loglink_eval(pars, {X1, X2, X3, X4, X5, X6}, yvec, dupid, dcate, offsets, 2, true).col(0));
}

// [[Rcpp::export]]
arma::mat zidegpd2d12(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat X4, arma::mat X5, arma::mat X6, arma::vec yvec, const arma::uvec dupid, int dcate, const Rcpp::List& offsets)
{
  return zidegpd_loglink_eval(pars, {X1, X2, X3, X4, X5, X6}, yvec, dupid, dcate, offsets, 2).cols(1, 27);
}

// [[Rcpp::export]]
double zidegpd3d0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2, const arma::mat& X3, const arma::mat& X4, arma::vec yvec, const arma::uvec& dupid, int dcate, const Rcpp::List& offsets)
{
  return arma::accu(zidegpd_loglink_eval(pars, {X1, X2, X3, X4}, yvec, dupid, dcate, offsets, 3, true).col(0));
}

// [[Rcpp::export]]
arma::mat zidegpd3d12(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat X4, arma::vec yvec, const arma::uvec dupid, int dcate, const Rcpp::List& offsets)
{
  return zidegpd_loglink_eval(pars, {X1, X2, X3, X4}, yvec, dupid, dcate, offsets, 3).cols(1, 14);
}

// [[Rcpp::export]]
double zidegpd4d0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2, const arma::mat& X3, const arma::mat& X4, const arma::mat& X5, arma::vec yvec, const arma::uvec& dupid, int dcate, const Rcpp::List& offsets)
{
  return arma::accu(zidegpd_loglink_eval(pars, {X1, X2, X3, X4, X5}, yvec, dupid, dcate, offsets, 4, true).col(0));
}

// [[Rcpp::export]]
arma::mat zidegpd4d12(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat X4, arma::mat X5, arma::vec yvec, const arma::uvec dupid, int dcate, const Rcpp::List& offsets)
{
  return zidegpd_loglink_eval(pars, {X1, X2, X3, X4, X5}, yvec, dupid, dcate, offsets, 4).cols(1, 20);
}

// [[Rcpp::export]]
double zidegpd5d0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2, const arma::mat& X3, const arma::mat& X4, arma::vec yvec, const arma::uvec& dupid, int dcate, const Rcpp::List& offsets)
{
  return arma::accu(zidegpd_loglink_eval(pars, {X1, X2, X3, X4}, yvec, dupid, dcate, offsets, 5, true).col(0));
}

// [[Rcpp::export]]
arma::mat zidegpd5d12(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat X4, arma::vec yvec, const arma::uvec& dupid, int dcate, const Rcpp::List& offsets)
{
  return zidegpd_loglink_eval(pars, {X1, X2, X3, X4}, yvec, dupid, dcate, offsets, 5).cols(1, 14);
}

// [[Rcpp::export]]
double zidegpd6d0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2, const arma::mat& X3, const arma::mat& X4, arma::vec yvec, const arma::uvec& dupid, int dcate, const Rcpp::List& offsets)
{
  return arma::accu(zidegpd_loglink_eval(pars, {X1, X2, X3, X4}, yvec, dupid, dcate, offsets, 6, true).col(0));
}

// [[Rcpp::export]]
arma::mat zidegpd6d12(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat X4, arma::vec yvec, const arma::uvec& dupid, int dcate, const Rcpp::List& offsets)
{
  return zidegpd_loglink_eval(pars, {X1, X2, X3, X4}, yvec, dupid, dcate, offsets, 6).cols(1, 14);
}
