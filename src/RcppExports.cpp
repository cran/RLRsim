// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// RLRsimCpp
List RLRsimCpp(int p, int k, int n, int nsim, int g, int q, Rcpp::NumericVector mu, Rcpp::NumericVector lambda, double lambda0, Rcpp::NumericVector xi, bool REML);
RcppExport SEXP RLRsim_RLRsimCpp(SEXP pSEXP, SEXP kSEXP, SEXP nSEXP, SEXP nsimSEXP, SEXP gSEXP, SEXP qSEXP, SEXP muSEXP, SEXP lambdaSEXP, SEXP lambda0SEXP, SEXP xiSEXP, SEXP REMLSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< int >::type g(gSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< bool >::type REML(REMLSEXP);
    rcpp_result_gen = Rcpp::wrap(RLRsimCpp(p, k, n, nsim, g, q, mu, lambda, lambda0, xi, REML));
    return rcpp_result_gen;
END_RCPP
}
