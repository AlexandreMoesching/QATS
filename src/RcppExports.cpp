// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// QATS_timer_cpp
List QATS_timer_cpp(int d0, int n_seeds, bool rotate, double n_rep, int n, int m, const arma::vec& logPi, const arma::mat& qq, const arma::mat& GG, arma::ivec SS, int UU);
RcppExport SEXP _QATS_QATS_timer_cpp(SEXP d0SEXP, SEXP n_seedsSEXP, SEXP rotateSEXP, SEXP n_repSEXP, SEXP nSEXP, SEXP mSEXP, SEXP logPiSEXP, SEXP qqSEXP, SEXP GGSEXP, SEXP SSSEXP, SEXP UUSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type d0(d0SEXP);
    Rcpp::traits::input_parameter< int >::type n_seeds(n_seedsSEXP);
    Rcpp::traits::input_parameter< bool >::type rotate(rotateSEXP);
    Rcpp::traits::input_parameter< double >::type n_rep(n_repSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type logPi(logPiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type qq(qqSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type GG(GGSEXP);
    Rcpp::traits::input_parameter< arma::ivec >::type SS(SSSEXP);
    Rcpp::traits::input_parameter< int >::type UU(UUSEXP);
    rcpp_result_gen = Rcpp::wrap(QATS_timer_cpp(d0, n_seeds, rotate, n_rep, n, m, logPi, qq, GG, SS, UU));
    return rcpp_result_gen;
END_RCPP
}
// Viterbi_timer_cpp
List Viterbi_timer_cpp(double n_rep, int n, int m, const arma::vec& logPi, const arma::mat& qq, const arma::mat& g_mseq);
RcppExport SEXP _QATS_Viterbi_timer_cpp(SEXP n_repSEXP, SEXP nSEXP, SEXP mSEXP, SEXP logPiSEXP, SEXP qqSEXP, SEXP g_mseqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type n_rep(n_repSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type logPi(logPiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type qq(qqSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type g_mseq(g_mseqSEXP);
    rcpp_result_gen = Rcpp::wrap(Viterbi_timer_cpp(n_rep, n, m, logPi, qq, g_mseq));
    return rcpp_result_gen;
END_RCPP
}
// Viterbi_cpp
void Viterbi_cpp(arma::ivec& xx, arma::imat& zeta, arma::mat& rho, const int& n, const int& m, const arma::vec& logPi, const arma::mat& qq, const arma::mat& g_mseq);
RcppExport SEXP _QATS_Viterbi_cpp(SEXP xxSEXP, SEXP zetaSEXP, SEXP rhoSEXP, SEXP nSEXP, SEXP mSEXP, SEXP logPiSEXP, SEXP qqSEXP, SEXP g_mseqSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::ivec& >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< arma::imat& >::type zeta(zetaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type logPi(logPiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type qq(qqSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type g_mseq(g_mseqSEXP);
    Viterbi_cpp(xx, zeta, rho, n, m, logPi, qq, g_mseq);
    return R_NilValue;
END_RCPP
}
// G_classifier_cpp
arma::vec G_classifier_cpp(double C1, double C2, double C3, double C4, int n, int m, const arma::vec& Pi, const arma::vec& logPi, const arma::mat& pp, const arma::mat& qq, const arma::mat& f_mseq, const arma::mat& g_mseq);
RcppExport SEXP _QATS_G_classifier_cpp(SEXP C1SEXP, SEXP C2SEXP, SEXP C3SEXP, SEXP C4SEXP, SEXP nSEXP, SEXP mSEXP, SEXP PiSEXP, SEXP logPiSEXP, SEXP ppSEXP, SEXP qqSEXP, SEXP f_mseqSEXP, SEXP g_mseqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type C1(C1SEXP);
    Rcpp::traits::input_parameter< double >::type C2(C2SEXP);
    Rcpp::traits::input_parameter< double >::type C3(C3SEXP);
    Rcpp::traits::input_parameter< double >::type C4(C4SEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type logPi(logPiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type pp(ppSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type qq(qqSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type f_mseq(f_mseqSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type g_mseq(g_mseqSEXP);
    rcpp_result_gen = Rcpp::wrap(G_classifier_cpp(C1, C2, C3, C4, n, m, Pi, logPi, pp, qq, f_mseq, g_mseq));
    return rcpp_result_gen;
END_RCPP
}
// PMAP_timer_cpp
List PMAP_timer_cpp(double n_rep, int n, int m, const arma::vec& Pi, const arma::mat& pp, const arma::mat& f_mseq);
RcppExport SEXP _QATS_PMAP_timer_cpp(SEXP n_repSEXP, SEXP nSEXP, SEXP mSEXP, SEXP PiSEXP, SEXP ppSEXP, SEXP f_mseqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type n_rep(n_repSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type pp(ppSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type f_mseq(f_mseqSEXP);
    rcpp_result_gen = Rcpp::wrap(PMAP_timer_cpp(n_rep, n, m, Pi, pp, f_mseq));
    return rcpp_result_gen;
END_RCPP
}
// PMAP_cpp
void PMAP_cpp(arma::ivec& xx, arma::mat& alpha_hat, arma::mat& alpha_bar, arma::mat& beta_hat, arma::mat& beta_bar, arma::vec& cc_inv, const int& n, const int& m, const arma::vec& Pi, const arma::mat& pp, const arma::mat& f_mseq);
RcppExport SEXP _QATS_PMAP_cpp(SEXP xxSEXP, SEXP alpha_hatSEXP, SEXP alpha_barSEXP, SEXP beta_hatSEXP, SEXP beta_barSEXP, SEXP cc_invSEXP, SEXP nSEXP, SEXP mSEXP, SEXP PiSEXP, SEXP ppSEXP, SEXP f_mseqSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::ivec& >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type alpha_hat(alpha_hatSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type alpha_bar(alpha_barSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type beta_hat(beta_hatSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type beta_bar(beta_barSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type cc_inv(cc_invSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type pp(ppSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type f_mseq(f_mseqSEXP);
    PMAP_cpp(xx, alpha_hat, alpha_bar, beta_hat, beta_bar, cc_inv, n, m, Pi, pp, f_mseq);
    return R_NilValue;
END_RCPP
}
// K_segmentation_cpp
arma::mat K_segmentation_cpp(int K_max, int n, int m, const arma::vec& logPi, const arma::mat& qq, const arma::mat& g_mseq);
RcppExport SEXP _QATS_K_segmentation_cpp(SEXP K_maxSEXP, SEXP nSEXP, SEXP mSEXP, SEXP logPiSEXP, SEXP qqSEXP, SEXP g_mseqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type K_max(K_maxSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type logPi(logPiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type qq(qqSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type g_mseq(g_mseqSEXP);
    rcpp_result_gen = Rcpp::wrap(K_segmentation_cpp(K_max, n, m, logPi, qq, g_mseq));
    return rcpp_result_gen;
END_RCPP
}
// H1_dbl_cpp
double H1_dbl_cpp(int l, int r, int x0, int m, const arma::vec& logPi, const arma::mat& qq, const arma::mat& GG);
RcppExport SEXP _QATS_H1_dbl_cpp(SEXP lSEXP, SEXP rSEXP, SEXP x0SEXP, SEXP mSEXP, SEXP logPiSEXP, SEXP qqSEXP, SEXP GGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type logPi(logPiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type qq(qqSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type GG(GGSEXP);
    rcpp_result_gen = Rcpp::wrap(H1_dbl_cpp(l, r, x0, m, logPi, qq, GG));
    return rcpp_result_gen;
END_RCPP
}
// H2_vec_cpp
arma::vec H2_vec_cpp(int l, int r, int x0, int m, const arma::vec& logPi, const arma::mat& qq, const arma::mat& GG);
RcppExport SEXP _QATS_H2_vec_cpp(SEXP lSEXP, SEXP rSEXP, SEXP x0SEXP, SEXP mSEXP, SEXP logPiSEXP, SEXP qqSEXP, SEXP GGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type logPi(logPiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type qq(qqSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type GG(GGSEXP);
    rcpp_result_gen = Rcpp::wrap(H2_vec_cpp(l, r, x0, m, logPi, qq, GG));
    return rcpp_result_gen;
END_RCPP
}
// H3_mat_cpp
arma::mat H3_mat_cpp(int l, int r, int x0, int m, const arma::vec& logPi, const arma::mat& qq, const arma::mat& GG);
RcppExport SEXP _QATS_H3_mat_cpp(SEXP lSEXP, SEXP rSEXP, SEXP x0SEXP, SEXP mSEXP, SEXP logPiSEXP, SEXP qqSEXP, SEXP GGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type logPi(logPiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type qq(qqSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type GG(GGSEXP);
    rcpp_result_gen = Rcpp::wrap(H3_mat_cpp(l, r, x0, m, logPi, qq, GG));
    return rcpp_result_gen;
END_RCPP
}
// sample_norm_HMM_export_cpp
List sample_norm_HMM_export_cpp(int n, int m, const arma::vec& Pi, const arma::mat& pp, double sigma);
RcppExport SEXP _QATS_sample_norm_HMM_export_cpp(SEXP nSEXP, SEXP mSEXP, SEXP PiSEXP, SEXP ppSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type pp(ppSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_norm_HMM_export_cpp(n, m, Pi, pp, sigma));
    return rcpp_result_gen;
END_RCPP
}
// QATS_nseeds_norm_cpp
arma::mat QATS_nseeds_norm_cpp(int n, int m, const arma::vec& Pi, const arma::mat& pp, double sigma, int d0, arma::ivec n_seeds, bool rotate, int n_rep, int n_sim);
RcppExport SEXP _QATS_QATS_nseeds_norm_cpp(SEXP nSEXP, SEXP mSEXP, SEXP PiSEXP, SEXP ppSEXP, SEXP sigmaSEXP, SEXP d0SEXP, SEXP n_seedsSEXP, SEXP rotateSEXP, SEXP n_repSEXP, SEXP n_simSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type pp(ppSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type d0(d0SEXP);
    Rcpp::traits::input_parameter< arma::ivec >::type n_seeds(n_seedsSEXP);
    Rcpp::traits::input_parameter< bool >::type rotate(rotateSEXP);
    Rcpp::traits::input_parameter< int >::type n_rep(n_repSEXP);
    Rcpp::traits::input_parameter< int >::type n_sim(n_simSEXP);
    rcpp_result_gen = Rcpp::wrap(QATS_nseeds_norm_cpp(n, m, Pi, pp, sigma, d0, n_seeds, rotate, n_rep, n_sim));
    return rcpp_result_gen;
END_RCPP
}
// QATS_vs_Viterbi_norm_cpp
List QATS_vs_Viterbi_norm_cpp(int n, int m, const arma::vec& Pi, const arma::mat& pp, double sigma, int d0, int n_seeds, bool rotate, int n_rep, int n_sim);
RcppExport SEXP _QATS_QATS_vs_Viterbi_norm_cpp(SEXP nSEXP, SEXP mSEXP, SEXP PiSEXP, SEXP ppSEXP, SEXP sigmaSEXP, SEXP d0SEXP, SEXP n_seedsSEXP, SEXP rotateSEXP, SEXP n_repSEXP, SEXP n_simSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type pp(ppSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type d0(d0SEXP);
    Rcpp::traits::input_parameter< int >::type n_seeds(n_seedsSEXP);
    Rcpp::traits::input_parameter< bool >::type rotate(rotateSEXP);
    Rcpp::traits::input_parameter< int >::type n_rep(n_repSEXP);
    Rcpp::traits::input_parameter< int >::type n_sim(n_simSEXP);
    rcpp_result_gen = Rcpp::wrap(QATS_vs_Viterbi_norm_cpp(n, m, Pi, pp, sigma, d0, n_seeds, rotate, n_rep, n_sim));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_QATS_QATS_timer_cpp", (DL_FUNC) &_QATS_QATS_timer_cpp, 11},
    {"_QATS_Viterbi_timer_cpp", (DL_FUNC) &_QATS_Viterbi_timer_cpp, 6},
    {"_QATS_Viterbi_cpp", (DL_FUNC) &_QATS_Viterbi_cpp, 8},
    {"_QATS_G_classifier_cpp", (DL_FUNC) &_QATS_G_classifier_cpp, 12},
    {"_QATS_PMAP_timer_cpp", (DL_FUNC) &_QATS_PMAP_timer_cpp, 6},
    {"_QATS_PMAP_cpp", (DL_FUNC) &_QATS_PMAP_cpp, 11},
    {"_QATS_K_segmentation_cpp", (DL_FUNC) &_QATS_K_segmentation_cpp, 6},
    {"_QATS_H1_dbl_cpp", (DL_FUNC) &_QATS_H1_dbl_cpp, 7},
    {"_QATS_H2_vec_cpp", (DL_FUNC) &_QATS_H2_vec_cpp, 7},
    {"_QATS_H3_mat_cpp", (DL_FUNC) &_QATS_H3_mat_cpp, 7},
    {"_QATS_sample_norm_HMM_export_cpp", (DL_FUNC) &_QATS_sample_norm_HMM_export_cpp, 5},
    {"_QATS_QATS_nseeds_norm_cpp", (DL_FUNC) &_QATS_QATS_nseeds_norm_cpp, 10},
    {"_QATS_QATS_vs_Viterbi_norm_cpp", (DL_FUNC) &_QATS_QATS_vs_Viterbi_norm_cpp, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_QATS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
