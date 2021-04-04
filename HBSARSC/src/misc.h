#include <Rmath.h>
#include <RcppArmadillo.h>

#ifndef MISC
#define MISC

arma::colvec rcpp_CumTrap(arma::colvec f, double delta);
double rcpp_trapint(arma::colvec f, double xdelta);
arma::colvec rcpp_GetSquish(double omega, double psi, arma::colvec xgrid);

double rcpp_rndtna(double mu, double sigma, double xtop);
double rcpp_rndtnb(double mu, double sigma, double xbot);
arma::colvec rcpp_rndtnb2(arma::colvec mu, double sigma, double xbot);
double rcpp_rndtnab(double mu, double sigma, double xbot, double xtop);
double rcpp_rndtnab(arma::colvec& mu, double sigma, double xbot, double xtop);
arma::colvec rcpp_rndtnab(arma::colvec& mu, arma::mat& sigma, double xbot, double xtop);

arma::colvec rcpp_Getfx(Rcpp::List fpars);
arma::colvec rcpp_GetMonofxgrid(Rcpp::List fpars);
arma::colvec rcpp_GetMonoSamefxgrid(Rcpp::List fpars);
arma::colvec rcpp_GetMonoDiffxgrid(Rcpp::List fpars);
arma::colvec rcpp_GetUfxgrid(Rcpp::List fpars);
arma::colvec rcpp_GetSfxgrid(Rcpp::List fpars, int d3);

arma::colvec rcpp_Getf0x(Rcpp::List f0pars);
arma::colvec rcpp_GetMonof0xgrid(Rcpp::List f0pars);
arma::colvec rcpp_GetMonoSamef0xgrid(Rcpp::List f0pars);
arma::colvec rcpp_GetMonoDiff0xgrid(Rcpp::List f0pars);
arma::colvec rcpp_GetUf0xgrid(Rcpp::List f0pars);
arma::colvec rcpp_GetSf0xgrid(Rcpp::List f0pars, int d3);


Rcpp::List rcpp_GetUpfxgrid(arma::colvec& fpara, arma::mat& phix, double xdelta, double xrange,
                            bool centering = true);
arma::colvec rcpp_GetSCfxobs(arma::colvec& xgrid, arma::uvec& xinx, arma::vec& xover, 
                             arma::colvec& fxgrid, int iflagpn);

arma::colvec rcpp_GetUpf0xgrid(arma::colvec& fpara, arma::colvec& vpara, arma::mat& phix, 
                               arma::mat& phix2, double xdelta, double xrange, bool centering=true);
arma::colvec rcpp_GetUpf0Biasxgrid(arma::colvec& fpara,arma::mat& phix, double xdelta, double xrange,
                                   bool centering=true);

arma::colvec rcpp_GetMonof(arma::colvec& theta, arma::mat& phixobs, arma::umat& uquadfacts, arma::mat& quadfacts, int fpm);
arma::colvec rcpp_QuadMult(arma::colvec& x, arma::mat& qvech, arma::umat& uquadfacts, arma::mat& quadfacts);

arma::colvec rcpp_cholSol (arma::colvec& b, arma::mat& U);

// vcbsar
arma::mat getZG_univ (arma::cube& zCube, arma::mat& vecXiOrdered);
arma::mat updateXi_univ_ind_simple (arma::vec& yResid, double Sigma, arma::mat& zMat, arma::cube& ztzs, arma::mat& sigRi, arma::uvec& id, arma::uvec& nivec);


// basis
arma::colvec IntCos2 (double x, arma::colvec kall, double xmin, double xrange, int nbasis);
double IntCosCrossProd (double x, int j, int k, double xmin, double xrange);
double IntConst2 (double x, double xmin, double xrange);
arma::colvec IntCos (double x, arma::colvec kall, double xmin, double xrange, int nbasis);
arma::colvec vech (arma::mat A);
arma::mat computeBasis (arma::colvec xobs, double xmin,double xrange, unsigned int nobs, unsigned int nbasis, unsigned int nr);

void update_lambda(arma::mat& lambda, arma::mat& lambdai, arma::mat resid_mat, arma::mat lambda_g0i, 
                   double lambda_fn, unsigned int nparw);
arma::mat update_phi(arma::mat lambdai, arma::mat ztz, arma::mat phi_v0i, arma::mat zdata, arma::mat betall,
                     unsigned int nparw, unsigned int nparz, unsigned int phidim);

// metropolis
arma::colvec GetMetVec(double w, double m0, double alpha, unsigned int ndim);
void SetMetVec(arma::colvec &Met, double w, double m0, double alpha, unsigned int ndim);
arma::colvec AdaptMetVar(arma::colvec x);
void AdaptMetUpdate(arma::colvec &x);
arma::colvec AdaptMetUpdate2(arma::colvec x);
void UpdateMet(arma::colvec &x, unsigned int nb, double fct);
arma::colvec UpdateMet2(arma::colvec x, unsigned int nb, double fct);

#endif
