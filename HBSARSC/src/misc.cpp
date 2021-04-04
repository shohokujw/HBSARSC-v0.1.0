#include "misc.h"
#define pi 3.141592653589793238462643383279502884197

using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

// computes cumulative integral of function
// [[Rcpp::export]]
arma::colvec rcpp_CumTrap(arma::colvec f, double delta) {
  unsigned int nr = f.n_elem;
  colvec fcum = arma::zeros<colvec> (nr);
  fcum.subvec(1,nr-1) = cumsum(f.subvec(0,nr-2)) + cumsum(f.subvec(1,nr-1));
  fcum *= (delta/2); 
  return fcum;
}


// trapidzoid integration
// [[Rcpp::export]]
double rcpp_trapint(arma::colvec f, double xdelta) {
  unsigned int nr = f.n_elem;
  double a = sum(f);
  a -= ((f[0]+f[nr-1])/2);
  return a*xdelta;
}


// Compute squish function
// [[Rcpp::export]]
arma::colvec rcpp_GetSquish(double omega, double psi, arma::colvec xgrid) {
  colvec c = psi*(xgrid-omega);
  colvec squish = zeros<colvec>(c.n_elem);
  // Worry about underflow and over flow
  // if c < -100, then hfun = 1
  // if c >  100, then hfun = -1
  uvec idz1 = find(c < -100);
  uvec idz2 = find(c >  100);
  
  if (idz1.n_elem > 0) {
    c.elem( idz1 ) = -100*ones<vec>(idz1.n_elem);  
  }
  
  if (idz2.n_elem > 0) {
    c.elem( idz2 ) = 100*ones<vec>(idz2.n_elem);  
  }
  
  c = exp(c);
  squish = (1-c)/(1+c);
  
  return squish;
}


// return -idx
// [[Rcpp::export]]
arma::uvec get_ex_idx(unsigned int len, arma::uvec idx) {
  if (len <= idx.n_elem) {
    return NULL;
  }
  
  uvec v = zeros<uvec>(len);
  for (unsigned int i=0; i<idx.n_elem; i++) {
    v(idx(i)) = 1;
  }
  
  uvec ex_vec (len - idx.n_elem);
  int j=0;
  for (unsigned int i=0; i<len; i++) {
    if (v(i) == 0) {
      ex_vec(j) = i;
      j++;
    }
  }
  
  return ex_vec;
}


// generates form a truncated above normal distribution
// [[Rcpp::export]]
double rcpp_rndtna(double mu, double sigma, double xtop) {
  double u = randu();
  double fa = 0.0;
  double fb = normcdf(xtop, mu, sigma);
  double p = fa + u*(fb-fa);
  // Worry about p = 0 or p = 1
  if (p < 0.00001) p = 0.00001;
  if (p > 0.99999) p = 0.99999;
  double x = R::qnorm(p, mu, sigma, true, false);
  
  if (x > xtop) x = xtop;
  
  return x;
}


// generates form a blow normal distribution
// [[Rcpp::export]]
double rcpp_rndtnb(double mu, double sigma, double xbot) {
  double u = arma::randu();
  double fa = arma::normcdf(xbot, mu, sigma);
  double fb = 1.0;
  double p = fa + u*(fb-fa);
  // Worry about p = 0 or p = 1
  if (p < 0.00001) p = 0.00001;
  if (p > 0.99999) p = 0.99999;
  double x = R::qnorm(p, mu, sigma, true, false);
  
  if (x < xbot) x = xbot;
  
  return x;
}


// generates form a blow normal distribution
// [[Rcpp::export]]
arma::colvec rcpp_rndtnb2(arma::colvec mu, double sigma, double xbot) {
  unsigned int n  = mu.n_elem;
  colvec fa = zeros<colvec> (n);
  colvec x  = zeros<colvec> (n);
  
  arma::colvec u = arma::randu(n);
  
  for (unsigned int i=0; i<n; i++) {
    double fa = normcdf(xbot, mu(i), sigma);
    double fb = 1.0;
    double p = fa + u(i)*(fb-fa);
    // Worry about p = 0 or p = 1
    if (p < 0.00001) p = 0.00001;
    if (p > 0.99999) p = 0.99999;
    x(i) = R::qnorm(p, mu(i), sigma, true, false);
    
    if (x(i) < xbot) x(i) = xbot;
    
  }
  
  
  return x;
}


// generates form a truncated normal distribution
// [[Rcpp::export]]
double rcpp_rndtnab(double mu, double sigma, double xbot, double xtop) {
  double u = randu();
  double fa = normcdf(xbot, mu, sigma);
  double fb = normcdf(xtop, mu, sigma);
  // Need to gaurd against bad falues for fa and fb
  double x0 = (xbot+xtop)/2; // Good guess that stasifies xbot < x < xtop
  double delta = fa < fb ? fb-fa : 0;
  double p = fa + u*delta;
  // Worry about p = 0 or p = 1
  if (p < 0.00001) p = 0.00001;
  if (p > 0.99999) p = 0.99999;
  double x = R::qnorm(p, mu, sigma, true, false);
  if (fa >= fb) x = x0;
  if (x < xbot) x = xbot;
  if (x > xtop) x = xtop;
  
  return x;
}

// [[Rcpp::export]]
arma::colvec rcpp_Getf0x(Rcpp::List f0pars) {
  colvec fx;
  
  int fshape = f0pars["iflagsc"];
  
  switch (fshape)
  {
  case 0: // free function
  {
    mat phix     = f0pars["phix"];
    colvec theta = f0pars["theta"];
    fx = phix * theta;
    break;
  }
  case 1: // monotone function
  {
    fx = rcpp_GetMonof0xgrid(f0pars);
    break;
  }
  case 2: // Increasing and convex or decreasing and concave
  {
    fx = rcpp_GetMonoSamef0xgrid(f0pars);
    break;
  }
  case 3: // Increasing and concave or decreasing and convex
  {
    fx = rcpp_GetMonoDiff0xgrid(f0pars);
    break;
  }
  case 4: // U-shaped function
  {
    fx = rcpp_GetUf0xgrid(f0pars);
    break;
  }
  case 5: // S-shaped (increasing,convex-to-concave) or (decreasing,concave-to-convex)
  {
    fx = rcpp_GetSf0xgrid(f0pars,1);
    break;
  }
  case 6: // S-shaped (increasing,concave-to-convex) or (decreasing,convex-to-concave)
  {
    fx = rcpp_GetSf0xgrid(f0pars,-1);
    break;
  }
  default:
    std::cout << "not available function type" << std::endl;
    break;
  }
  
  return fx;
}


// [[Rcpp::export]]
arma::colvec rcpp_Getfx(Rcpp::List fpars) {
  colvec fx;
  
  int fshape = fpars["iflagsc"];
  
  switch (fshape)
  {
  case 0: // free function
  {
    mat phix     = fpars["phix"];
    colvec theta = fpars["theta"];
    fx = phix * theta;
    break;
  }
  case 1: // monotone function
  {  
    fx = rcpp_GetMonofxgrid(fpars);
    break;
  }
  case 2: // Increasing convex or decreasing and concave
  {  
    fx = rcpp_GetMonoSamefxgrid(fpars);
    break;
  }
  case 3: // Increasing and concave or decreasing and convex
  {
    fx = rcpp_GetMonoDiffxgrid(fpars);
    break;
  }
  case 4: // U-shaped function
  {
    fx = rcpp_GetUfxgrid(fpars);
    break;
  }
  case 5: // S-shaped (increasing,convex-to-concave) or (decreasing,concave-to-convex)
  {
    fx = rcpp_GetSfxgrid(fpars,1);
    break;
  }
  case 6: // S-shaped (increasing,concave-to-convex) or (decreasing,convex-to-concave)
  {
    fx = rcpp_GetSfxgrid(fpars,-1);
    break;
  }
  default:
    std::cout << "not available function type" << std::endl;
  break;
  }
  
  return fx;
}


// [[Rcpp::export]]
arma::colvec rcpp_GetMonofxgrid(Rcpp::List fpars) {
  mat phix         = fpars["phix"];
  colvec theta     = fpars["theta"];
  int iflagpn      = fpars["iflagpn"];
  double delta     = fpars["delta"];
  bool iflagCenter = fpars["iflagCenter"];
  double range     = fpars["range"];
  bool iflagZ      = fpars["iflagZ"];
  colvec dfx;
  
  colvec z  = phix * theta;
  if (iflagZ) {
    dfx = exp(z);
  } else {
    dfx = square(z);
  }
  colvec fx = iflagpn*rcpp_CumTrap(dfx,delta);
  if(iflagCenter) { // center f
    // Compute mean
    double c = rcpp_trapint(fx,delta);
    fx -= c/range;
  }
  
  return fx;
}


// [[Rcpp::export]]
arma::colvec rcpp_GetMonoSamefxgrid(Rcpp::List fpars) {
  mat phix         = fpars["phix"];
  colvec theta     = fpars["theta"];
  int iflagpn      = fpars["iflagpn"];
  double delta     = fpars["delta"];
  bool iflagCenter = fpars["iflagCenter"];
  double range     = fpars["range"];
  bool iflagZ      = fpars["iflagZ"];
  colvec d2fx;
  
  colvec z  = phix * theta;
  if (iflagZ) {
    d2fx = exp(z);
  } else {
    d2fx = square(z);
  }
  colvec dfx = iflagpn*rcpp_CumTrap(d2fx,delta); 
  colvec fx  = rcpp_CumTrap(dfx, delta); 
  if(iflagCenter) { // center f
    // Compute mean
    double c = rcpp_trapint(fx,delta);
    fx -= c/range;
  }
  
  return fx;
}


// [[Rcpp::export]]
arma::colvec rcpp_GetMonoDiffxgrid(Rcpp::List fpars) {
  mat phix         = fpars["phix"];
  colvec theta     = fpars["theta"];
  int iflagpn      = fpars["iflagpn"];
  double delta     = fpars["delta"];
  bool iflagCenter = fpars["iflagCenter"];
  double range     = fpars["range"];
  bool iflagZ      = fpars["iflagZ"];
  colvec d2fx;
  
  colvec z  = phix * theta;
  if (iflagZ) {
    d2fx = exp(z);
  } else {
    d2fx = square(z);
  }
  colvec dfx = iflagpn*rcpp_CumTrap(d2fx,delta); 
  colvec fx  = rcpp_CumTrap(reverse(dfx), delta); 
  if(iflagCenter) { // center f
    // Compute mean
    double c = rcpp_trapint(fx,delta);
    fx -= c/range;
  }
  
  return fx;
}


// [[Rcpp::export]]
arma::colvec rcpp_GetUfxgrid(Rcpp::List fpars) {
  mat phix         = fpars["phix"];
  colvec theta     = fpars["theta"];
  colvec hfun      = fpars["hfun"];
  int iflagpn      = fpars["iflagpn"];
  double delta     = fpars["delta"];
  bool iflagCenter = fpars["iflagCenter"];
  double range     = fpars["range"];
  bool iflagZ      = fpars["iflagZ"];
  colvec dfx;
  
  colvec z  = phix * theta;
  if (iflagZ) {
    dfx = exp(z);
  } else {
    dfx = square(z);
  }
  dfx = hfun%dfx;
  colvec fx = iflagpn * rcpp_CumTrap(dfx,delta);
  
  if(iflagCenter) { // center f
    // Compute mean 
    double c = rcpp_trapint(fx,delta);
    fx -= c/range;
  }
  
  return fx;
}


// [[Rcpp::export]]
arma::colvec rcpp_GetSfxgrid(Rcpp::List fpars, int d3) {
  mat phix         = fpars["phix"];
  colvec theta     = fpars["theta"];
  colvec hfun      = fpars["hfun"];
  int iflagpn      = fpars["iflagpn"];
  double delta     = fpars["delta"];
  bool iflagCenter = fpars["iflagCenter"];
  double range     = fpars["range"];
  bool iflagZ      = fpars["iflagZ"];
  colvec xgrid     = fpars["xgrid"];
  double xmin      = fpars["xmin"];
  colvec d2fx;
  colvec dfx;
  colvec fx;
  double c2;
  
  colvec z  = phix * theta;
  if (iflagZ) {
    d2fx = exp(z);
  } else {
    d2fx = square(z);
  }
  //d3 =  1 -> hfun is "squish down" for "S" 
  // d3 = -1 -> hfun is "squish up" for "rotated S"
  d2fx = d3*hfun%d2fx;
  dfx  = rcpp_CumTrap(d2fx, delta);
  c2   = std::min(0.0, dfx.min());
  fx   = rcpp_CumTrap(dfx,delta);
  fx   = iflagpn*(fx - c2*(xgrid-xmin));
  
  if(iflagCenter) { // center f
    // Compute mean 
    double c = rcpp_trapint(fx,delta);
    fx -= c/range;
  }
  
  return fx;
}


// [[Rcpp::export]]
arma::colvec rcpp_GetMonof0xgrid(Rcpp::List f0pars) {
  mat phix         = f0pars["phix"];
  colvec theta     = f0pars["theta"];
  mat phix2        = f0pars["phix2"];
  colvec vpara     = f0pars["vpara"];
  int iflagpn      = f0pars["iflagpn"];
  double delta     = f0pars["delta"];
  bool iflagCenter = f0pars["iflagCenter"];
  double range     = f0pars["range"];
  bool iflagZ      = f0pars["iflagZ"];
  colvec dfx;
  
  colvec z   = phix * theta;
  colvec z2  = square(z);
  colvec v2  = phix2 * vpara;
  if (iflagZ) {
    dfx = exp(z + v2/2);
  } else {
    dfx = z2 + v2;
  }
  
  colvec fx  = iflagpn*rcpp_CumTrap(dfx,delta);
  if(iflagCenter){ // center f
    // Compute mean
    double c = rcpp_trapint(fx,delta);
    fx -= c/range;
  }
  
  return fx;
}


// [[Rcpp::export]]
arma::colvec rcpp_GetMonoSamef0xgrid(Rcpp::List f0pars) {
  mat phix         = f0pars["phix"];
  colvec theta     = f0pars["theta"];
  mat phix2        = f0pars["phix2"];
  colvec vpara     = f0pars["vpara"];
  int iflagpn      = f0pars["iflagpn"];
  double delta     = f0pars["delta"];
  bool iflagCenter = f0pars["iflagCenter"];
  double range     = f0pars["range"];
  bool iflagZ      = f0pars["iflagZ"];
  colvec d2fx;
  
  colvec z  = phix * theta;
  colvec z2 = square(z);
  colvec v2 = phix2 * vpara;
  if (iflagZ) {
    d2fx = exp(z + v2/2);
  } else {
    d2fx = z2 + v2;
  }
  colvec dfx = iflagpn*rcpp_CumTrap(d2fx,delta);
  colvec fx  = rcpp_CumTrap(dfx,delta);
  if(iflagCenter){ // center f
    // Compute mean
    double c = rcpp_trapint(fx,delta);
    fx -= c/range;
  }
  
  return fx;
}


// [[Rcpp::export]]
arma::colvec rcpp_GetMonoDiff0xgrid(Rcpp::List f0pars) {
  mat phix         = f0pars["phix"];
  colvec theta     = f0pars["theta"];
  mat phix2        = f0pars["phix2"];
  colvec vpara     = f0pars["vpara"];
  int iflagpn      = f0pars["iflagpn"];
  double delta     = f0pars["delta"];
  bool iflagCenter = f0pars["iflagCenter"];
  double range     = f0pars["range"];
  bool iflagZ      = f0pars["iflagZ"];
  colvec d2fx;
  
  colvec z  = phix * theta;
  colvec z2 = square(z);
  colvec v2 = phix2 * vpara;
  if (iflagZ) {
    d2fx = exp(z + v2/2);
  } else {
    d2fx = z2 + v2;
  }
  colvec dfx = iflagpn*rcpp_CumTrap(d2fx,delta);
  colvec fx  = rcpp_CumTrap(reverse(dfx),delta);
  if(iflagCenter){ // center f
    // Compute mean
    double c = rcpp_trapint(fx,delta);
    fx -= c/range;
  }
  
  return fx;
}


// [[Rcpp::export]]
arma::colvec rcpp_GetUf0xgrid(Rcpp::List f0pars) {
  mat phix         = f0pars["phix"];
  colvec theta     = f0pars["theta"];
  mat phix2        = f0pars["phix2"];
  colvec vpara     = f0pars["vpara"];
  colvec hfun      = f0pars["hfun"];
  int iflagpn      = f0pars["iflagpn"];
  double delta     = f0pars["delta"];
  bool iflagCenter = f0pars["iflagCenter"];
  double range     = f0pars["range"];
  bool iflagZ      = f0pars["iflagZ"];
  colvec dfx;
  
  colvec z  = phix * theta;
  colvec z2 = square(z);
  colvec v2 = phix2 * vpara;
  if (iflagZ) {
    dfx = exp(z + v2/2);
  } else {
    dfx = z2 + v2;
  }
  dfx = hfun%dfx;
  
  colvec fx = iflagpn*rcpp_CumTrap(dfx,delta);
  if(iflagCenter){ // center f
    // Compute mean
    double c = rcpp_trapint(fx,delta);
    fx -= c/range;
  }
  
  return fx;
}


// [[Rcpp::export]]
arma::colvec rcpp_GetSf0xgrid(Rcpp::List f0pars, int d3) {
  mat phix         = f0pars["phix"];
  colvec theta     = f0pars["theta"];
  mat phix2        = f0pars["phix2"];
  colvec vpara     = f0pars["vpara"];
  colvec hfun      = f0pars["hfun"];
  int iflagpn      = f0pars["iflagpn"];
  double delta     = f0pars["delta"];
  bool iflagCenter = f0pars["iflagCenter"];
  double range     = f0pars["range"];
  bool iflagZ      = f0pars["iflagZ"];
  colvec xgrid     = f0pars["xgrid"];
  double xmin      = f0pars["xmin"];
  colvec d2fx;
  
  colvec z  = phix * theta;
  colvec z2 = square(z);
  colvec v2 = phix2 * vpara;
  if (iflagZ) {
    d2fx = exp(z + v2/2);
  } else {
    d2fx = z2 + v2;
  }
  d2fx        = d3*hfun%d2fx;
  // d3 =  1 -> hfun is "squish down"
  // d3 = -1 -> hfun is "squish up"
  colvec dfx = rcpp_CumTrap(d2fx, delta);
  double c2  = std::min(0.0, dfx.min());
  colvec fx  = rcpp_CumTrap(dfx, delta);
  fx         = iflagpn*(fx - c2*(xgrid-xmin));
  if(iflagCenter){ // center f
    // Compute mean
    double c = rcpp_trapint(fx,delta);
    fx -= c/range;
  }
  
  return fx;
}


// [[Rcpp::export]]
Rcpp::List rcpp_GetUpfxgrid(arma::colvec& fpara, arma::mat& phix, double xdelta, double xrange, bool centering) {
  colvec z = phix * fpara;
  colvec dfx = z%z;
  colvec fx = rcpp_CumTrap(dfx, xdelta);
  if (centering) { // Center f(x): int f(x) = 0
    // Compute mean
    double c = rcpp_trapint(fx, xdelta);  
    fx -= (c/xrange);
  }
  
  return Rcpp::List::create(Rcpp::Named("fx")=fx,
                            Rcpp::Named("dfx")=dfx);
}


// Compute f at observations
// [[Rcpp::export]]
arma::colvec rcpp_GetSCfxobs(arma::colvec& xgrid, arma::uvec& xinx, arma::vec& xover, arma::colvec& fxgrid, int iflagpn) {
  double xdelta  = xgrid(1) - xgrid(0);
  unsigned int n = xinx.n_elem;
  colvec fxobs   = zeros<colvec> (n);
  uvec idn0      = find(xover > 0);
  if (idn0.n_elem < n) { // Got data on grid points
    uvec idx = get_ex_idx(n, idn0);
    fxobs.elem(idx) = fxgrid.elem(xinx.elem(idx));
  } 
  // overshoot
  colvec y0 = fxgrid.elem(xinx.elem(idn0));
  colvec y1 = fxgrid.elem((xinx.elem(idn0)+1));
  fxobs.elem(idn0) = y0 + (y1-y0)%xover.elem(idn0)/xdelta;  // Linearly interpolate  
  
  fxobs = iflagpn*fxobs;
  
  return fxobs;
}



// Compute increasing f = int_a^x Z(s)^2 ds on xgrid
// [[Rcpp::export]]
arma::colvec rcpp_GetUpf0xgrid(arma::colvec& fpara, arma::colvec& vpara, arma::mat& phix, arma::mat& phix2,
                               double xdelta, double xrange, bool centering) {
  colvec z  = phix * fpara;
  colvec v  = phix2 * vpara;
  colvec fx = rcpp_CumTrap(square(z)+v, xdelta);
  if (!centering) {  // do not center f
    return fx;
  } else {  // center f
    // Compute mean
    return fx - rcpp_trapint(fx,xdelta)/xrange;
  }
}


// Compute increasing f = int_a^x Z(s)^2 ds on xgrid
// [[Rcpp::export]]
arma::colvec rcpp_GetUpf0Biasxgrid(arma::colvec& fpara, arma::mat& phix, double xdelta, double xrange, 
                                   bool centering) {
  colvec z   = phix * fpara;
  colvec zv2 = square(z);
  colvec fx  = rcpp_CumTrap(zv2,xdelta); // v is not added so biased!
  if (!centering) {  // do not center f
    return fx;
  }else{  // center f
    // Compute mean
    return fx - rcpp_trapint(fx,xdelta)/xrange;
  }
}


// [[Rcpp::export]]
arma::colvec rcpp_GetMonof(arma::colvec& theta, arma::mat& phixobs, arma::umat& uquadfacts, arma::mat& quadfacts, int fpm){
  return(fpm*rcpp_QuadMult(theta,phixobs,uquadfacts,quadfacts));
}


// [[Rcpp::export]]
arma::colvec rcpp_QuadMult(arma::colvec& x, arma::mat& qvech, arma::umat& uquadfacts, arma::mat& quadfacts) {
  int nr = qvech.n_rows;
  int nc = qvech.n_cols;
  
  mat tmp = zeros<mat>(nr, nc);
  for (int j=0; j<nc; j++) {
    uvec vec1 = uquadfacts.col(1);
    uvec vec2 = uquadfacts.col(2);
    tmp.col(j) = quadfacts.col(0) % x.elem(vec1) % qvech.col(j) % x.elem(vec2);
  }
  colvec xQx = sum(tmp, 0).t();
  return xQx;
}


// [[Rcpp::export]]
arma::colvec rcpp_cholSol (arma::colvec& b, arma::mat& U) {
  return (solve(trimatu(U), solve(trimatl(U.t()), b)));
}


// [[Rcpp::export]]
arma::colvec IntCos2 (double x, arma::colvec kall, double xmin, double xrange, int nbasis) {
  colvec xout(nbasis);
  double z = (x-xmin)/xrange;
  xout = sin(kall*2.0*pi*z)/(2.0*pi*kall)+z-1.0/2.0;
  return xout;
}


// [[Rcpp::export]]
double IntCosCrossProd (double x, int j, int k, double xmin, double xrange) {
  double z = (x-xmin)/xrange;
  double xout = sin((j+k)*pi*z)/(pi*(j+k))+sin((k-j)*pi*z)/(pi*(k-j))-(1.0-cos(pi*(j+k)))/(pow(pi*(j+k),2))-(1.0-cos(pi*(j-k)))/(pow(pi*(j-k),2));
  return xout;
}


// [[Rcpp::export]]
double IntConst2 (double x, double xmin, double xrange) {
  double xout = (x-xmin)/xrange-1.0/2.0;
  return xout;
}


// [[Rcpp::export]]
arma::colvec IntCos (double x, arma::colvec kall, double xmin, double xrange, int nbasis) {
  colvec xout(nbasis);
  double z = (x-xmin)/xrange;
  
  xout=sqrt(2.0)*sin(kall*pi*z)/(pi*kall)- sqrt(2.0)*(1.0-cos(pi*kall))/(pow(pi*kall,2));
  
  return xout;
}


// [[Rcpp::export]]
arma::colvec vech (arma::mat A) {
  unsigned int nr = A.n_rows;
  unsigned int nc = A.n_cols;
  colvec vec(nr * (nc + 1)/2, fill::zeros);
  
  unsigned int k = 0;
  for (unsigned int ir=0; ir<nr; ir++) {
    for (unsigned int ic=0; ic<ir+1; ic++) {
      vec(k) = A(ir,ic);
      k++;
    }  
  }
  
  return vec;
}


// GetPhi(xobs,xmin,xrange,IntCos2,IntCosCrossProd,IntConst2,IntCos,nbasis,nobs,phixobs)
// [[Rcpp::export]]
arma::mat computeBasis (arma::colvec xobs, double xmin,double xrange, unsigned int nobs, unsigned int nbasis, unsigned int nr) {
  mat phixobs(nr,nobs, fill::zeros);
  colvec kall = regspace<colvec>(1,  1,  nbasis);
  
  double xi, b, c;
  colvec a(nbasis);
  colvec d(nbasis);
  mat phi(nbasis,nbasis);
  mat phi1(nbasis+1,nbasis+1);
  
  for (unsigned int i=0; i<nobs; i++) {
    xi = xobs(i);
    a = IntCos2(xi, kall, xmin, xrange, nbasis); // phijj
    phi.diag() = a;
    for (unsigned int j=1; j<=nbasis-1; j++) {
      for (unsigned int k=(j+1); k<=nbasis; k++) {
        b = IntCosCrossProd(xi,k,j,xmin,xrange); // phijk
        phi(j-1,k-1) = b;
        phi(k-1,j-1) = b;
      }  
    }
    c = IntConst2(xi,xmin,xrange); // phi00
    d = IntCos(xi,kall,xmin,xrange,nbasis); // phi0k
    phi1(0,0)=c;
    phi1(span(0), span(1,nbasis))=d.t();
    phi1(span(1,nbasis), span(0))=d;
    phi1(span(1,nbasis), span(1,nbasis))=phi;
    phixobs.col(i) = vech(phi1);
  }
  
  return phixobs;
}


// [[Rcpp::export]]
void update_lambda(arma::mat& lambda, arma::mat& lambdai, arma::mat resid_mat, arma::mat lambda_g0i, 
              double lambda_fn, unsigned int nparw) {
  
  double sse = 0;
  mat lambda_gni;
  mat lambda_gn;
  
  if (nparw > 1) {
    lambda_gni = lambda_g0i + resid_mat.t()*resid_mat;
    lambda_gn  = inv(lambda_gni);
    lambdai    = wishrnd(lambda_gn, lambda_fn); // random draw from a Wishart distribution
    lambda     = inv(lambdai);
  } else {
    sse       = accu(square(resid_mat));
    lambda_gn = lambda_g0i + sse;
    lambda    = lambda_gn + (2 * randg(distr_param(lambda_fn / 2, 1.0)));
    lambdai   = inv(lambda);
  }
  
  return;
}


// [[Rcpp::export]]
arma::mat update_phi(arma::mat lambdai, arma::mat ztz, arma::mat phi_v0i, arma::mat zdata, arma::mat betall,
                     unsigned int nparw, unsigned int nparz, unsigned int phidim) {
  mat phi        = zeros<mat>(nparz,nparw);
  mat phi_vec    = zeros<mat>(nparz,nparw);
  mat phi_vbni   = zeros<mat>(nparz*nparw,nparz*nparw);
  mat phi_vbni12 = zeros<mat>(nparz*nparw,nparz*nparw);
  mat ztb        = zeros<mat>(nparz,nparz);
  colvec a       = zeros<colvec>(nparz);
  colvec phi_bn  = zeros<colvec>(nparz*nparw);
  mat rMat;
  
  phi_vbni   = kron(lambdai, ztz) + phi_v0i; // nparz*nparw x nparz*nparw
  phi_vbni12 = chol(phi_vbni);
  ztb        = zdata.t() * betall;
  for (unsigned int j = 0; j < nparw; j++) {
    if (nparw > 1) {
      a = ztb * lambdai.col(j);
      phi_bn.subvec(nparz*j, nparz*(j+1)-1) = a;
    } else {
      phi_bn = ztb * lambdai(0, 0);
    }
    
  }
  rMat    = randn(phidim, 1);
  phi_vec = solve(phi_vbni, phi_bn + phi_vbni12.t()*rMat);
  if (nparw > 1) {
    for (unsigned int i = 0; i < nparw; i++) {
      phi.col(i) = phi_vec.rows(nparz*i, nparz*(i+1)-1);
    }  
  } else {
    phi = phi_vec;
  }
  
  
  return phi;
}


// w=0.5, m0.001, alpha=4.0, ndim=1
// [[Rcpp::export]]
arma::colvec GetMetVec(double w, double m0, double alpha, unsigned int ndim) {
  double beta    = m0*(alpha-1);
  int pmet       = 0;
  int icount     = 0;
  colvec beta_AM = beta*ones<colvec>(ndim);
  colvec m_AM    = m0*ones<colvec>(ndim);
  colvec met_var = m0*ones<colvec>(ndim);
  colvec Met     = ones<colvec>(7+3*ndim);
  
  Met(0) = ndim;
  Met(1) = pmet;
  Met(2) = icount;
  Met(3) = alpha;
  Met(4) = beta;
  Met(5) = w;
  Met(6) = m0;
  Met.subvec(7,          7 +   ndim-1) = beta_AM;
  Met.subvec(7 +   ndim, 7 + 2*ndim-1) = m_AM;
  Met.subvec(7 + 2*ndim, 7 + 3*ndim-1) = met_var;
    
  return Met;
}


// w=0.5, m0.001, alpha=4.0, ndim=1
// [[Rcpp::export]]
void SetMetVec(arma::colvec &Met, double w, double m0, double alpha, unsigned int ndim) {
  double beta    = m0*(alpha-1);
  int pmet       = 0;
  int icount     = 0;
  colvec beta_AM = beta*ones<colvec>(ndim);
  colvec m_AM    = m0*ones<colvec>(ndim);
  colvec met_var = m0*ones<colvec>(ndim);
  
  Met(0) = ndim;
  Met(1) = pmet;
  Met(2) = icount;
  Met(3) = alpha;
  Met(4) = beta;
  Met(5) = w;
  Met(6) = m0;
  Met.subvec(7,          7 +   ndim-1) = beta_AM;
  Met.subvec(7 +   ndim, 7 + 2*ndim-1) = m_AM;
  Met.subvec(7 + 2*ndim, 7 + 3*ndim-1) = met_var;
}



// [[Rcpp::export]]
arma::colvec AdaptMetVar(arma::colvec x) {
  unsigned int ndim = x(0);
  colvec beta_AM    = x.subvec(7, 7+ndim-1);
  double alpha      = x(3);
  colvec met_var    = beta_AM/randg(ndim, distr_param(alpha, 1.0));
  return met_var;
}


// [[Rcpp::export]]
void AdaptMetUpdate(arma::colvec &x) {
  unsigned int ndim = x(0);
  int icount        = x(2);
  colvec m_AM       = x.subvec(7 +   ndim, 7 + 2*ndim-1);
  colvec met_var    = x.subvec(7 + 2*ndim, 7 + 3*ndim-1);
  double w          = x(5);
  double alpha      = x(3);
  colvec mnew;
  
  x(2)                        = icount+1;
  x.subvec(7+ndim,7+2*ndim-1) = m_AM + (met_var - m_AM)/x(2); // met_var
  mnew                        = w*met_var + (1-w)*x.subvec(7+ndim,7+2*ndim-1);
  x.subvec(7,7+ndim-1)        = (alpha-1)*mnew; // beta_AM
}


// [[Rcpp::export]]
arma::colvec AdaptMetUpdate2(arma::colvec x) {
  colvec Met = x;
  unsigned int ndim = Met(0);
  int icount        = Met(2);
  colvec m_AM       = Met.subvec(7 +   ndim, 7 + 2*ndim-1);
  colvec met_var    = Met.subvec(7 + 2*ndim, 7 + 3*ndim-1);
  double w          = Met(5);
  double alpha      = Met(3);
  colvec mnew;
  
  Met(2)                        = icount+1;
  Met.subvec(7+ndim,7+2*ndim-1) = m_AM + (met_var - m_AM)/Met(2); // met_var
  mnew                          = w*met_var + (1-w)*Met.subvec(7+ndim,7+2*ndim-1);
  Met.subvec(7,7+ndim-1)        = (alpha-1)*mnew; // beta_AM
  return Met;
}


// fct=10
// [[Rcpp::export]]
void UpdateMet(arma::colvec& x, unsigned int nb, double fct) {
  double       pmet = x(1);
  unsigned int ndim = x(0);
  double       m0   = x(6);
  
  pmet /= nb;
  if (pmet > 0.6) { // Too few acceptances, decrease metm 
    m0 *= fct;
    SetMetVec(x, 0.5, m0, 4.0, ndim);
  } else if (pmet < 0.3) {  // Too few acceptances, decrease metm 
    m0 /= fct;
    if (m0 < 1e-8) m0 = 1e-8;
    SetMetVec(x, 0.5, m0, 4.0, ndim);
  } else {
    x(1) = 0;
  }
}


// fct=10
// [[Rcpp::export]]
arma::colvec UpdateMet2(arma::colvec x, unsigned int nb, double fct) {
  double       pmet = x(1);
  unsigned int ndim = x(0);
  double       m0   = x(6);
  colvec Met = x;

  pmet /= nb;
  if (pmet > 0.6) { // Too few acceptances, decrease metm 
    m0 *= fct;
    SetMetVec(Met, 0.5, m0, 4.0, ndim);
  } else if (pmet < 0.3) {  // Too few acceptances, decrease metm 
    m0 /= fct;
    if (m0 < 1e-8) m0 = 1e-8;
    SetMetVec(Met, 0.5, m0, 4.0, ndim);
  } else {
    Met(1) = 0;
  }
  
  return Met;
}





  