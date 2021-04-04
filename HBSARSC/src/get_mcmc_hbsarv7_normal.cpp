// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include "misc.h"

#define use_progress

#ifdef use_progress
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#endif


using namespace arma;

// [[Rcpp::export]]
Rcpp::List get_mcmc_hbsarv7(Rcpp::List& data_all,
                            arma::vec& id_group,
                            double xdelta,
                            arma::colvec& xgrid,
                            double xrange,
                            arma::uvec& xinx,
                            arma::vec& xover,
                            arma::vec& idn0,
                            arma::vec& nobs1,
                            Rcpp::List& iptHB1,
                            double xmin,
                            arma::mat& ztz,
                            arma::cube& vtv,
                            arma::cube& wtw,
                            arma::cube& phi2,
                            
                            arma::vec&  cpara,
                            Rcpp::List& probit_para,
                            Rcpp::List& gamma_para,
                            Rcpp::List& tau_para,
                            Rcpp::List& eta_para,
                            Rcpp::List& met_theta_para,
                            Rcpp::List& met_gamma0_para,                            
                            Rcpp::List& met_psi_para,
                            Rcpp::List& theta0_para,
                            arma::vec&  squish_para1,
                            Rcpp::List& squish_para2,
                            Rcpp::List& sigma_para,
                            Rcpp::List& v_para,
                            Rcpp::List& z_para,
                            Rcpp::List& phi_para,
                            Rcpp::List& fx_para,
                            Rcpp::List& basis_para,
                            
                            arma::colvec& alpha,
                            arma::mat& betall,
                            arma::mat& phi,
                            arma::mat& lambda,
                            arma::mat& lambdai,
                            arma::mat& thetall,
                            arma::colvec& theta0,
                            arma::colvec& xiparall,
                            double xipar0,
                            arma::colvec& gammall,
                            arma::colvec& tau2all,
                            arma::colvec& tauall,
                            arma::colvec& sigma2,
                            arma::colvec& sigma,
                            
                            bool iflagpsi,
                            bool iflagCenter,
                            bool iflaglm,
                            bool bFlagZ,
                            int  iflagpn,
                            int  iflagsc,
                            bool disp) {
  
  unsigned int ntot      = cpara(0);
  unsigned int ngroup    = cpara(1);
  unsigned int nparv     = cpara(2);
  unsigned int nparw     = cpara(3);
  unsigned int nparw2    = cpara(4);
  unsigned int nparz     = cpara(5);
  unsigned int phidim    = cpara(6);
  unsigned int nbasis    = cpara(7);
  unsigned int nblow     = cpara(8);
  unsigned int nblow0    = cpara(9);
  unsigned int maxmodmet = cpara(10);
  unsigned int nskip     = cpara(11);
  unsigned int nmcmcall  = cpara(12);
  unsigned int smcmc     = cpara(13);
  unsigned int nint      = cpara(14);
  unsigned int ntheta    = cpara(15);
  unsigned int nparx     = cpara(16);
  unsigned int hfun_nsim = cpara(17);
  int ingroup = ngroup;
  
  colvec ydata  = data_all["ydata"];
  colvec valpha = data_all["valpha"];
  colvec wbeta  = data_all["wbeta"];
  mat vdata     = data_all["vdata"];
  mat wdata     = data_all["wdata"];
  mat zdata     = data_all["zdata"];
  
  uvec id_probit1 = probit_para["id_probit1"];
  uvec id_probit0 = probit_para["id_probit0"];
  
  vec metw        = met_theta_para["metw"];
  vec metm        = met_theta_para["metm"];
  vec met_alpha   = met_theta_para["met_alpha"];
  vec mets        = met_theta_para["mets"];
  vec metv        = met_theta_para["metv"];
  vec met_mAM     = met_theta_para["met_mAM"];
  vec met_vAM     = met_theta_para["met_vAM"];
  vec met_beta    = met_theta_para["met_beta"];
  vec met_beta_AM = met_theta_para["met_beta_AM"];
  vec met_var_all = met_theta_para["met_var_all"];
  bool iflagAM    = met_theta_para["iflagAM"];
  vec pmet        = met_theta_para["pmet"];
  vec icount_AM   = met_theta_para["icount_AM"];
  
  double met_alpha_g0 = met_gamma0_para["met_alpha_g0"];
  double met_beta_g0  = met_gamma0_para["met_beta_g0"];
  
  //double met_alpha_zeta = met_zeta_para["met_alpha_zeta"];
  //double met_beta_zeta  = met_zeta_para["met_beta_zeta"];
  
  double met_alpha_psi = met_psi_para["met_alpha_psi"];
  double met_beta_psi  = met_psi_para["met_beta_psi"];
  
  double gmax            = gamma_para["gmax"];
  double gamma_mu        = gamma_para["gamma_mu"];
  double gamma_sigma     = gamma_para["gamma_sigma"];
  double gamma_beta      = gamma_para["gamma_beta"];
  double gamma_alpha     = gamma_para["gamma_alpha"];
  double gamma_prob      = gamma_para["gamma_prob"];
  double gamma_mu_m0     = gamma_para["gamma_mu_m0"];
  double gamma_mu_v02    = gamma_para["gamma_mu_v02"];
  double mu_bot          = gamma_para["mu_bot"];
  double mu_top          = gamma_para["mu_top"];
  double gamma_sigma_m0  = gamma_para["gamma_sigma_m0"];
  double gamma_sigma_v02 = gamma_para["gamma_sigma_v02"];
  double sigma_bot       = gamma_para["sigma_bot"];
  double sigma_top       = gamma_para["sigma_top"];
  colvec gam0vec         = gamma_para["gam0vec"];
  double wk              = gamma_para["wk"];
  
  double tau2_s0 = tau_para["tau2_s0"];
  double tau2_rn = tau_para["tau2_rn"];
  
  double eta02_s0 = eta_para["eta02_s0"];
  double eta02_rn = eta_para["eta02_rn"];
  double eta02    = eta_para["eta02"];
  double eta0     = eta_para["eta0"];
  colvec th0v     = eta_para["th0v"];
    
  double theta0_v02 = theta0_para["theta0_v02"];
  
  
  double zeta_mu        = squish_para1[0];
  double zeta_sigma2    = squish_para1[1];
  double zeta_sigma     = squish_para1[2];
  double zeta_m0        = squish_para1[3];
  double zeta_v02       = squish_para1[4];
  double zeta_s0        = squish_para1[5];
  double zeta_rn        = squish_para1[6];
  double zeta0          = squish_para1[7];  
  double omega0         = squish_para1[8];
  double psi_fixed      = squish_para1[9];
  double psi_log_fixed  = squish_para1[10];
  double psi0           = squish_para1[11];  
  double psi_log_mu     = squish_para1[12];
  double psi_log_sigma2 = squish_para1[13];
  double psi_log_sigma  = squish_para1[14];
  double psi_log_m0     = squish_para1[15];
  double psi_log_v02    = squish_para1[16];
  double psi_log_s0     = squish_para1[17];
  double psi_log_rn     = squish_para1[18];
  
  colvec zetall         = squish_para2["zetall"];
  colvec ezall          = squish_para2["ezall"];
  colvec psiall         = squish_para2["psiall"];
  colvec psiall_log     = squish_para2["psiall_log"];
  colvec omegall        = squish_para2["omegall"];
  mat hfunall           = squish_para2["hfunall"];
  
  double sigma2_alpha     = sigma_para["sigma2_alpha"];
  double sigma2_beta      = sigma_para["sigma2_beta"];
  double sigma2_alpha_m0  = sigma_para["sigma2_alpha_m0"];
  double sigma2_alpha_v02 = sigma_para["sigma2_alpha_v02"];
  double sigma2_beta_s0   = sigma_para["sigma2_beta_s0"];
  double sigma2_beta_rn   = sigma_para["sigma2_beta_rn"];
  
  mat alpha_v0i      = v_para["alpha_v0i"];
  colvec alpha_v0im0 = v_para["alpha_v0im0"];
  
  mat phi_v0i      = z_para["phi_v0i"];
  double lambda_fn = z_para["lambda_fn"];
  mat lambda_g0i   = z_para["lambda_g0i"];
  
  colvec kall  = basis_para["kall"];
  vec ktheta   = basis_para["ktheta"];
  
  mat phixgrid     = phi_para["phixgrid"];
  mat phi2xgrid    = phi_para["phi2xgrid"];
  mat phixdata     = phi_para["phixdata"];
  
  mat fxgridall   = fx_para["fxgridall"];
  colvec f0xgrid  = fx_para["f0xgrid"];
  colvec f0Bxgrid = fx_para["f0Bxgrid"];
  colvec fxobs    = fx_para["fxobs"];
  colvec fxobsm   = fx_para["fxobsm"];
  
  
  // Matrics for saving MCMC iterations
  unsigned int nparv_save = nparv > 0 ? nparv : 1;
  mat alphag = zeros<mat>(smcmc, nparv_save);
  
  
  mat sigmag        = zeros<mat>(smcmc,ngroup);    // Error variance
  vec sigma2_alphag = zeros<vec>(smcmc);
  vec sigma2_betag  = zeros<vec>(smcmc);
  mat tauallg       = zeros<mat>(smcmc,ntheta);    // StdDev for spectral coefficients
  vec eta0g         = zeros<vec>(smcmc);           // STDEV for upper level spectral coefficients
  vec gamma_mugg    = zeros<vec>(smcmc);           // Smoothing parameter for upper level model
  mat gammallg      = zeros<mat>(smcmc,ngroup);    // Smoothing parameters for lower level model
  mat thetam        = zeros<mat>(ntheta,ngroup);  
  mat thetas        = zeros<mat>(ntheta,ngroup);
  mat theta0g       = zeros<mat>(smcmc,ntheta);
  mat phig          = zeros<mat>(smcmc,phidim);
  mat lambdag       = zeros<mat>(smcmc,nparw2);
  mat betam         = zeros<mat>(size(betall));
  mat betas         = zeros<mat>(size(betall));
  
  vec gamma_alphag  = zeros<vec>(smcmc);
  vec gamma_betag   = zeros<vec>(smcmc);
  
  // squish
  mat zetag      = zeros<mat>(smcmc, ngroup);
  mat omegag     = zeros<mat>(smcmc, ngroup);
  colvec zeta0g  = zeros<colvec>(smcmc);
  colvec omega0g = zeros<colvec>(smcmc);
  mat psiallg    = zeros<mat>(smcmc, ngroup);
  colvec psi0g   = zeros<colvec>(smcmc);
  
  mat fxgridm       = zeros<mat>(nint+1,ngroup);
  mat fxgrids       = zeros<mat>(nint+1,ngroup);
  colvec f0xgridm   = zeros<colvec>(nint+1);
  colvec f0xgrids   = zeros<colvec>(nint+1);
  colvec f0Bxgridm  = zeros<colvec>(nint+1);
  colvec f0Bxgrids  = zeros<colvec>(nint+1);
  colvec fxobss     = zeros<colvec>(ntot);
  
  
  // local
  colvec yresid = zeros<colvec>(ntot);
  colvec yresidj;
  colvec rj;
  mat rMat;
  mat ranMat;
  colvec yj;
  mat vbni; 
  mat vbni12;
  colvec bn;
  colvec bj;
  double sse;
  
  // alpha
  mat vtvj;
  colvec vj;
  
  // beta
  mat wj;
  mat zbreg;
  colvec zbregj;
  colvec v0ib0;
  mat wtwj;
  colvec wjbj;
  
  // sigma
  double s2;
  double sigma2_rn;
  double sigma2_sn;
  
  // sigma_alpha, sigma_beta 
  double sigma2_beta_sn;
  double testpr;
  double var_met_sigma;
  double std_met_sigma;
  double sigma2_alpha_new;
  double lsum;
  double pnew;
  double pold;
  
  // phi
  mat phi_vbni   = zeros<mat> (nparz*nparw,nparz*nparw);
  mat phi_vbni12 = zeros<mat> (nparz*nparw,nparz*nparw);
  mat ztb        = zeros<mat> (nparz,nparz);
  colvec phi_vec = zeros<colvec> (nparz*nparw);
  colvec phi_bn  = zeros<colvec> (nparz*nparw);
  colvec a       = zeros<colvec> (nparz);
  mat resid_mat  = zeros<mat> (ngroup,nparw);
  mat lambda_gni = zeros<mat> (nparw, nparw);
  mat lambda_gn  = zeros<mat> (nparw, nparw);
  
  // theta
  mat vbn;
  mat vbn12;
  double sse_new;
  double sse_old;
  mat amat (ntheta,ntheta);
  mat amati (ntheta,ntheta);
  double g;
  mat phixj;
  mat phi2j;  
  colvec thetaj; 
  double gammaj; 
  double gpar;
  colvec gamjvec;
  colvec thv;
  colvec thvi;
  mat bmat;
  uvec zc;
  uvec zvi;
  unsigned int nz;
  double t0resid2_new;
  colvec tkresid2_new = zeros<colvec>(ntheta-1);
  colvec tresid2_new  = zeros<colvec>(ntheta);
  double t0resid2_old;
  colvec tkresid2_old = zeros<colvec>(ntheta-1);
  colvec tresid2_old  = zeros<colvec>(ntheta);
  double testp;
  Rcpp::List fpars;
  colvec resid_new;
  colvec resid_old;
  colvec fxgrid_new;
  colvec thvk;
  colvec fxj;
  colvec theta0k;
  colvec fxobs_old;
  colvec theta_old;
  double xiparj_old;
  colvec thetak_old;
  colvec gamvec;
  colvec ths;
  uvec zv;
  double met_var_all_new;
  //double met_beta_new;
  //double met_var0;
  double met_var;
  //double met_std0;
  double met_std;
  double xiparj_new;
  colvec thetak_new = zeros<colvec>(ntheta-1);
  colvec theta_new  = zeros<colvec>(ntheta);
  colvec hfunj_new;
  
  // squish
  double zeta_vn2;
  double zeta_vn;
  double zeta_mn;
  colvec resid;
  //double psij_old;
  double psij_log_old;
  double var_met_psi=0.0;
  double std_met_psi;
  double psij_log_new=0.0;
  double psij_new=0.0;
  double zetaj_old;
  //double omegaj_old;
  //double var_met_zetal;
  //double var_met_zeta;
  double std_met_zeta;
  double zetaj_new=0.0;
  double ezeta_new;
  double omegaj_new=0.0;
  colvec hfunm;
  double psi_log_v2;
  double psi_log_v; 
  double psi_log_mn;
  colvec zeta_sim;
  colvec ezim;
  colvec omega_sim;
  colvec psi_log_sim;
  colvec psi_sim;
  colvec hfj;
  double zeta_sn;
  double psi_log_sn;
  
  
  colvec fpara;
  unsigned int kbot;
  colvec vpar;
  colvec thvb;
  Rcpp::List f0pars;
  
  // theta0
  colvec t0var;
  double xi_vn;
  double xi_bn;
  
  colvec ymall;
  unsigned int tk;
  
  colvec resid2;
  colvec gamk;
  double tau2_sn;
  unsigned int tbot;
  
  
  // tau2
  colvec resid_tau  = zeros<colvec> (ngroup);
  colvec resid2_tau = zeros<colvec> (ngroup);
  colvec tau2spec;
  
  // eta02
  colvec resid2_eta;
  double eta02_sn;
  double etv0;
  
  // gamma
  colvec u1;
  double u2;
  double bmin;
  double bmax;
  double gamma_sum;
  double lgamma_sum;
  double var_met_g0;
  double std_met_g0;
  double gamma_mu_new;
  double gamma_sigma_new;
  double gamma_sigma2_new;
  double gamma_beta_new;
  double gamma_alpha_new;
  double gamma_prob_new;
  double testpa;
  colvec gam0vec_new;
  colvec th0v_new;
  colvec theta02;
  //double gamma0;
  
  uvec z;
  uvec zi;
  
  colvec fxobs_new;
  
  // adpative metroplis
  vec pmetm = zeros<colvec> (nbasis);
  unsigned int pok=0;
  
#ifdef use_progress
  Progress p(nmcmcall, true);
#endif
  
  unsigned int isave = 0;
  for (unsigned int imcmc=0; imcmc<nmcmcall; imcmc++) {
    
#ifdef use_progress
    if ( Progress::check_abort() ) { // check for user abort
      break;
    }
    p.increment();
#else
    R_CheckUserInterrupt(); // void return type (not bool!!)
#endif
    
    // Generate fixed effect alpha
    if (nparv > 0) {
      yresid = ydata - wbeta - fxobs;   // residuals
      vbni   = alpha_v0i;
      bn     = alpha_v0im0;
      for (unsigned int j = 0; j < ngroup; j++) {
        uvec gidx = iptHB1[j];
        vtvj = vtv.slice(j);
        vj   = vdata.elem(gidx);
        rj   = yresid.elem(gidx);
        vbni += vtvj/sigma2(j);
        bn   += (vj.t() * rj )/sigma2(j);
      }
      
      vbni12 = chol(vbni);
      rMat   = randn(nparv, 1);
      alpha  = solve(vbni, bn + vbni12.t()*rMat);
      valpha = vdata * alpha;
    } // End generate alpha
    
    // Do HB parametric model
    yresid = ydata - fxobs - valpha;
    //--------------------------------------------
    // Generate Beta_j
    zbreg = zdata * phi;
    sse = 0;
    for (unsigned int j = 0; j < ngroup; j++) {
      uvec gidx   = iptHB1[j];
      zbregj = zbreg.row(j).t();
      v0ib0  = lambdai * zbregj;
      wj     = wdata.rows(gidx);
      yj     = yresid.elem(gidx);
      wtwj   = wtw.slice(j);
      vbni   = wtwj / sigma2(j) + lambdai;
      vbni12 = chol(vbni);   // Upper triangluar: t(vbni12) %*% vbni12 = vbni
      bn     = (wj.t() * yj) / sigma2(j) + v0ib0;
      
      // Generate draw by using solve: avoides inverse of vni
      ranMat        = randn(nparw, 1);
      bj            = solve(vbni, bn + vbni12.t()*ranMat);
      betall.row(j) = bj.t();
      // Compute residual
      wjbj             = wj * bj;
      rj               = yj - wjbj;
      wbeta.rows(gidx) = wjbj;
      
      //--------------------------------------------
      // Generate Sigma2, the error variance
      sigma2_rn = sigma2_alpha + nobs1(j);
      sigma2_sn = sigma2_beta  + accu(square(rj));
      s2        = iflaglm ? 1 : sigma2_sn / (2 * randg(distr_param(sigma2_rn / 2, 1.0))); // Probit likelihood: variance = 1
      sigma2(j) = s2;
      sigma(j)  = sqrt(s2);
    }
    
    
    //-----------------------------------------------------------------
    // Generate sigma2_alpha and sigma2_beta
    // HB model for variance: sigma_j^2 ~ IG(alpha/2,beta/2)
    // Hyperpriors:   alpha ~ N(m0,v02)I(alpha>2) and beta~G(r0/2,s0/2)
    // Generate sigma2_beta
    
    
    sigma2_beta_sn = sigma2_beta_s0 + accu(1/sigma2);
    sigma2_beta    = randg(distr_param(sigma2_beta_rn/2,2/sigma2_beta_sn));
    // Generate sigma2_alpha
    var_met_sigma    = met_beta_g0/randg(distr_param(met_alpha_g0, 1.0));
    std_met_sigma    = sqrt(var_met_sigma);
    sigma2_alpha_new = rcpp_rndtnb(sigma2_alpha,std_met_sigma,2); // Generate candidate from truncated normal
    
    lsum        = accu(log(sigma2));
    //  IG Likelihood
    testpr      = ingroup*log(sigma2_beta/2)*(sigma2_alpha_new - sigma2_alpha)/2;
    testpr      = testpr - ingroup*(lgamma(sigma2_alpha_new/2)-lgamma(sigma2_alpha/2));
    testpr      = testpr - lsum*(sigma2_alpha_new-sigma2_alpha)/2;
    // Prior
    testpr      = testpr - pow(sigma2_alpha_new - sigma2_alpha_m0, 2)/(2*sigma2_alpha_v02);
    testpr      = testpr + pow(sigma2_alpha     - sigma2_alpha_m0, 2)/(2*sigma2_alpha_v02);
    // Generating function
    pnew        = 1-normcdf(2.0,sigma2_alpha,std_met_sigma);
    pold        = 1-normcdf(2.0,sigma2_alpha_new,std_met_sigma);
    testpr      = testpr + log(pnew) - log(pold);
    if (log(randu())<testpr) sigma2_alpha = sigma2_alpha_new;
    
    
    //------------------------------------------
    // Generate Phi (Beta = Zdata*Phi + delta)
    phi_vbni   = kron(lambdai, ztz) + phi_v0i;
    phi_vbni12 = chol(phi_vbni);
    ztb        = zdata.t() * betall;
    for (unsigned int j = 0; j < nparw; j++) {
      if (nparw > 1) {
        a = ztb * lambdai.col(j);
      }
      else {
        a = ztb * lambdai(0, 0);
      }
      
      phi_bn.subvec(nparz*j, nparz*(j+1)-1) = a;
      
    }
    rMat = randn(phidim, 1);
    phi_vec  = solve(phi_vbni, phi_bn + phi_vbni12.t()*rMat);
    for (unsigned int i = 0; i < nparw; i++) {
      phi.col(i) = phi_vec.subvec(nparz*i, nparz*(i+1)-1);
    }
    
    //------------------------------------------
    // Generate Lambda from Beta = Zdata*breg + delta  
    resid_mat = betall - zdata * phi;
    if (nparw > 1) {
      lambda_gni = lambda_g0i + resid_mat.t()*resid_mat;
      lambda_gn  = inv(lambda_gni);
      lambdai    = wishrnd(lambda_gn, lambda_fn); // random draw from a Wishart distribution
      lambda     = inv(lambdai);
    }
    else {
      sse       = accu(square(resid_mat));
      lambda_gn = lambda_g0i + sse;
      lambda    = lambda_gn + (2 * randg(distr_param(lambda_fn / 2, 1.0)));
      lambdai   = 1/lambda;
    }
    
    //------------------------------------------
    // Do HB nonparametric model
    yresid = ydata - valpha - wbeta;   // Residuals
    
    // Free model: do not need numerical integration
    if(iflagsc==0) {
      // Loop over groups to generate theta_j
      for (unsigned int j = 0; j < ngroup; j++) { // Loop over groups to generate theta_j
        uvec gidx = iptHB1[j];
        rj      = yresid.elem(gidx);
        phixj   = phixdata.rows(gidx);
        phi2j   = phi2.slice(j); // Get phi_j'phi_j
        thetaj  = thetall.col(j);
        gammaj  = gammall(j);
        if (nparx > 0) {
          gamjvec = zeros<colvec>(kall.n_elem + nparx);
          gamjvec.head(nparx) = ones<colvec>(nparx);
          gamjvec.tail(kall.n_elem) = exp(-kall*gammaj);
        } else {
          gamjvec = exp(-kall*gammaj);
        }
        colvec thv = tau2all%gamjvec;
        uvec id = find(thv < 1e-10);
        if (id.n_elem > 0) {
          thv.elem( id ) = 1e-10*ones<colvec>(id.n_elem);
        }
        colvec thvi = 1/thv;
        mat bmat    = phi2j/sigma2(j);
        mat vbni    = bmat;
        vbni.diag() += thvi;
        // Issue with poorly conditioned vbni
        if (rcond(vbni) > 1e-10) {
          vbni12 = chol(vbni); // Upper triangluar: t(vbni12) %*% vbni12 = vbni
          bn     = (phixj.t() * rj)/sigma2(j) + theta0/thv;
          // Generate draw by using solve: avoides inverse of vni
          rMat   = randn(ntheta, 1);
          thetaj = solve(vbni, bn + vbni12.t()*rMat);
        } else {
          // "solve" failes due to poorly conditioned number
          // Try lemma https://math.stackexchange.com/questions/17776/inverse-of-the-sum-of-matrices
          // inv(A+B) = inv(A) - inv(A)*B*inv(A)/(1+g)
          // g = trace(B*inv(A)) if g not 1
          // Set A = thvi and B = phi2j/sigma2[j]
          amat  = zeros<mat>(ntheta,ntheta);
          amati = zeros<mat>(ntheta,ntheta);
          amat.diag()  = thvi;
          amati.diag() = thv;
          mat bamati   = bmat * amati;
          g            = accu(bamati.diag());
          vbn          = amati - amati * bmat * amati/(1+g);
          bn           = vbn * ((phixj.t() * rj)/sigma2(j) + theta0/thv);
          vbn12        = chol(vbn);
          thetaj       = bn + vbn12.t()*randn(ntheta, 1);
        }
        
        thetall.col(j) = thetaj;
        
        // Compute new fj
        fxgridall.col(j) = phixgrid * thetaj;  // fj computed on xgrid
        fxj              = phixj * thetaj;     // fj computed at observations
        fxobs.rows(gidx) = fxj;
      } // End loop over groups to generate theta_j
    } else {
      //theta00 = theta0(0);                   // Upper-level Fourier coefficient for constant
      theta0k = theta0.tail(theta0.n_elem-1);  // Upper-level Fourier coefficient for consines
      
      // Loop over groups to generate theta_j
      for (unsigned int j = 0; j < ngroup; j++) {
        
        uvec gidx  = iptHB1[j];
        yresidj    = yresid(gidx);                       // Parametric residuals for group j
        fxobs_old  = fxobs(gidx);                        // fJ(x) at current values of theta
        theta_old  = thetall.col(j);                     // Old value for theta_j
        xiparj_old = xiparall(j);                        // theta_j0 = exp(xi_j)
        thetak_old = theta_old.tail(theta_old.n_elem-1); //Skip theta_j0
        gammaj     = gammall(j);
        if (nparx > 0) {
          gamvec                   = zeros<colvec>(kall.n_elem+nparx);
          gamvec.head(nparx)       = ones<colvec>(nparx); // pad non-frequency variances with ones
          gamvec.tail(kall.n_elem) = exp(-gammaj*kall);
        } else {
          gamvec                   = exp(-gammaj*kall);
        }
        
        thv                        = tau2all % gamvec;
        ths                        = sqrt(thv);
        // Worry about zero variances: which creates testp=NaN.
        // If Var(theta_jk) is too small, set theta_jk = theta_0k
        zv      = find(ths < 1e-10);
        zvi     = find(ths >= 1e-10);
        nz      = zv.n_elem;
        if (nz > 0) {
          thv(zv) = 1e-20*ones<colvec>(nz);
          ths(zv) = 1e-10*ones<colvec>(nz);  
        }
        
        //-----------------------------------------------------------------
        // get variances for t-dist random walk
        met_var_all_new  = met_beta_AM(j)/randg(distr_param(met_alpha(j), 1.0));
        double ck        = metw(j)*met_var_all_new + (1-metw(j))*met_mAM(j);
        //met_beta_new     = (met_alpha(j)-1)*ck;
        
        ck        = 5.66;   // Constant from Harito
        //met_var0  = ck*met_var_all_new;
        met_var   = ck*met_var_all_new;
        //met_std0  = sqrt(met_var0);
        met_std   = sqrt(met_var);
        // Somtimes met_std goes to 0.  Shake it up
        if((met_std < 1e-5)&&(randu()<0.1)){
          met_std = 1e-5;
          met_var = pow(met_std, 2);
        }
        
        // Random walk from for normals
        xiparj_new                        = xiparj_old + met_std*ths(0)*randn(); 
        thetak_new                        = thetak_old + met_std*ths.tail(ths.n_elem-1)%randn(ntheta-1,1);
        theta_new(0)                      = exp(xiparj_new);
        theta_new.tail(thetak_new.n_elem) = thetak_new;
        // If variance is too small, replace theta_jk with theta_0k
        if (nz>0) theta_new(zv) = theta0(zv);
        t0resid2_new                          = pow(xiparj_new - xipar0,  2);
        tkresid2_new                          = pow(thetak_new - theta0k, 2);
        tresid2_new                           = zeros<colvec>(ntheta);
        tresid2_new(0)                        = t0resid2_new;
        tresid2_new.tail(tkresid2_new.n_elem) = tkresid2_new;
        
        t0resid2_old                          = pow(xiparj_old - xipar0,  2); // squared residual  k=0
        tkresid2_old                          = pow(thetak_old - theta0k, 2); // squared residuals k>0
        tresid2_old                           = zeros<colvec>(ntheta);
        tresid2_old(0)                        = t0resid2_old;
        tresid2_old.tail(tkresid2_old.n_elem) = tkresid2_old;
        
        // HB Prior
        // Remove entries where thv = 0 & thetak_new = theta_0k
        if (nz>0) {
          tresid2_new = tresid2_new(zvi);
          tresid2_old = tresid2_old(zvi);
          thv         = thv(zvi);
        }
        
        testp = 0;
        testp = testp - accu(tresid2_new/thv)/2 + accu(tresid2_old/thv)/2;
        
        //-------------------------------------------
        // Generate Squish function
        if (iflagsc == 4) {
          // HB model for log(psi)
          // Generate new psi_log
          if (iflagpsi) {
            //psij_old      = psiall(j);
            psij_log_old  = psiall_log(j);
            // Random walk metropolis on zeta and psi_log
            var_met_psi  = met_beta_psi/randg(distr_param(met_alpha_psi, 1.0));
            std_met_psi  = sqrt(var_met_psi);
            psij_log_new = psij_log_old + std_met_psi*randn();
            psij_new     = exp(psij_log_new);
            testp        = testp - (pow(psij_log_new-psi_log_mu, 2) - pow(psij_log_old-psi_log_mu, 2))/(2*psi_log_sigma2);
          } else {
            // Constant value for psi = psi_fixed
            psij_log_new   = psi_log_fixed;
            psij_new       = psi_fixed;
          }
          // End generate new psi_log
          //------------------------------------------------
          // Generate zeta and omega = a + (b-a)*exp(zeta)/(1+exp(zeta)) where a < x < b
          zetaj_old      = zetall(j);
          //omegaj_old     = omegall(j);
          //var_met_zeta   = met_beta_zeta/randg(distr_param(met_alpha_zeta, 1.0));
          std_met_zeta   = sqrt(var_met_psi);
          zetaj_new      = zetaj_old + std_met_zeta*randn();
          ezeta_new      = exp(zetaj_new);
          omegaj_new     = xmin + xrange*ezeta_new/(1+ezeta_new);
          hfunj_new      = rcpp_GetSquish(omegaj_new,psij_new,xgrid);
          
          // HB prior for zetaj ~ N(zeta0,zeta_v2)
          testp = testp - (pow(zetaj_new-zeta_mu,2) - pow(zetaj_old-zeta_mu,2))/(2*zeta_sigma2);
        }  // End Generate Squish Function parameters
        
        
        // Get fJ(x) for theta_new
        // Get function for group j
        fpars = Rcpp::List::create(Rcpp::Named("iflagsc")     = iflagsc,
                                   Rcpp::Named("iflagpn")     = iflagpn,
                                   Rcpp::Named("iflagCenter") = iflagCenter,
                                   Rcpp::Named("iflagZ")      = bFlagZ,
                                   Rcpp::Named("theta")       = theta_new,
                                   Rcpp::Named("phix")        = phixgrid,
                                   Rcpp::Named("delta")       = xdelta,
                                   Rcpp::Named("range")       = xrange,
                                   Rcpp::Named("xmin")        = xmin,
                                   Rcpp::Named("xgrid")       = xgrid,
                                   Rcpp::Named("hfun")        = hfunj_new);
        
        fxgrid_new = rcpp_Getfx(fpars);
        
        //---------------------------------------------------
        // Compute fj(x_obs)
        // fxobs = GetUpfxobs(xgrid,xinx,xover,fxgrid)
        uvec a = xinx(gidx);
        vec  b = xover(gidx);
        fxobs_new = rcpp_GetSCfxobs(xgrid, a, b, fxgrid_new, iflagpn);
        
        //---------------------------------------------------
        // Metropolis Test
        // Compute Likelihood
        resid_new = yresidj - fxobs_new;  // Note: fxobs_new_j only for group j
        sse_new = accu(square(resid_new));
        
        resid_old = yresidj - fxobs_old;
        sse_old = accu(square(resid_old));
        
        // Compute log test probabability
        // Likelihood
        testp = testp - (sse_new - sse_old)/(2*sigma2(j)); // Likelihood    
        
        
        if (log(randu()) < testp) { // Accept candidate
          xiparall(j)      = xiparj_new;
          thetall.col(j)   = theta_new;
          fxobs.rows(gidx) = fxobs_new;
          fxgridall.col(j) = fxgrid_new;
          if (iflagsc==2) {
            hfunall.col(j) = hfunj_new;
            omegall(j)     = omegaj_new;
            zetall(j)      = zetaj_new;
            psiall(j)      = psij_new;
            psiall_log(j)  = psij_log_new;
          }
          met_var_all(j)   = met_var_all_new;
          pmet(j)++;
        }
        
        if (iflagAM) { // Update mean, variance, alpha, and beta for metropolis
          icount_AM(j)  = icount_AM(j)+1;
          met_mAM(j)    = met_mAM(j)+(met_var_all(j)-met_mAM(j))/(icount_AM(j));
          met_vAM(j)    = ((icount_AM(j)-1)/icount_AM(j))*met_vAM(j)+(pow(met_var_all(j)-met_mAM(j), 2))/icount_AM(j);
        }
        
        ck              = metw(j)*met_var_all(j)+(1-metw(j))*met_mAM(j);
        met_beta_AM(j)  = (met_alpha(j)-1)*ck;
        
      } //  End loop over groups
    } // End Generate theta with shape constraints
    
    // ---------------------------------------------------------------
    //  Generate upper-level spectral coefficients
    // -----------------------------------------------
    //  Generate upper-level model theta_0k
    //  theta_jk ~ N(theta_0k,tau_k^2*eta_j^2*exp(-k*gamma_j)) for k>0
    //  theta_0k ~ N(0,eta0^2*(1+k/gamma_beta)^(-gamma_alpha))
    //  Note: ktheta = c(rep(0,nparx),1:nbasis)
    if (nparx > 0) {
      t0var = zeros<colvec>(ntheta);
      t0var.head(nparx)       = theta0_v02*ones<colvec>(nparx); // Pad front with v02  
      t0var.tail(th0v.n_elem) = th0v;
    } else {
      t0var = th0v;           // Prior variances of theta0  
    }
    // Generate theta_0k.  Loop over k
    for (unsigned int k = 0; k < ntheta; k++) {
      // Get variance of theta_jk
      // Note: ktheta = c(rep(0,nparx),1:nbasis)
      // theta_jk ~ N(theta_0k,tau_k^2) if 1 <= nparx
      // theta_jk ~ N(theta_0k,tau_k^2*eta_j^2*exp(-k*gamma_j) if k > nparx
      thvk = tau2all(k)*ones<colvec>(ngroup);
      if (k >= nparx) { 
        thvk = thvk%exp(-ktheta(k)*gammall); // Var(theta_jk) for j = 1, ..., J
      } 
      
      if ((iflagsc>0)&&(k==0)) {
        // Got intercept, which is positive
        xi_vn     = 1/(ingroup/tau2all(k) + 1/theta0_v02);
        xi_bn     = xi_vn*accu(xiparall)/tau2all(k);    // mean
        xipar0    = xi_bn + sqrt(xi_vn)*randn();
        theta0(0) = exp(xipar0);
      }else{
        xi_vn     = 1/(accu(1/thvk) + 1/t0var(k));  // variance
        xi_bn     = xi_vn*accu(thetall.row(k).t()/thvk);    // mean
        theta0(k) = xi_bn + sqrt(xi_vn)*randn();
      }
    } // End loop to generate theta_0k
    
    //--------------------------------------------------------
    // Generate upper-level model for squish function paramters
    if (iflagsc==4) {
      //---------------------------------------------------------
      // U-shaped has squish function 
      // h(x) = (1-exp(psi*(x-omega)))/(1+exp(psi*(x-omega))
      // HB model for omega_j
      // omegaj = xmin + xrange*exp(zetaj)
      // HB model for zeta
      // zetaj       ~ N(zeta_mu,zeta_sigma2)
      // zeta_mu     ~ N(m0,v2)
      // zeta_sigma2 ~ IG(r0/2,s0/2)
      //---------------------------------------------------------
      // Generate zeta_mu
      zeta_vn2 = 1/(ingroup/zeta_sigma2 + 1/zeta_v02);
      zeta_vn  = sqrt(zeta_vn2);
      zeta_mn  = zeta_vn2*(accu(zetall)/zeta_sigma2 + zeta_m0/zeta_v02);
      zeta_mu  = zeta_mn + zeta_vn*randn();
      zeta0    = zeta_mu;
      //---------------------------------------------------------
      // Generate zeta_sigma2
      resid       = zetall - zeta_mu;
      zeta_sn     = zeta_s0 + accu(square(resid));
      zeta_sigma2 = zeta_sn/(2*randg(distr_param(zeta_rn/2, 1.0)));
      
      hfunm = mean(hfunall, 1); // apply(hfunall,1,mean)  // initialize 
      //------------------------------------------------
      // Generate psi for upper-leve model
      if (iflagpsi>0) {
        // Generate psi_log_mu
        psi_log_v2  = 1/(ingroup/psi_log_sigma2 + 1/psi_log_v02);
        psi_log_v   = sqrt(psi_log_v2);
        psi_log_mn  = psi_log_v2*(accu(psiall_log)/psi_log_sigma2 + psi_log_m0/psi_log_v02);
        psi_log_mu  = psi_log_mn + psi_log_v*randn();
        //psi0_log    = psi_log_mu;
        // Generate psi_log_sigma2
        resid          = psiall_log - psi_log_mu;
        sse            = accu(square(resid));
        psi_log_sn     = psi_log_s0 + sse;
        psi_log_sigma2 = psi_log_sn/(2*randg(distr_param(psi_log_rn/2, 1.0)));
        psi_log_sigma  = sqrt(psi_log_sigma2);
        
        psi0        = exp(psi_log_mu);
        //psi0_log    = psi_log_mu;
      }
      
      //-----------------------------------------------------------
      // Need to approximate E[h(x|omega,psi)|Data] by simulation
      // Generate zeta for simulation
      zeta_sim  = zeta_mu + zeta_sigma*randn(hfun_nsim, 1);
      ezim      = exp(zeta_sim);
      omega_sim = xmin + xrange*ezim/(1+ezim);
      // Generate psi for simulation
      if (iflagpsi) {
        psi_log_sim = psi_log_mu + psi_log_sigma*randn(hfun_nsim, 1);
        psi_sim     = exp(psi_log_sim);
      }else{
        psi_sim     = psi_fixed*ones<colvec>(hfun_nsim);
      }
      hfunm = zeros<colvec>(nint+1);
      // Simulate E[h(x|omega,psi)] omgea and psi from HB model for zeta and log(psi)
      for (unsigned int jh=0; jh<hfun_nsim; jh++){
        // Get squish function at simulated omega and psi
        hfj    = rcpp_GetSquish(omega_sim(jh),psi_sim(jh),xgrid);
        hfunm  = hfunm + hfj;
      }
      hfunm   = hfunm/hfun_nsim;
    } // End generate upper-level model for squish function parameters
    
    
    //--------------------------------------------------------------------------
    // Compute upper-level f0
    if (iflagsc==0) {
      f0xgrid  = phixgrid * theta0;
      f0Bxgrid = f0xgrid;
    } else {  // Shape constraints: Need to do bias correction
      // Note: some issues with theta_j0 = exp(xipar_j) is log normal
      fpara = theta0;    // Upper level spectral coefficients
      //fpara[1]         = fpara[1]*exp(tau2all[1]/2)  # log normal theta_{0,j}
      
      // Integrate exp(-k*gamma_j) where gamma_j ~ G(alpha,beta)
      kbot = nparx;
      // Note: E[eta_j^2] = 1 and E[exp(-gamma_j*k)] = (1+k/beta)^(-alpha)
      
      if (nparx>0) {
        vpar = zeros<colvec>(ntheta);
        vpar.head(nparx) = ones<colvec>(nparx);    // pad front of variance with ones 
        vpar.tail(gam0vec.n_elem) = gam0vec;
      } else {
        vpar = gam0vec;  
      }
      thvb    = tau2all%vpar;
      etv0    = exp(tau2all(0));
      thvb(0) = exp(2*xipar0  + tau2all(0))*(etv0-1);  // lognormal variance V(theta_j0)
      
      if(iflagsc < 4) hfunm = zeros<colvec>(nint+1);
      //thvb   = apply(thetall,1,sd)^2
      f0pars = Rcpp::List::create(Rcpp::Named("iflagsc")     = iflagsc,
                                  Rcpp::Named("iflagpn")     = iflagpn,
                                  Rcpp::Named("iflagCenter") = iflagCenter,
                                  Rcpp::Named("iflagZ")      = bFlagZ,
                                  Rcpp::Named("theta")       = fpara,
                                  Rcpp::Named("vpara")       = thvb,
                                  Rcpp::Named("phix")        = phixgrid,
                                  Rcpp::Named("phix2")       = phi2xgrid,
                                  Rcpp::Named("range")       = xrange,
                                  Rcpp::Named("delta")       = xdelta,
                                  Rcpp::Named("xmin")        = xmin,
                                  Rcpp::Named("xgrid")       = xgrid,
                                  Rcpp::Named("hfun")        = hfunm);
      
      f0xgrid  = rcpp_Getf0x(f0pars);
      
      // Get biased f0
      fpars = Rcpp::List::create(Rcpp::Named("iflagsc")     = iflagsc,
                                 Rcpp::Named("iflagpn")     = iflagpn,
                                 Rcpp::Named("iflagCenter") = iflagCenter,
                                 Rcpp::Named("iflagZ")      = bFlagZ,
                                 Rcpp::Named("theta")       = theta0,
                                 Rcpp::Named("phix")        = phixgrid,
                                 Rcpp::Named("delta")       = xdelta,
                                 Rcpp::Named("range")       = xrange,
                                 Rcpp::Named("xmin")        = xmin,
                                 Rcpp::Named("xgrid")       = xgrid,
                                 Rcpp::Named("hfun")        = hfunm);
      
      f0Bxgrid = rcpp_Getfx(fpars);  // f using theta0, which is biased
    } // End compute f0
    
    
    //-----------------------------------------------------------------
    // End nonparametric model
    //-----------------------------------------------------
    ymall = valpha + wbeta + fxobs;  // mean of Y
    
    
    // Generate tau2all for HB Smoothing Prior
    tk = 0;   // Index to handle theta_j0 = exp(xi_j) with shape constraints
    // If shape constraint, theta_j0 = exp(xipar_j) & xipar_j ~ N(xipar_0,tau_0^2)
    if (iflagsc>0) {
      resid_tau  = xiparall - xipar0;
      resid2_tau = square(resid_tau);
      tau2_sn    = tau2_s0 + accu(resid2_tau);
      tau2all(0) = tau2_sn/(2*randg(distr_param(tau2_rn/2, 1.0)));
      tk         = 1;   //  Increment index to skip tau2_0
    }
    for (unsigned int k=tk; k<ntheta; k++) {
      resid_tau  = thetall.row(k).t() - theta0(k);
      resid2_tau = square(resid_tau);
      gamk       = ones<colvec>(resid2_tau.n_elem);
      if (k >= nparx) gamk = exp(-ktheta(k)*gammall);
      tau2_sn    = tau2_s0 + accu(resid2_tau/gamk);
      tau2all(k) = tau2_sn/(2*randg(distr_param(tau2_rn/2, 1.0)));
    }
    // Get std dev
    tauall      = sqrt(tau2all);
    
    //--------------------------------------------------------------------
    // Generate eta0^2 where theta_0k ~ N(0,eta0^2*(1+kall/gamma_beta)^(-gamma_alpha))
    kbot       = nparx;
    resid2_eta = square(theta0.subvec(kbot,ntheta-1));
    sse        = accu(resid2_eta/gam0vec);
    eta02_sn   = eta02_s0 + sse;
    eta02      = eta02_sn/(2*randg(distr_param(eta02_rn/2, 1.0)));
    eta0       = sqrt(eta02);
    th0v       = eta02*gam0vec;
    
    
    //--------------------------------------------------------------------
    // Generate gamma for each population
    tbot     = nparx;
    tau2spec = tau2all.subvec(tbot,ntheta-1);  // Variance parameters
    for (unsigned int j=0; j<ngroup; j++) {  // Loop over groups to generate gammaj with slice sampling
      gammaj     = gammall(j);
      gamvec     = exp(-kall*gammaj);
      resid      = thetall(span(tbot,ntheta-1),j) - theta0.subvec(tbot,ntheta-1);
      resid2     = square(resid);
      colvec ck  = resid2/(2*tau2spec);
      
      // Worry about ck = 0
      z          = find(ck == 0);   // Find index where ck == 0
      zi         = find(ck != 0);   // Find index where ck == 0
      nz         = z.n_elem;
      if (nz == nbasis) {  // all theta's are zeros!
        // Reset gamma_j 
        gammaj     = 1;
        gammall(j) = gammaj;
      }else{
        if(nz>0) ck(z) = ones<colvec>(nz);   // Set zeros to 1 for the time being
        u1   = randu(nbasis, 1);
        bn   = gammaj + (log(ck - log(u1)%gamvec) - log(ck))/kall;
        if (nz>0) bn = bn(zi);   // drop the z's that are 0
        
        // Slice sampling for gamma^(alpha-1)
        u2   = randu();
        bmin = gammaj*(pow(u2, 1/(gamma_alpha - 1)));
        colvec bg = zeros<colvec>(bn.n_elem+1);
        bg.head(bn.n_elem) = bn;
        bg.tail(1) = gmax;
        bmax = min(bg);  // Sometimes gamma wanders off to a large number.  Limit it to 5
        if(bmin>bmax){
          gammaj  = 1.0;
        }else{
          u2   = randu();
          gpar = wk - gamma_beta;
          gammaj = bmax + log(u2 + (1.0-u2)*exp((bmin-bmax)*gpar))/gpar;
        }  
        gammall(j) = gammaj;
        
      } // End if ck all zeros
    } // End generate gammaj for each poulation
    // Summary stats used in generating hyperparameters
    gamma_sum  = accu(gammall);
    lgamma_sum = accu(log(gammall));
    
    //------------------------------------------------------------------------
    // Generate alpha and beta for gamma_j ~ G(alpha,beta)
    // Generate gamma_mu = gamma_alpha/gamma_beta and gamma_sigma2 = gamma_mu/gamma_beta
    // from truncated normal distribution
    kbot             = nparx;
    var_met_g0       = met_beta_g0/randg(distr_param(met_alpha_g0, 1.0)); // Variance for Metropolis
    std_met_g0       = sqrt(var_met_g0);                                  // STD DEV for Metropolis
    
    gamma_mu_new     = rcpp_rndtnab(gamma_mu,std_met_g0,mu_bot,mu_top);
    gamma_sigma_new  = rcpp_rndtnab(gamma_sigma,std_met_g0,sigma_bot,sigma_top);
    gamma_sigma2_new = pow(gamma_sigma_new,2);
    gamma_beta_new   = gamma_mu_new/gamma_sigma2_new;
    gamma_alpha_new  = gamma_beta_new*gamma_mu_new;
    gamma_prob_new   = R::pgamma(gmax,gamma_alpha_new,1/gamma_beta_new, 1, 0);
    
    // Likelihood for gamma_j ~ G(alpha,beta)I(gamma_j < gmax)
    testpa = 0;
    testpa = testpa - ingroup*(log(gamma_prob_new) - log(gamma_prob));  // Normalizing constant
    testpa = testpa + ingroup*gamma_alpha_new*log(gamma_beta_new);
    testpa = testpa - ingroup*gamma_alpha*log(gamma_beta);
    testpa = testpa - ingroup*lgamma(gamma_alpha_new);
    testpa = testpa + ingroup*lgamma(gamma_alpha);
    testpa = testpa + (gamma_alpha_new - gamma_alpha)*lgamma_sum;
    testpa = testpa - (gamma_beta_new  - gamma_beta)*gamma_sum;
    // Likelihood theta_0j ~ N(0,eta0^2*(1+k/gamma_beta)^(-gamma_alpha))
    gam0vec_new    = pow(1+kall/gamma_beta_new, -gamma_alpha_new);
    th0v_new       = eta02*gam0vec_new;
    theta02        = square(theta0.subvec(kbot, ntheta-1));
    testpa         = testpa - accu(theta02/(th0v_new))/2 + accu(theta02/(th0v))/2;
    testpa         = testpa - accu(log(th0v_new))/2 + accu(log(th0v))/2;
    
    // Prior for gamma_mu is truncated normal
    testpa = testpa - pow(gamma_mu_new - gamma_mu_m0, 2)/(2*gamma_mu_v02);
    testpa = testpa + pow(gamma_mu     - gamma_mu_m0, 2)/(2*gamma_mu_v02);
    // Prior for sigma is truncated normal
    testpa = testpa - pow(gamma_sigma_new - gamma_sigma_m0, 2)/(2*gamma_sigma_v02);
    testpa = testpa + pow(gamma_sigma     - gamma_sigma_m0, 2)/(2*gamma_sigma_v02);
    
    // Truncated random walk for mu
    pnew = normcdf(mu_top,gamma_mu,std_met_g0)     - normcdf(mu_bot,gamma_mu,std_met_g0);      // Normalizing constant for gamma_mu_new
    pold = normcdf(mu_top,gamma_mu_new,std_met_g0) - normcdf(mu_bot,gamma_mu_new,std_met_g0);  // Normalizing constant for gamma_mu_old
    testpa = testpa + log(pnew) - log(pold);
    // Truncated random walk for sigma
    pnew = normcdf(sigma_top,gamma_sigma,std_met_g0)     - normcdf(sigma_bot,gamma_sigma,std_met_g0);            // Normalizing constant for gamma_sigma_new
    pold = normcdf(sigma_top,gamma_sigma_new,std_met_g0) - normcdf(sigma_bot,gamma_sigma_new,std_met_g0);    // Normalizing constant for gamma_sigma_old
    testpa = testpa + log(pnew) - log(pold);  
    if(log(randu())<testpa){
      gamma_alpha  = gamma_alpha_new;
      gamma_beta   = gamma_beta_new;
      gamma_mu     = gamma_mu_new;
      gamma_sigma  = gamma_sigma_new;
      //gamma_sigma2 = gamma_sigma2_new;
      gamma_prob   = gamma_prob_new;
      gam0vec      = gam0vec_new;
      th0v         = th0v_new;
      //gamma0       = gamma_mu;
    }  // End generate gamma_alpha and gamma_beta
    // End HB Smoothing Prior
    
    //-------------------------------------------------------
    // If probit model, generate latent ydata where ydata0 are observed 0/1
    if(iflaglm){
      ydata(id_probit1)     = rcpp_rndtnb2(ymall(id_probit1),1.0,0.0);  // Y > 0
      ydata(id_probit0)     = rcpp_rndtnb2(ymall(id_probit0),1.0,0.0);  // Y < 0
    }
    
    if(iflagsc>0){
      // shape constrained: modify metropolis
      if ((imcmc < nblow0*maxmodmet) && (imcmc == floor((double)imcmc/(double)nblow0)*nblow0) && imcmc > 0) {
        pmetm = pmet/nblow0;
        pok=0;
        for (unsigned int j = 0; j < ngroup; j++) {
          if (pmetm(j) > 0.6){  // Too many acceptances, increase metm  
            metm(j) *= 10;
            if (disp) cout << j << "th group pmet = " << pmet(j)/nblow0 << " > 0.6. Increase metm and redo MCMC loop" << endl;
            mets(j)           = metm(j)/sqrt(met_alpha(j)-2);
            metv(j)           = pow(mets(j), 2);     // Variance of inverse gamma
            met_beta(j)       = (met_alpha(j)-1)*metm(j);
            met_beta_AM(j)    = met_beta(j); // updated IG parameter for adpative metroplis
            met_var_all(j)    = metm(j); // IG values                                    
            met_mAM(j)        = metm(j);
            met_vAM(j)        = metv(j);
            
          } else if (pmetm(j) < 0.3) {  // Too few acceptances, decrease metm 
            metm(j) /= 10;
            if (disp) cout << j << "th group pmet = " << pmet(j)/nblow0 << " < 0.3. Reduce metm and redo MCMC loop" << endl;
            mets(j)           = metm(j)/sqrt(met_alpha(j)-2);
            metv(j)           = pow(mets(j), 2);     // Variance of inverse gamma
            met_beta(j)       = (met_alpha(j)-1)*metm(j);
            met_beta_AM(j)    = met_beta(j); // updated IG parameter for adpative metroplis
            met_var_all(j)    = metm(j); // IG values
            met_mAM(j)        = metm(j);
            met_vAM(j)        = metv(j);
          }else{
            pok++;
          }
        }
        
        pmet.zeros();
        icount_AM.zeros();
        if ( pok == ngroup) {
          if (disp)  cout << "Got good acceptance rates: break out of loop" << endl; 
          imcmc = nblow0*maxmodmet;
#ifdef use_progress
          if (nblow0*maxmodmet - imcmc > 0) {
            p.increment(nblow0*maxmodmet - imcmc); 
          }
#endif
        }
      }
    }
    
    
    
    // Save MCMC iterations
    if ((imcmc >= nblow+nblow0*maxmodmet) && (imcmc % nskip == 0)) {
      if (nparv > 0) alphag.row(isave) = alpha.t();
      sigmag.row(isave)    = sigma.t();
      sigma2_alphag(isave) = sigma2_alpha;
      sigma2_betag(isave)  = sigma2_beta;
      tauallg.row(isave)   = tauall.t();
      eta0g(isave)         = eta0;
      gamma_mugg(isave)    = gamma_mu;
      gamma_alphag(isave)  = gamma_alpha;
      gamma_betag(isave)   = gamma_beta;
      gammallg.row(isave)  = gammall.t();
      thetam               += thetall;
      thetas               += square(thetall);
      theta0g.row(isave)   = theta0.t();
      fxobsm               += fxobs;
      fxobss               += square(fxobs);
      fxgridm              += fxgridall;
      fxgrids              += square(fxgridall);
      f0xgridm             += f0xgrid;
      f0xgrids             += square(f0xgrid);
      f0Bxgridm            += f0Bxgrid;
      f0Bxgrids            += square(f0Bxgrid);
      phig.row(isave)      = phi_vec.t();
      lambdag.row(isave)   = vectorise(lambda).t();
      betam                += betall;
      betas                += square(betall);
      
      // U-shaped has squish function
      if(iflagsc == 4){
        zetag.row(isave)  = zetall.t();
        omegag.row(isave) = omegall.t();
        zeta0g(isave)     = zeta0;
        omega0g(isave)    = omega0;
        if(iflagpsi){
          psiallg.row(isave) = psiall.t();
          psi0g(isave)       = psi0;
          
        }  // Estimated psi for squish function
      }  // Got squish function
      
      isave++;
    } // End save MCMC
  } // end of mcmc loop
  
  cout << "MCMC is done!" << endl;
  
  Rcpp::List squishg = Rcpp::List::create(Rcpp::Named("zetag")   = zetag,
                                          Rcpp::Named("omegag")  = omegag,
                                          Rcpp::Named("zeta0g")  = zeta0g,
                                          Rcpp::Named("omega0g") = omega0g,
                                          Rcpp::Named("psiallg") = psiallg,
                                          Rcpp::Named("psi0g")   = psi0g);
  
  mat fgridg (nint+1, 2*ngroup+4);
  fgridg.cols(0,        ngroup-1) = fxgridm;
  fgridg.cols(ngroup, 2*ngroup-1) = fxgrids;
  fgridg.col(2*ngroup)   = f0xgridm;
  fgridg.col(2*ngroup+1) = f0xgrids;
  fgridg.col(2*ngroup+2) = f0Bxgridm;
  fgridg.col(2*ngroup+3) = f0Bxgrids;
  
  mat fxobsg (ntot, 2);
  fxobsg.col(0) = fxobsm;
  fxobsg.col(1) = fxobss;
  
  
  mat betag (ngroup, 2*nparw);
  betag.cols(0,       nparw-1) = betam;
  betag.cols(nparw, 2*nparw-1) = betas;
  
  
  // maximum element of list is 20.
  return Rcpp::List::create(Rcpp::Named("pmet")          = pmet,
                            Rcpp::Named("alphag")        = alphag,
                            Rcpp::Named("sigmag")        = sigmag,
                            Rcpp::Named("sigma2_alphag") = sigma2_alphag,
                            Rcpp::Named("sigma2_betag")  = sigma2_betag,
                            Rcpp::Named("tauallg")       = tauallg,
                            Rcpp::Named("eta0g")         = eta0g,
                            Rcpp::Named("gamma_mugg")    = gamma_mugg,
                            Rcpp::Named("gamma_alphag")  = gamma_alphag,
                            Rcpp::Named("gamma_betag")   = gamma_betag,
                            Rcpp::Named("gammallg")      = gammallg,
                            Rcpp::Named("theta0g")       = theta0g,
                            Rcpp::Named("phig")          = phig,
                            Rcpp::Named("lambdag")       = lambdag,
                            Rcpp::Named("squishg")       = squishg,
                            Rcpp::Named("betag")         = betag,
                            Rcpp::Named("fgridg ")       = fgridg,
                            Rcpp::Named("fxobsg")        = fxobsg,
                            Rcpp::Named("thetam")        = thetam,
                            Rcpp::Named("thetas")        = thetas);
}



