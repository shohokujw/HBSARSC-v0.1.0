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
Rcpp::List get_mcmc_hbsarv5_nromal(arma::mat& data_all,
                                   arma::mat& vdata,
                                   arma::mat& wdata,
                                   arma::mat& zdata,
                                   arma::vec& id_group,
                                   arma::mat& iptw2,
                                   arma::mat& wtw,
                                   arma::mat& ztz,
                                   double xdelta,
                                   arma::colvec& xgrid,
                                   double xrange,
                                   arma::uvec& xinx,
                                   arma::vec& xover,
                                   arma::vec& idn0,
                                   arma::mat& met_vec_para,
                                   arma::vec& met_para,
                                   
                                   arma::colvec& phi0xgrid,
                                   arma::mat& phixgrid,
                                   arma::mat& phi2xgrid,
                                   
                                   arma::mat& fndfxgridall,
                                   arma::mat& f0,
                                   
                                   arma::mat& vtv,
                                   arma::mat& betall,
                                   arma::mat& phi,
                                   arma::mat& lambda,
                                   arma::mat& lambdai,
                                   arma::mat& thetall,
                                   arma::colvec& xiparll,
                                   arma::colvec& theta0,
                                   double xipar0,
                                   arma::colvec& gammall,
                                   arma::vec& gamma_para,
                                   arma::colvec& tau2all,
                                   arma::colvec& tauall,
                                   
                                   arma::colvec& alpha_m0,
                                   arma::mat& alpha_v0,
                                   arma::mat& alpha_v0i,
                                   arma::colvec& alpha_v0im0,
                                   arma::vec& sigma_para,
                                   arma::colvec& phi_b0,
                                   arma::mat& phi_v0i,
                                   arma::colvec& phi_v0ib0,
                                   double lambda_f0,
                                   double lambda_fn,
                                   arma::mat& lambda_g0,
                                   arma::mat& lambda_g0i,
                                   arma::vec& tau_para,
                                   arma::vec& eta_para,
                                   
                                   double abot,
                                   double labot,
                                   
                                   arma::vec& nobs1,
                                   Rcpp::List& iptHB1,
                                   arma::vec& cpara,
                                   
                                   bool iflagCenter,
                                   int  iflagpn,
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
  int icntgroup = ngroup;
  
  arma::vec metw        = met_vec_para.col(0);
  arma::vec metm        = met_vec_para.col(1);
  arma::vec met_alpha   = met_vec_para.col(2);
  arma::vec mets        = met_vec_para.col(3);
  arma::vec metv        = met_vec_para.col(4);
  arma::vec met_mAM     = met_vec_para.col(5);
  arma::vec met_vAM     = met_vec_para.col(6);
  arma::vec met_beta    = met_vec_para.col(7);
  arma::vec met_beta_AM = met_vec_para.col(8);
  arma::vec met_var_all = met_vec_para.col(9);
  
  bool iflagAM        = met_para(0);
  double met_alpha_g0 = met_para(1);
  double met_beta_g0  = met_para(2);
  
  arma::mat fxgridall  = fndfxgridall.cols(0, ngroup-1);
  arma::mat dfxgridall = fndfxgridall.cols(ngroup, fndfxgridall.n_cols-1);
  
  arma::colvec f0xgrid  = f0.col(0);
  arma::colvec f0Bxgrid = f0.col(1);
  
  arma::colvec kall  = linspace<colvec>(1, nbasis, nbasis);
  arma::colvec kall0 = linspace<colvec>(0, nbasis, nbasis+1);
  
  double gamma_alpha     = gamma_para(0);
  double gamma_beta      = gamma_para(1);
  double gamma_sigma     = gamma_para(2);
  double gamma_lmu       = gamma_para(3);
  double gamma_lsigma    = gamma_para(4);
  double gamma0          = gamma_para(5);
  double gamma_lmu_m0    = gamma_para(6);
  double gamma_lmu_v0    = gamma_para(7);
  double gamma_lsigma_m0 = gamma_para(8);
  double gamma_lsigma_v0 = gamma_para(9);
  
  double sigma2_s0 = sigma_para(0);
  double sigma2_rn = sigma_para(1);
  double sigma2    = sigma_para(2);
  double sigma     = sigma_para(3);
  
  double tau2_s0 = tau_para(0);
  double tau2_rn = tau_para(1);
  
  double eta_s0 = eta_para(0);
  double eta_rn = eta_para(1);
  double eta2   = eta_para(2);
  double eta    = eta_para(3);
  
  arma::colvec ydata  = data_all.col(0);
  arma::colvec valpha = data_all.col(1);
  arma::colvec wbeta  = data_all.col(2);
  arma::colvec fxobs  = data_all.col(3);
 
   // save mcmc result
  vec pmet             = zeros<vec> (ngroup);
  vec icount_AM        = zeros<vec> (ngroup);
  mat alphag           = zeros<mat> (smcmc,nparv);
  colvec sigmag        = zeros<colvec> (smcmc);
  mat tauallg          = zeros<mat> (smcmc,nbasis+1);
  colvec etag          = zeros<colvec> (smcmc);
  colvec gamma0g       = zeros<colvec> (smcmc);
  mat gammallg         = zeros<mat> (smcmc,ngroup);
  mat thetam           = zeros<mat> (nbasis+1,ngroup);
  mat thetas           = zeros<mat> (nbasis+1,ngroup);
  mat theta0g          = zeros<mat> (smcmc,nbasis+1);
  colvec fxobsm        = zeros<colvec> (ntot);
  colvec fxobss        = zeros<colvec> (ntot);
  mat fxgridm          = zeros<mat> (nint+1,ngroup);
  mat fxgrids          = zeros<mat> (nint+1,ngroup);
  mat dfxgridm         = zeros<mat> (nint+1,ngroup);
  colvec f0xgridm      = zeros<colvec> (nint+1);
  colvec f0xgrids      = zeros<colvec> (nint+1);
  colvec f0Bxgridm     = zeros<colvec> (nint+1);
  colvec f0Bxgrids     = zeros<colvec> (nint+1);
  mat phig             = zeros<mat> (smcmc,phidim);
  mat lambdag          = zeros<mat> (smcmc,nparw2);
  mat betam            = zeros<mat> (ngroup,nparw);
  mat betas            = zeros<mat> (ngroup,nparw);
  colvec gamma_alphag  = zeros<colvec> (smcmc);
  colvec gamma_betag   = zeros<colvec> (smcmc);
  colvec gamma_lmug    = zeros<colvec> (smcmc);
  colvec gamma_mug     = zeros<colvec> (smcmc);
  colvec gamma_lsigmag = zeros<colvec> (smcmc);
  colvec yresid        = zeros<colvec> (ntot);
  
  // parametric
  colvec alpha   = zeros<colvec> (nparv);
  mat vbni       = zeros<mat> (nparv,nparv);
  mat vbni12     = zeros<mat> (nparv,nparv);
  colvec v0ib0   = zeros<colvec> (nparw);
  double sse = 0.0;
  mat wtwj       = zeros<mat> (nparw, nparw);
  mat zbreg      = zeros<mat> (ngroup, nparz);
  colvec zbregj  = zeros<colvec> (nparz);
  double sigma2_sn = 0.0;
  mat phi_vec    = zeros<mat> (nparz,nparw);
  mat phi_vbni   = zeros<mat> (nparz*nparw,nparz*nparw);
  mat phi_vbni12 = zeros<mat> (nparz*nparw,nparz*nparw);
  mat ztb        = zeros<mat> (nparz,nparz);
  colvec a       = zeros<colvec> (nparz);
  colvec phi_bn  = zeros<colvec> (nparz*nparw);
  mat resid      = zeros<mat> (ngroup,nparw);
  mat lambda_gni = zeros<mat> (nparw, nparw);
  mat lambda_gn  = zeros<mat> (nparw, nparw);
  
  // theta & theta0
  // double theta00;
  colvec theta0k         = zeros<colvec> (nbasis);
  colvec thv             = zeros<colvec> (nbasis+1);
  double thv0            = 0.0;
  colvec thvk            = zeros<colvec> (nbasis);
  double ths0            = 0.0;
  colvec thsk            = zeros<colvec> (nbasis);
  colvec th0_thv0        = zeros<colvec> (nbasis+1);
  colvec th0_thvk        = zeros<colvec> (ngroup);
  colvec theta_old       = zeros<colvec> (nbasis+1);
  //double theta0_old;
  double xiparj_old      = 0.0;
  colvec thetak_old      = zeros<colvec> (nbasis);
  double t0resid2_old    = 0.0;
  colvec tkresid2_old    = zeros<colvec> (nbasis);
  double t0resid2_new    = 0.0;
  colvec tkresid2_new    = zeros<colvec> (nbasis);
  double met_var_all_new = 0.0;
  double ck              = 0.0;
  //double met_beta_new;
  double xiparj_new = 0.0;
  double theta0_new = 0.0;
  colvec thetak_new = zeros<colvec> (nbasis);
  colvec theta_new  = zeros<colvec> (nbasis+1);
  Rcpp::List outfx;
  double testp   = 0.0;
  double sse_new = 0.0;
  double sse_old = 0.0;
  double eta_sn  = 0.0;
  
  // compute f0
  colvec tvar0_new = zeros<colvec> (nbasis);
  colvec tvar0_old = zeros<colvec> (nbasis);
  colvec t02       = zeros<colvec> (nbasis); 
  colvec fpara     = zeros<colvec> (nbasis+1);
  colvec vpar      = zeros<colvec> (nbasis+1);
  colvec thvb      = zeros<colvec> (nbasis+1);
  
  // xi0
  double xi_vni   = 0.0;
  double xi_vn    = 0.0;
  double xi_bn    = 0.0;
  colvec resid_xi = zeros<colvec> (ngroup);
  
  // tau
  colvec resid_tau  = zeros<colvec> (ngroup);
  colvec resid2_tau = zeros<colvec> (ngroup);
  colvec tk         = zeros<colvec> (ngroup);
  double tau2_sn    = 0.0;;
  
  // eta
  colvec resid_eta = zeros<colvec> (nbasis+1);
  colvec vec0      = zeros<colvec> (nbasis+1);
  double etv0      = 0.0;
  
  // gamma
  double gamma_mu=0.0;
  colvec gamvec  = zeros<colvec> (nbasis);
  colvec gam0vec = zeros<colvec> (nbasis+1);
  double var_met_g0=0.0;
  double std_met_g0=0.0;
  double gamma_old=0.0;
  double lgamma_old=0.0;
  double gbot=0.0;
  double gtop=0.0;
  double lgamma_new=0.0;
  double gamma_new=0.0;
  colvec tvold        = zeros<colvec> (nbasis);
  colvec tvnew        = zeros<colvec> (nbasis);
  colvec resid_gamma  = zeros<colvec> (nbasis);
  colvec resid2_gamma = zeros<colvec> (nbasis);
  double p1=0.0;
  double p2=0.0;
  colvec lngammall = zeros<colvec> (ngroup);
  double s1=0.0;
  double s2=0.0;
  double lmbot=0.0;
  double gamma_lmu_new=0.0;
  double lstop=0.0;
  double gamma_lsigma_new=0.0;
  double gamma_mu_new=0.0;
  double gamma_sigma_new=0.0;
  double gamma_sigma2_new=0.0;
  double gamma_alpha_new=0.0;
  double gamma_beta_new=0.0;
  double gamma0_new=0.0;
  double lmbot_new=0.0;
  double lstop_new=0.0;
    
  double met_var0=0.0;
  double met_var=0.0; 
  double met_std0=0.0;
  double met_std=0.0; 
  
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
      colvec bn (nparv);
      yresid = ydata - wbeta - fxobs;   // residuals
      vbni = vtv / sigma2 + alpha_v0i;  // precision
      vbni12 = chol(vbni);
      bn = (vdata.t() * yresid) / sigma2 + alpha_v0im0;
      mat rMat = randn(nparv, 1);
      alpha = solve(vbni, bn + vbni12.t()*rMat);
      valpha = vdata * alpha;
    } // End generate alpha
    
    // Do HB parametric model
    yresid = ydata - fxobs - valpha;
    //--------------------------------------------
    // Generate Beta_j
    zbreg = zdata * phi;
    sse = 0;
    for (unsigned int j = 0; j < ngroup; j++) {
      uvec gidx = iptHB1[j];
      zbregj = zbreg.row(j).t();
      v0ib0 = lambdai * zbregj;
      mat wj (nobs1(j), nparw);
      for (int k = 0; k < nobs1(j); k++) {
        wj.row(k) = wdata.row(gidx(k));
      }
      colvec yj = yresid.elem(gidx);
      wtwj =  wtw.rows(iptw2(j, 0), iptw2(j, 1));
      vbni = wtwj / sigma2 + lambdai;
      vbni12 = chol(vbni);   // Upper triangluar: t(vbni12) %*% vbni12 = vbni
      colvec bn = (wj.t() * yj) / sigma2 + v0ib0;
      
      // Generate draw by using solve: avoides inverse of vni
      mat ranMat = randn(nparw, 1);
      colvec bj = solve(vbni, bn + vbni12.t()*ranMat);
      betall.row(j) = bj.t();
      // Compute residual
      colvec wjbj = wj * bj;
      colvec rj = yj - wjbj;
      sse = sse + accu(square(rj));
      for (int k = 0; k < nobs1(j); k++) {
        wbeta(gidx(k)) = wjbj(k);
      }
      //--------------------------------------------
    }
    
    //--------------------------------------------
    // Generate Sigma2, the error variance
    sigma2_sn = sigma2_s0 + sse;
    sigma2 = sigma2_sn / (2 * randg(distr_param(sigma2_rn / 2, 1.0)));
    sigma = sqrt(sigma2);
    
    //------------------------------------------
    // Generate Phi (Beta = Zdata*Phi + delta)
    phi_vbni = kron(lambdai, ztz) + phi_v0i;
    phi_vbni12 = chol(phi_vbni);
    ztb = zdata.t() * betall;
    for (unsigned int j = 0; j < nparw; j++) {
      if (nparw > 1) {
        a = ztb * lambdai.col(j);
      }
      else {
        a = ztb * lambdai(0, 0);
      }
      
      phi_bn.subvec(nparz*j, nparz*(j+1)-1) = a;
      
    }
    mat rMat = randn(phidim, 1);
    phi_vec = solve(phi_vbni, phi_bn + phi_vbni12.t()*rMat);
    for (unsigned int i = 0; i < nparw; i++) {
      phi.col(i) = phi_vec.rows(nparz*i, nparz*(i+1)-1);
    }
    
    //------------------------------------------
    // Generate Lambda from Beta = Zdata*breg + delta  
    resid = betall - zdata * phi;
    if (nparw > 1) {
      lambda_gni = lambda_g0i + resid.t()*resid;
      lambda_gn = inv(lambda_gni);
      lambdai = wishrnd(lambda_gn, lambda_fn); // random draw from a Wishart distribution
      lambda = inv(lambdai);
    }
    else {
      sse = accu(square(resid));
      lambda_gn = lambda_g0i + sse;
      lambda = lambda_gn + (2 * randg(distr_param(lambda_fn / 2, 1.0)));
      lambdai = 1 / lambda;
    }
    
    //------------------------------------------
    // Do HB nonparametric model
    yresid = ydata - valpha - wbeta;   // Residuals
    
    //theta00 = theta0(0);                 // Upper-level Fourier coefficient for constant
    theta0k = theta0.subvec(1, nbasis);  // Upper-level Fourier coefficient for consines
    
    // Loop over groups to generate theta_j
    for (unsigned int j = 0; j < ngroup; j++) {
      gamvec = exp(-gammall(j)*kall0);
      thv = tau2all % gamvec; // Schur product: element-wise multiplication of two objects
      
      // Need to separate theta_{j,0} from theta_{j,k} since theta_{g,0} > 0!
      thv0 = thv(0);                 // variance for xipar_j = log(theta_j0)
      thvk = thv.subvec(1, nbasis);  // variance for theta_k where k > 0
      ths0 = sqrt(thv0);             // std dev  for ln.theta_00 = xipar_j
      thsk = sqrt(thvk);             // std dev  for theta_k where k > 0
      
      int nj = nobs1(j);
      vec gidx = iptHB1[j];
      colvec yresidj(nj);
      colvec fxobs_old(nj);
      for (int k = 0; k < nj; k++) {
        yresidj(k)   = yresid(gidx(k));
        fxobs_old(k) = fxobs(gidx(k));
      }
      
      theta_old = thetall.col(j);                  // Old value for theta_j
      //theta0_old = theta_old(0);                   // Old value for theta_{j,0}
      xiparj_old = xiparll(j);                     // xipar_j = log(theta_j0)
      thetak_old = theta_old.subvec(1, nbasis);    // Old values of theta_{j,k>0}
      t0resid2_old = pow(xiparj_old - xipar0, 2);  // squared residual  k=0
      tkresid2_old = pow(thetak_old - theta0k, 2); // squared residuals k>0
      
      // get variances for t-dist random walk
      met_var_all_new = met_beta_AM(j) / randg(distr_param(met_alpha(j), 1.0));
      ck = metw(j)*met_var_all_new + (1 - metw(j))*met_mAM(j);
      //met_beta_new = (met_alpha(j) - 1)*ck;
      
      ck = 5.66;   // Constant from Harito
      met_var0 = ck * met_var_all_new;
      met_var  = ck * met_var_all_new;
      met_std0 = sqrt(met_var0);
      met_std  = sqrt(met_var);
      
      // Random walk from for normals: theta_j0 = exp(xiparj) for constant
      xiparj_new = xiparj_old + met_std0*ths0*randn();
      theta0_new = exp(xiparj_new);
      thetak_new = thetak_old + met_std*thsk%randn(nbasis);
      theta_new(0) = theta0_new;
      theta_new.subvec(1,nbasis) = thetak_new;
      
      // Get fJ(x) for theta_new
      outfx = rcpp_GetUpfxgrid(theta_new, phixgrid, xdelta, xrange, iflagCenter);
      
      colvec fxgrid_new = outfx["fx"];      // f computed on xgrid
      colvec dfxgrid_new = outfx["dfx"];     // derivative of f computed on xgrid
      
      t0resid2_new  = pow(xiparj_new - xipar0,2);
      tkresid2_new  = square(thetak_new - theta0k);
      
      // Compute fj(x_obs)
      // fxobs = GetUpfxobs(xgrid,xinx,xover,fxgrid)
      uvec a(nj);
      vec b(nj);
      for (int k = 0; k < nj; k++) {
        a(k) = xinx(gidx(k));
        b(k) = xover(gidx(k));
      }
      colvec fxobs_new = rcpp_GetSCfxobs(xgrid, a, b, fxgrid_new, iflagpn);
      
      // Metropolis Test 
      colvec resid_new = yresidj - fxobs_new;  // Note: fxobs_new_j only for group j
      sse_new = accu(square(resid_new));
      
      colvec resid_old = yresidj - fxobs_old;
      sse_old = accu(square(resid_old));
      
      // Compute log test probabability
      testp = -(sse_new - sse_old) / (2*sigma2);               // Likelihood         
      testp -= (t0resid2_new - t0resid2_old) / (2*thv0);       // Priors for xipar_j
      testp -= accu((tkresid2_new - tkresid2_old) / (2*thvk)); // Priors for theta_{j,k>0}
      
      if (log(randu()) < testp) { // Accept candidate
        xiparll(j)     = xiparj_new;
        thetall.col(j) = theta_new;
        for (int k = 0; k < nj; k++) {
          fxobs(gidx(k)) = fxobs_new(k);
        }
        fxgridall.col(j)  = fxgrid_new;
        dfxgridall.col(j) = dfxgrid_new;
        met_var_all(j)    = met_var_all_new;
        pmet(j)++;
      }
      
      if (iflagAM) { // Update mean, variance, alpha, and beta for metropolis
        icount_AM(j)++;
        met_mAM(j) = met_mAM(j) + (met_var_all(j) - met_mAM(j)) / (icount_AM(j));
        met_vAM(j) = ((icount_AM(j) - 1) / icount_AM(j))*met_vAM(j) + (pow(met_var_all(j) - met_mAM(j), 2)) / icount_AM(j);
      }
      
      ck = metw(j)*met_var_all(j) + (1 - metw(j))*met_mAM(j);
      met_beta_AM(j) = (met_alpha(j) - 1)*ck;
      
    } //  End loop over groups
    
    // -----------------------------------------------
    // Generate upper-level model theta_0k
    // theta_jk ~ N(theta_0k,tau_k^2*exp(-k*gamma_j)) for k>0
    // theta_0k ~ N(0,eta0^2*exp(-k*gamma_0))
    gam0vec  = pow(1 + kall0 / gamma_beta, -gamma_alpha);
    th0_thv0 = eta2 * gam0vec;           // Variance of theta_0k
    
    // Generate xipar_j = ln(theta_00) ~ N(xipar_0,tau_0^2)
    // xipar_0 ~ N(0,tau_0^2)
    xi_vni    = icntgroup / tau2all(0) + 1 / th0_thv0(0);
    xi_vn     = 1 / xi_vni;
    xi_bn     = xi_vn * (accu(xiparll) / tau2all(0));
    xipar0    = xi_bn + sqrt(xi_vn)*randn();
    theta0(0) = exp(xipar0);
    // Loop to generate theta_0k for frequency k >= 1
    for (unsigned int k = 0; k < nbasis; k++) {
      th0_thvk      = tau2all(k + 1)*exp(-(double)(k + 1)*gammall);   // Var(theta_jk) = tau_k^2*exp(-k*gamma_j)
      xi_vn         = 1 / (accu(1 / th0_thvk) + 1 / th0_thv0(k + 1));  // variance
      xi_bn         = xi_vn*accu(thetall.row(k + 1).t() / th0_thvk); // mean
      theta0(k + 1) = xi_bn + sqrt(xi_vn)*randn();
    } // End loop to generate theta_0k for k >= 1
    
    //-----------------------------------------------
    // Update smoothing parameters: variances
    // theta_0 =exp(xipar) and xipar is Normal 
    resid_xi   = xiparll - xipar0;
    tau2_sn    = tau2_s0 + accu(square(resid_xi));
    tau2all(0) = tau2_sn / (2 * randg(distr_param(tau2_rn/2, 1.0)));
    
    // Generate tau_k^2:   theta_{j,k} ~ N(theta_{0,k},tau_k^2*exp(-k*gamma_j))
    for (unsigned int k = 0; k < nbasis; k++) {
      tk             = thetall.row(k + 1).t();
      resid_tau      = tk - theta0(k + 1);
      resid2_tau     = square(resid_tau);
      tau2_sn        = tau2_s0 + accu(resid2_tau / exp(-(double)(k + 1)*gammall));
      tau2all(k + 1) = tau2_sn / (2*randg(distr_param(tau2_rn/2, 1.0)));
    } // End loop over frequencies for   
    
    // Get std dev
    tauall = sqrt(tau2all);
    
    // Variance of theta_0k
    vec0                        = pow(1 + kall0 / gamma_beta, -gamma_alpha);
    resid_eta(0)                = xipar0;
    resid_eta.subvec(1, nbasis) = theta0.subvec(1, nbasis);
    eta_sn                      = eta_s0 + accu(square(resid_eta) / vec0);
    eta2                        = eta_sn / (2*randg(distr_param(eta_rn/2, 1.0)));
    eta                         = sqrt(eta2);
    
    // Generate gamma for each population
    for (unsigned int j = 0; j < ngroup; j++) { // Loop over groups to generate gammaj
      var_met_g0 = met_beta_g0 / randg(distr_param(met_alpha_g0, 1.0));  // Variance for Metropolis
      std_met_g0 = sqrt(var_met_g0);                                // STD DEV for Metropolis
      gamma_old  = gammall(j);
      lgamma_old = log(gamma_old);
      //----------------------------------------------------
      // Generate gamma for group j
      // Caution! Sometimes gamma wandard off and we get NaNs
      // Bound gamma.  Hope it is not needed very often!
      gbot = -2.3;
      gtop = 2.3;
      lgamma_new      = rcpp_rndtnab(lgamma_old, std_met_g0, gbot, gtop);
      gamma_new       = exp(lgamma_new);
      tvold           = tau2all.subvec(1, nbasis)%exp(-gamma_old*kall);
      tvnew           = tau2all.subvec(1, nbasis)%exp(-gamma_new*kall);
      colvec thetallj = thetall.col(j);
      resid_gamma     = thetallj.subvec(1, nbasis) - theta0.subvec(1, nbasis);
      resid2_gamma    = square(resid_gamma);
      // Test probability
      // theta_jk ~ N(theta_0k,tau_k^2*exp(-k*gamma_j))
      testp = -accu(resid2_gamma / tvnew) / 2 - accu(log(tvnew)) / 2;
      testp += accu(resid2_gamma / tvold) / 2 + accu(log(tvold)) / 2;
      // "Prior" gamma_j ~ G(alpha,beta)
      testp += (gamma_alpha - 1)*(lgamma_new - lgamma_old) - gamma_beta*(gamma_new - gamma_old);
      // Include bounding probs
      p1 = normcdf((gtop - lgamma_old) / std_met_g0) - normcdf((gbot - lgamma_old) / std_met_g0);
      p2 = normcdf((gtop - lgamma_new) / std_met_g0) - normcdf((gbot - lgamma_new) / std_met_g0);
      testp += (log(p1) - log(p2));
      if (log(randu()) < testp) {
        gammall(j) = gamma_new;
      }
      // Ugh: if lgamma is close to gtop, MCMC is stuck.  Restart gamma
      if (gtop - lgamma_new < .1) gammall(j) = 1;
      
    }
    lngammall = log(gammall);
    // Summary stats
    s1 = accu(lngammall);
    s2 = accu(gammall);
    // ------------------------------------------------------------------------
    // Generate alpha and beta for gamma_j ~ G(alpha,beta)
    // Change parameters to mu = alpha/beta and sigma = sqrt(alpha)/beta
    // Random walk on log(mu) and log(sigma)
    // Bound alpha > abot.  So mu > abot/beta and sigma > sqrt(abot)/beta
    
    var_met_g0       = met_beta_g0/randg(distr_param(met_alpha_g0, 1.0));  // Variance for Metropolis
    std_met_g0       = sqrt(var_met_g0);                          // STD DEV for Metropolis
    // Generate mu given sigma old
    // alpha > abot means mu/sigma > sqrt(abot) or log(mu) > log(sigma) + log(sqrt(abot))
    lmbot            = gamma_lsigma + labot;
    gamma_lmu_new    = rcpp_rndtnb(gamma_lmu,std_met_g0,lmbot);    // candidate for log(mu)
    // Generate sigma given mu_new
    // sigma < mu/sqrt(a)
    lstop            = gamma_lmu - labot;
    gamma_lsigma_new = rcpp_rndtna(gamma_lsigma,std_met_g0,lstop); // candidate for log(sigma)
    gamma_mu_new     = exp(gamma_lmu_new);
    gamma_sigma_new  = exp(gamma_lsigma_new);
    gamma_sigma2_new = pow(gamma_sigma, 2);
    gamma_alpha_new  = pow(gamma_mu_new, 2)/gamma_sigma2_new;
    gamma_beta_new   = gamma_mu_new/gamma_sigma2_new;
    gamma0_new       = gamma_mu_new;
    lmbot_new        = gamma_lsigma_new + labot;
    lstop_new        = gamma_lmu_new    - labot;
    
    //  Compute test probability
    // loglikelihood ratio for gamma_j ~ G(alpha,beta)
    testp = icntgroup*(gamma_alpha_new*log(gamma_beta_new) - gamma_alpha*log(gamma_beta) 
                      - lgamma(gamma_alpha_new) + lgamma(gamma_alpha));
    testp += (s1*(gamma_alpha_new - gamma_alpha) - s2*(gamma_beta_new-gamma_beta));
    
    // Prior distribution: log(mu) and log(sigma) are normal
    testp -= (pow(gamma_lmu_new   -gamma_lmu_m0, 2)    - pow(gamma_lmu   -gamma_lmu_m0, 2))/(2*gamma_lmu_v0);
    testp -= (pow(gamma_lsigma_new-gamma_lsigma_m0, 2) - pow(gamma_lsigma-gamma_lsigma_m0, 2))/(2*gamma_lsigma_v0);
    // Normalizing constant for bounded random walk for gamma_lmu
    p1    = normcdf((gamma_lmu     -lmbot)/std_met_g0);
    p2    = normcdf((gamma_lmu_new -lmbot_new)/std_met_g0);
    if (p1 < 0.00001) p1 = 0.00001;
    if (p2 < 0.00001) p2 = 0.00001;
    testp += (log(p1) - log(p2));
    
    // Normalizing constant for bounded random walk for gamma_lsigma
    p1    = normcdf((lstop     - gamma_lsigma)/std_met_g0);
    p2    = normcdf((lstop_new - gamma_lsigma_new)/std_met_g0);
    if(p1 < 0.00001) p1 = 0.00001;
    if(p2 < 0.00001) p2 = 0.00001;
    testp += (log(p1) - log(p2));
    // Density for upper-level model of theta_0k
    // tvar0_new  = eta2*(1/1+kall/gamma_beta_new)^(gamma_alpha_new)
    // tvar0_old  = eta2*(1/1+kall/gamma_beta)^(gamma_alpha)
    tvar0_new  = eta2*pow(1+kall/gamma_beta_new, -gamma_alpha_new);
    tvar0_old  = eta2*pow(1+kall/gamma_beta,     -gamma_alpha);
    t02        = square(theta0.subvec(1, nbasis));
    testp      -= (accu(t02/(tvar0_new))/2 + accu(log(tvar0_new))/2);
    testp      += (accu(t02/(tvar0_old))/2 + accu(log(tvar0_old))/2);
    
    if (log(randu()) < testp) {
      gamma_alpha  = gamma_alpha_new;
      gamma_beta   = gamma_beta_new;
      gamma0       = gamma0_new;
      gamma_mu     = gamma_mu_new;
      gamma_sigma  = gamma_sigma_new;
      gamma_lmu    = gamma_lmu_new;
      gamma_lsigma = gamma_lsigma_new;
    }  // End generate gamma_alpha and gamma_beta
    
    //--------------------------------------------------------------------------
    // Compute f0.  Note: some issues with theta_j0 = exp(xipar_j) is log normal
    fpara      = theta0;                      // Upper level spectral coefficients
    fpara(0)   = fpara(0)*exp(tau2all(0)/2);  // log normal theta_{0,j}
    
    // Integrate exp(-k*gamma_j) where gamma_j ~ G(alpha,beta)
    vpar             = pow(1/(1+kall0/gamma_beta), gamma_alpha);
    thvb             = tau2all%vpar;
    etv0             = exp(tau2all(0));
    thvb(0)          = exp(2*xipar0 + tau2all(0))*(etv0-1);  // lognormal variance V(theta_j0);
    f0xgrid          = rcpp_GetUpf0xgrid(fpara,thvb,phixgrid,phi2xgrid,xdelta,xrange,iflagCenter);
    f0Bxgrid         = rcpp_GetUpf0Biasxgrid(fpara,phixgrid,xdelta,xrange,iflagCenter); // added
    
    // If during "adjustment period" check pmet
    
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
    
    // Save MCMC iterations
    if ((imcmc >= nblow+nblow0*maxmodmet) && (imcmc % nskip == 0)) {
      if (nparv > 0) alphag.row(isave) = alpha.t();
      sigmag(isave)        = sigma;
      tauallg.row(isave)   = tauall.t();           
      etag(isave)          = eta;
      gamma0g(isave)       = gamma0;
      gamma_alphag(isave)  = gamma_alpha;
      gamma_betag(isave)   = gamma_beta;
      gamma_lmug(isave)    = gamma_lmu;
      gamma_mug(isave)     = gamma_mu;
      gamma_lsigmag(isave) = gamma_lsigma;
      gammallg.row(isave)  = gammall.t();
      thetam               += thetall;
      thetas               += square(thetall);
      theta0g.row(isave)   = theta0.t();
      fxobsm               += fxobs;
      fxobss               += square(fxobs);
      fxgridm              += fxgridall;
      fxgrids              += square(fxgridall);
      dfxgridm             += dfxgridall;
      f0xgridm             += f0xgrid;
      f0xgrids             += square(f0xgrid);
      f0Bxgridm            += f0Bxgrid;
      f0Bxgrids            += square(f0Bxgrid);
      phig.row(isave)      = phi_vec.t();
      lambdag.row(isave)   = vectorise(lambda).t();  
      
      betam           += betall;
      betas           += square(betall);
      isave++;
    } // End save MCMC
  } // end of mcmc loop
  
  cout << "MCMC is done!" << endl;
  
  mat gammag (smcmc, 6+ngroup);
  gammag.col(0)             = gamma0g;
  gammag.col(1)             = gamma_alphag;
  gammag.col(2)             = gamma_betag;
  gammag.col(3)             = gamma_lmug;
  gammag.col(4)             = gamma_mug;
  gammag.col(5)             = gamma_lsigmag;
  gammag.cols(6,6+ngroup-1) = gammallg;
  
  mat fgridg (nint+1, 3*ngroup+4);
  fgridg.cols(0,          ngroup-1) = fxgridm;
  fgridg.cols(ngroup,   2*ngroup-1) = fxgrids;
  fgridg.cols(2*ngroup, 3*ngroup-1) = dfxgridm;
  fgridg.col(3*ngroup)   = f0xgridm;
  fgridg.col(3*ngroup+1) = f0xgrids;
  fgridg.col(3*ngroup+2) = f0Bxgridm;
  fgridg.col(3*ngroup+3) = f0Bxgrids;
  
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
                            Rcpp::Named("tauallg")       = tauallg,
                            Rcpp::Named("etag")          = etag,
                            Rcpp::Named("gammag")        = gammag,
                            Rcpp::Named("thetam")        = thetam,
                            Rcpp::Named("thetas")        = thetas,
                            Rcpp::Named("theta0g")       = theta0g,
                            Rcpp::Named("fgridg ")       = fgridg,
                            Rcpp::Named("fxobsg")        = fxobsg,
                            Rcpp::Named("phig")          = phig,
                            Rcpp::Named("lambdag")       = lambdag,
                            Rcpp::Named("betag")         = betag);
}