"hbsarv7" <- function(y, w, x, v, z, nbasis=30,id_group,
                      nint=200, mcmc = list(), prior = list(),xmin,xmax, 
                      shape = "", iflagCenter=FALSE, iflagpsi=TRUE, disp=FALSE) {
  
  # check arguments
  if (missing(y)) stop("y is required")
  if (missing(x)) stop("x is required")
  
  if (missing(w)) w = NULL
  if (missing(v)) v = NULL
  if (missing(z)) z = NULL
  
  
  shape_type = c("Free", "Increasing", "Decreasing", "InvertedU", "Ushape")
  shape <- match.arg(shape,shape_type,several.ok=TRUE)
  
  
  # iflagsc gives free or shaped constraint
  # 0 is free or unconstrainded
  # 1 is monotone
  # 2 is U shaped
  #------------------------------------------------------------------------------
  # iflagpn gives increaseing or decreasing
  #  1 for increasing
  # -1 for decreasing
  if (shape == "Free") {
    iflagsc = 0; iflagpn = 1
  } else if (shape == "Increasing") {
    iflagsc = 1; iflagpn = 1
  } else if (shape == "Decreasing") {
    iflagsc = 1; iflagpn = -1
  } else if (shape == "InvertedU") {
    iflagsc = 4; iflagpn = 1
  } else if (shape == "Ushape") {
    iflagsc = 4; iflagpn = -1
  }
  
  iflagZ = FALSE;
  
  #-------------------------------------------------------------------------------
  # Iflaglm gives different glm models
  # 0 is normal regression
  # 1 is binary probit
  #------------------------------------------------------------------------------
  # flaglm_type = c("normal ", "binary")
  # flaglm <- match.arg(flaglm,shape_type,several.ok=TRUE)
  iflaglm = 0
  
  nparx = 0
  #------------------------------------------------------------------------------  
  if (missing(xmin)) xmin = min(x)    # xmin < X < xmax.  Change as needed
  if (missing(xmax)) xmax = max(x)    # Change as needed
  if(iflagsc>0 & nparx==0) nparx = 1  # make sure you get the constant
  
  mcvals=list(nblow0=5000,nblow=10000,smcmc=1000,nskip=10,ndisp=1000,maxmodmet=10)
  mcvals[names(mcmc)]=mcmc
  nblow0    = mcvals$nblow0
  nblow     = mcvals$nblow
  smcmc     = mcvals$smcmc
  nskip     = mcvals$nskip
  ndisp     = mcvals$ndisp
  maxmodmet = mcvals$maxmodmet
  
  nmcmc     = nblow + nskip*smcmc
  
  if(iflagsc==0){
    # shape constraints: use adaptive Metropolis
    maxmodmet = 0  # Maximum times to modify Metropolis
    nblow0    = 0  # Initial number of initial MCMC parmeters
  }
  
  ntheta    = nbasis + nparx
  hfun_nsim = 500  # Simulated E[h(x|omega)]
  nmcmcall  = maxmodmet*nblow0 + nmcmc
  
  #-----------------------------------------------------
  ydata  = y
  xdata  = x
  wdata  = w
  vdata  = v
  zdata  = z
  ntot   = length(y)
  
  yname = colnames(ydata)
  xname = colnames(xdata)
  wname = colnames(wdata)
  vname = colnames(vdata)
  zname = colnames(zdata)
  
  # Break it up
  if(iflaglm == 1){ # Probit Model
    ydata0 = ydata           # ydata0 are the observed 0/1 variables
    # Initialize the latent y, which is N(mean,1)
    ydata = -1 + 2*ydata0   # ydata are the latent variables
    id_probit1 = which(ydata0==1)  # index for ones
    id_probit0 = which(ydata0==0)  # index for zeros
  } else {
    id_probit1 = rep(0, ntot)
    id_probit0 = rep(0, ntot)
  }
  
  if (is.null(wdata)) {
    wdata = matrix(1, ntot, 1)
  } else {
    wdata = cbind(1, wdata)
  }
  nparw  = ncol(wdata)
  nparw2 = nparw*nparw
  
  if (missing(vdata) || is.null(vdata)) {
    vdata  = matrix(0, ntot, 1)
    vnames = "None"
    vtv    = 1
    nparv  = 0
  } else {
    vdata  = as.matrix(vdata)
    vtv    = crossprod(vdata)
    nparv  = ncol(vdata) 
  }
  
  if (!is.numeric(id_group)) {
    id_group.org = id_group
    id_group     = as.integer(as.factor(id_group))
    ugroups.org  = sort(unique(id_group.org))
  } else {
    id_group.org = id_group
    ugroups.org  = sort(unique(id_group))
  }
  
  ugroups    = sort(unique(id_group))  # unique groups. Note: unique does not sort
  ngroup     = length(ugroups)
  
  if (missing(zdata) || is.null(zdata)) {
    zdata=matrix(1.0,nrow=ngroup,ncol=1)
    znames='CNST'
    colnames(zdata)=znames
    ztz = crossprod(zdata)
    # ndimz=0 # number of Z variables, excluding the intercept
  }else{
    if(!is.matrix(zdata)) zdata=as.matrix(zdata)
    znames=colnames(zdata)
    if(is.null(znames)) znames=paste('z',1:ncol(zdata),sep='')
    zdata=cbind(1.0,zdata)
    znames=c('CNST',znames)
    colnames(zdata)=znames
    ztz = crossprod(zdata)
    # ndimz=ncol(zdata)-1
  }
  nparz      = NCOL(zdata)
  
  if(ngroup != length(ugroups)){
    cat("\nError with group ids in Zdata and 2nd column of datall\n")
  }
  
  # id_group may not be sequential from 1 ... ngroup
  # Need an index to get correct betam
  id_beta  = matrix(0,ntot,1)
  for(j in 1:ngroup){
    id = which(id_group == ugroups[j])
    id_beta[id] = matrix(j,length(id),1)
  }
  
  #---------------------------------------------------
  outpt  = MakeHBPointer(id_group)
  nobs1  = outpt$nobs
  iptHB1 = outpt$iptd
  
  iptHB2 = vector("list", length = length(iptHB1)) # index for c++
  for (i in 1:length(iptHB1)) {
    iptHB2[[i]] = iptHB1[[i]]-1  
  }
  
  #----------------------------------------------------------
  # Grid on [xmim,xmax] for plotting functions
  # Also, grid for integration
  a      = floor(min(xdata))
  if(a<xmin) xmin = a
  b      = ceiling(max(xdata))
  if(b>xmax) xmax = b
  xdelta = (xmax-xmin)/nint
  xgrid  = seq(xmin,xmax,xdelta)
  xgrid  = matrix(xgrid,nint+1,1)
  xrange = xmax-xmin
  
  #------------------------------------------------
  # Computational strategy with shape constraints
  #   Instead of computing f(x_i) at observations,
  #   find where x_i falls in xgrid and use fxgrid
  #   Find location of xdata in xgrid
  #   xgrid[xinx[i]] < xi <= xgrid[xinx[i]+1]
  xinx   = matrix(1,ntot,1)
  for(i in 1:ntot){  # Loop over observations
    xi   = xdata[i]
    if(xi <= xgrid[2]){
      xinx[i] = 1
    } else if(xi == xmax){
      xinx[i] = nint+1
    }else{
      xinx[i] = max(which(xi>xgrid))
    }
  }
  # Get excess for that xdata is over boundy of xgrid based on xinx
  # Used in numerical integration
  xover = xdata - xgrid[xinx]
  idn0  = which(xover>0)  # which xdata is not at grid
  
  # *************************************************************
  # Use Metropolis-Hastings with shape constraints
  # * Random walk Metropolis-Hastings with t-errors for fmodel > 1
  # *  pmet = P(accept candidate theta)
  # * If pmet is > 0.7, then step size is too small: increase metm 
  # * If pmet is < 0.3, then step size is too large: decrease metm 
  # * Keep mets > metm  (std dev > mean)
  # * Why?  If mets > metm, then alpha < 3
  # * Also IG has a dead zone near zeros if alpha > 2.25
  # * So, you want to keep 2 < alpha < 2.25.
  # * Fix alpha.  Focus on metm to get mets and beta
  # **************************************************************
  metw      = rep(.5,ngroup);       #  weight*(IG factor) + (1-weight)*(Running Mean)
  metm      = rep(0.001,ngroup);    #  Mean of inverse gamma
  met_alpha = rep(4,ngroup);        #  Keep alpha fixed
  
  mets      = metm/sqrt(met_alpha-2);
  metv      = mets^2;         # Variance of inverse gamma 
  
  met_mAM     = metm;        # Adaptive Metropolis mean
  met_vAM     = metv;        # Adaptive Metropolis variance
  
  met_beta    = (met_alpha-1)*metm    # Beta from mean parameter
  met_beta_AM = met_beta              # Beta for adaptive metropolis
  met_var_all = metm;   #@ IG values                                    @
  iflagAM     = 1                     # 0/1 flag for adpative metropolis
  pmet        = rep(0,ngroup);
  icount_AM   = rep(0,ngroup);
  
  # Metropolis for gamma_0 in upper level model
  # STD DEV or "step size" is IG
  metm_g0      = 0.001              # mean of metropolis step size
  mets_g0      = 0.01                # std dev of metroplis step size
  met_alpha_g0 = 2 + (metm_g0/mets_g0)^2
  met_beta_g0  = metm_g0*(met_alpha_g0 -1)
  
  # Metropolis for zeta 
  metm_zeta      = 0.001              # mean of metropolis step size
  mets_zeta      = 0.01                # std dev of metroplis step size
  met_alpha_zeta = 2 + (metm_zeta/mets_zeta)^2
  met_beta_zeta  = metm_zeta*(met_alpha_zeta -1)
  
  # Metropolis for psi 
  metm_psi      = 0.5            # mean of metropolis step size
  mets_psi      = 2              # std dev of metroplis step size
  met_alpha_psi = 2 + (metm_psi/mets_psi)^2
  met_beta_psi  = metm_psi*(met_alpha_psi -1)
  
  
  # Metropolis for eta2 
  metm_eta2      = 0.001            # mean of metropolis step size
  mets_eta2      = 0.01              # std dev of metroplis step size
  met_alpha_eta2 = 2 + (metm_eta2/mets_eta2)^2
  met_beta_eta2  = metm_eta2*(met_alpha_eta2 -1)
  
  # Metropolis for sigma 
  metm_sigma      = 0.001            # mean of metropolis step size
  mets_sigma      = 0.1              # std dev of metroplis step size
  met_alpha_sigma = 2 + (metm_sigma/mets_sigma)^2
  met_beta_sigma  = metm_sigma*(met_alpha_sigma -1)
  
  # Vectors of frequencies for basis functions
  kall      = matrix(1:nbasis)   # frequences 1, 2, ..., nbasis for cosine functions
  lkall     = log(1+kall)
  ktheta    = kall   
  # Pad first entries with zeros
  if(nparx>0){
    ktheta = c(rep(0,nparx),ktheta)      # includes 0 for shape constraint
  } 
  
  #---------------------------------------------------
  # Compute Cosine functions on xgrid and xdata 
  phi0xgrid = matrix(1/sqrt(xrange),nint+1,1)
  phixgrid  = sqrt(2/xrange)*cos(pi*xgrid%*%(t(kall)-xmin)/xrange)
  # Add parametric functions to the cosine basis
  if(iflagsc>0){  # Got shape constraints
    # Add constant function to phi with shape constraints
    parxgrid = phi0xgrid
    if(nparx > 1){  # Got other parametric functions
      for(k in 1:(nparx-1)){
        parxgrid = cbind(parxgrid,xgrid^(k/2))
      }
    }
    phixgrid  = cbind(parxgrid,phixgrid)
    phi2xgrid = phixgrid^2    # Cosine squared functions on xgrid
    xipar0    = 0
  }else{
    # Free.  No constant, but may have parametric functions for mu(x)
    if(nparx > 0) {
      for(k in 1:nparx){
        phixgrid = cbind(zgrid^k,phixgrid)
      }
    } 
    phi2xgrid = phixgrid^2  
  }
  #----------------------------------------------------------------
  # If free function, compute phi at xdata
  if(iflagsc==0){
    phixdata  = sqrt(2/xrange)*cos(pi*xdata%*%(t(kall)-xmin)/xrange)
    if(nparx > 0){
      for(k in 1:nparx){
        phixdata = cbind(xdata^k,phixdata)
      }
    }
  }else{
    phixdata = as.matrix(0)
  }
  
  #---------------------------------------------------------
  # Compute v'v, w'w, and phi'phi for each population
  if(nparv>0){	
	  vtv = array(0,dim=c(nparv,nparv,ngroup))
  } else {
	  vtv = array(0,dim=c(1,1,ngroup))
  }
  
  wtw  = array(0,dim=c(nparw,nparw,ngroup))
  phi2 = array(0,dim=c(nbasis,nbasis,ngroup))
  for(j in 1:ngroup){
    wj    = wdata[iptHB1[[j]],]
    wj    = as.matrix(wj)
    wtw[,,j] = crossprod(wj)
    if(nparv>0){
      vj    = as.matrix(vdata[iptHB1[[j]],])
      vtv[,,j] = crossprod(vj)      
    } else {
	    vtv[,,j] = matrix(0, 1,1)	
	  }
    if(iflagsc==0){
      phixj     = phixdata[iptHB1[[j]],]
      phi2[,,j] = crossprod(phixj)
    }
    
  }
  
  #-----------------------------------------------------
  # Initialize priors
  # Fixed effect alpha ~ N(m0,V0)
  if(nparv>0){
    alpha_m0    = matrix(0,nparv,1)
    alpha_v0    = 1000*diag(nparv)
    alpha_v0i   = diag(nparv)/1000
    alpha_v0im0 = alpha_v0i %*% alpha_m0
  }else{  # no fixed effects: use placeholders
    alpha_m0    = matrix(0, 1,1)
    alpha_v0    = matrix(1000, 1,1)
    alpha_v0i   = matrix(1/1000, 1,1)
    alpha_v0im0 = matrix(0, 1,1)
  }
  
  # HB model for variance: sigma_j^2 ~ IG(alpha/2,beta/2)
  # Hyperpriors:   alpha ~ N(m0,v02)I(alpha>2) and beta~G(r0/2,s0/2)
  sigma2_alpha     = 8
  sigma2_beta      = 2
  sigma2_alpha_m0  = 8
  sigma2_alpha_v02 = 25
  sigma2_beta_r0   = 8
  sigma2_beta_s0   = 2
  sigma2_beta_rn   = sigma2_beta_r0  + ngroup
  sigma2           = matrix(1,ngroup,1)
  sigma            = sqrt(sigma2)
  
  
  # HB Regression beta = zdata*Phi + delta
  phidim    = nparw*nparz
  phi_b0    = matrix(0,phidim,1)    # prior mean 
  phi_v0i   = diag(phidim)/10000    # prior precision
  phi_v0ib0 = matrix(0,phidim,1)
  
  # Var(delta) = Lambda ~ IW
  lambda_f0  = 5
  lambda_fn  = lambda_f0 + ngroup
  lambda_g0  = 10*diag(nparw)
  lambda_g0i = diag(nparw)/10
  
  # theta_jk ~ N(theta_0k, tau_k^2*exp(-k*gamma))
  # tau_k^2  ~ IG(r0/2,s0/2)
  tau2_mu        = 1   # Note: Keep tau2_mu = 1
  tau2_sigma     = 10
  tau2_r0        = 2*((tau2_mu/tau2_sigma)^2 + 2)
  tau2_s0        = 2*tau2_mu*(tau2_r0/2 - 1)
  
  tau2_rn    = tau2_r0 + ngroup
  
  #-----------------------------------------------------
  w0   = 1
  wk   = sum(kall)/2
  wlk  = sum(lkall)/2
  
  
  #------------------------------------------------------
  # Initialize parameters
  if (nparv > 0) {
    alpha         = matrix(0,nparv,1) 
  } else {
    alpha         = matrix(0,1,1)
  }
  betall          = matrix(0,ngroup,nparw)
  phi             = matrix(0,nparz,nparw)
  lambda          = diag(1,nparw,nparw)
  lambdai         = diag(1,nparw,nparw)
  thetall         = matrix(0,ntheta,ngroup)    # lower level spectral coefficients
  theta0          = matrix(0,ntheta,1)     # upper level spectral coefficients
  xiparall        = 0
  xipar0          = 0
  if(iflagsc > 0){
    xiparall        = matrix(-.5,ngroup,1)       # xipar_j0 = log(theta_j0)
    xipar0          = -.5                        # xipar_00 = log(theta_00)
    thetall[1,]     = matrix(exp(xiparall),1,ngroup) 
    theta0[1,]      = exp(xipar0)
  }
  
  hfunall = matrix(1,2,ngroup)   # Squish Function
  hfunj_new = 1
  hfunj_old = 1
  
  # Parameters of squish function
  # if(iflagsc == 4)
  #---------------------------------------------------------
  # U-shaped has squish function 
  # h(x) = (1-exp(psi*(x-omega)))/(1+exp(psi*(x-omega))
  # HB model for omega_j
  # omegaj = xmin + xrange*exp(zetaj)/(1+exp(zetaj))
  # zetaj       ~ N(zeta_mu,zeta_sigma2)
  # zeta_mu     ~ N(zeta_m0,zeta_v02)
  # zeta_sigma2 ~ IG(zeta_r0/2,zeta_s0/2)
  zeta_mu       = log(.4/.6)
  zeta_sigma2   = 1
  zeta_sigma    = sqrt(zeta_sigma2)
  zeta_m0       = log(.4/.6)
  zeta_v02      = 9
  zeta_v0       = sqrt(zeta_v02)
  zeta_r0       = 4
  zeta_s0       = 2
  zeta_rn       = zeta_r0 + ngroup
  zetall        = matrix(zeta_m0,ngroup,1)
  ezall         = exp(zetall)
  omegall       = xmin + xrange*ezall/(1+ezall)
  zeta0         = zeta_mu
  ez0           = exp(zeta0)
  omega0        = xmin + xrange*ez0/(1+ez0)
  #---------------------------------------------------------
  # HB model for log(psi_j) (if iflagpsi = 1)
  # psij_log       ~ N(psi_log_mu,psi_log_sigma2)
  # psi_log_mu     ~ N(psi_log_m0,psi_log_v02)
  # psi_log_sigma2 ~ IG(psi_log_r0/2,psi_log_s0/2)
  psi_fixed      = 100                         # Use psi_fixed if iflagpsi = 0
  psi_log_fixed  = log(psi_fixed)
  psi0           = psi_fixed
  psi0_log       = psi_log_fixed
  psiall         = matrix(psi_fixed,ngroup,1)  # All values of psi in lower leve model
  psiall_log     = log(psiall)                 # Do HB & metropolis on log(psi)
  psi_log_mu     = psi_log_fixed
  psi_log_sigma2 = 1
  psi_log_sigma  = sqrt(psi_log_sigma2)
  psi_log_m0     = psi_log_fixed
  psi_log_v02    = 9
  psi_log_v0     = sqrt(psi_log_v02)
  psi_log_r0     = 5
  psi_log_s0     = 1
  psi_log_rn     = psi_log_r0 + ngroup
  
  
  hfunall  = matrix(1,nint+1,ngroup)
  for(j in 1:ngroup){
    hfunall[,j] = rcpp_GetSquish(omegall[j],psiall[j],xgrid)
  }  
  
  
  #---------------------------------------------------------------
  # Smoothing parameter  gamma
  gammall         = matrix(1,ngroup,1)       # lower level smoothing
  
  # Initialize parameters for gamma_j ~ G(alpha,beta)I(gamma < gmax)
  # Reparametrized to achieve better mixing
  # mu = alpha/beta and sigma2 = alpha/beta2
  # Priors: beta ~ Gamma(r0,s0) and alpha ~ N(alpha_m0,alpha_v02)I(0<alpha)
  gmax             = 3       # Maximum value for gamma
  gamma_mu         = 1
  gamma_sigma      = 1
  gamma_sigma2     = gamma_sigma^2
  gamma_beta       = gamma_mu/gamma_sigma2
  gamma_alpha      = gamma_mu*gamma_beta
  gamma_prob       = pgamma(gmax,shape=gamma_alpha,rate=gamma_beta)
  # hyper-parameters for gamma_mu and gamma_sigma
  # gamma_mu ~ Truncated N(m0,v02)
  gamma_mu_m0      = 1
  gamma_mu_v02     = 4
  mu_bot           = .1
  mu_top           = gmax
  
  # gamma_sigma ~ Truncated N(m0,v02)
  gamma_sigma_m0   = 1
  gamma_sigma_v02  = 16
  sigma_bot        = .5
  sigma_top        = 10
  
  #------------------------------------------------------------------
  tau2all         = matrix(1,ntheta,1)   # upper level variances for theta_jk
  tauall          = sqrt(tau2all)
  
  #------------------------------------------------------------------
  # Upper-level smoothing parameters
  # theta_0k ~ N(0,eta0^2*exp(-gamma0*log(1+k))) for k > nparx
  eta02           = 1                         # variance parameter for theta_0k, k > 0
  eta0            = sqrt(eta02)
  eta02_r0        = 4
  eta02_s0        = 4
  eta02_rn        = eta02 + nbasis
  # gamma0   ~ Gamma(alpha,beta)
  gamma0          = 1                                       
  gam0vec         = (1+kall/gamma_beta)^(-gamma_alpha)
  th0v            = eta02*gam0vec
  
  
  # Prior for theta_0k ~ N(0,v02) for k <= nparx
  theta0_v0       = 1
  theta0_v02      = theta0_v0^2
  
  wbeta           = matrix(0,ntot,1)
  for(j in 1:ngroup){   # Compute w*beta
    wj    = as.matrix(wdata[iptHB1[[j]],])
    bj    = as.matrix((betall[j,]))
    wbeta[iptHB1[[j]],1] = wj %*% bj
  }
  
  if(nparv>0){
    valpha = vdata%*%alpha
  } else {
    valpha  = matrix(0,ntot,1)
  }
  
  fxgridall       = matrix(0,nint+1,ngroup)
  f0xgrid         = matrix(0,nint+1,1)
  f0Bxgrid        = f0xgrid
  fxobs           = matrix(0,ntot,1)
  fxobsm          = matrix(0,ntot,1)  

  ## convert
  data_all = list(ydata=ydata, valpha=valpha, wbeta=wbeta, vdata=vdata, wdata=wdata, zdata=zdata)
  cpara = c(ntot, ngroup, nparv, nparw, nparw2, nparz, phidim, nbasis, nblow, nblow0, maxmodmet,
            nskip, nmcmcall, smcmc, nint, ntheta, nparx, hfun_nsim)
  
  probit_para = list(id_probit1 = id_probit1, id_probit0 = id_probit0)
  
  met_theta_para = list(metw=metw, metm=metm, met_alpha=met_alpha, mets=mets, metv=metv, met_mAM=met_mAM,
                        met_vAM=met_vAM, met_beta=met_beta, met_beta_AM=met_beta_AM, met_var_all=met_var_all,
                        iflagAM=iflagAM, pmet=pmet, icount_AM=icount_AM)
  
  met_gamma0_para = list(met_alpha_g0=met_alpha_g0, met_beta_g0=met_beta_g0)
  
  met_psi_para    = list(met_alpha_psi=met_alpha_psi, met_beta_psi=met_beta_psi)
  
  gamma_para = list(gmax=gmax, gamma_mu=gamma_mu, gamma_sigma=gamma_sigma,
                    gamma_beta=gamma_beta, gamma_alpha=gamma_alpha, gamma_prob=gamma_prob,
                    gamma_mu_m0=gamma_mu_m0, gamma_mu_v02=gamma_mu_v02, mu_bot=mu_bot, mu_top=mu_top,
                    gamma_sigma_m0=gamma_sigma_m0, gamma_sigma_v02=gamma_sigma_v02, sigma_bot=sigma_bot,
                    sigma_top=sigma_top, gam0vec=gam0vec, wk=wk)
  
  tau_para = list(tau2_s0=tau2_s0, tau2_rn=tau2_rn)
  eta_para = list(eta02_s0=eta02_s0, eta02_rn=eta02_rn, eta02=eta02, eta0=eta0, th0v=th0v)
  
  theta0_para = list(theta0_v02=theta0_v02)
  
  squish_para1 = c(zeta_mu, zeta_sigma2, zeta_sigma, zeta_m0, zeta_v02, zeta_s0, zeta_rn, zeta0, omega0, 
                   psi_fixed, psi_log_fixed, psi0, psi_log_mu, psi_log_sigma2, psi_log_sigma, 
                   psi_log_m0, psi_log_v02, psi_log_s0, psi_log_rn)
  
  squish_para2 = list(zetall=zetall, ezall=ezall, psiall=psiall, psiall_log=psiall_log, omegall=omegall, 
                      hfunall=hfunall)
  
  
  sigma_para = list(sigma2_alpha=sigma2_alpha, sigma2_beta=sigma2_beta, sigma2_alpha_m0=sigma2_alpha_m0,
                    sigma2_alpha_v02=sigma2_alpha_v02, sigma2_beta_s0=sigma2_beta_s0,
                    sigma2_beta_s0=sigma2_beta_s0, sigma2_beta_rn=sigma2_beta_rn)
  
  v_para = list(alpha_v0i=alpha_v0i, alpha_v0im0=alpha_v0im0)
  z_para = list(phi_v0i=phi_v0i, lambda_fn=lambda_fn, lambda_g0i=lambda_g0i)
  
  basis_para = list(kall=kall, ktheta=ktheta)
  phi_para = list(phixgrid=phixgrid, phi2xgrid=phi2xgrid, phixdata=phixdata)
  
  fx_para = list(fxgridall=fxgridall, f0xgrid=f0xgrid, f0Bxgrid=f0Bxgrid, fxobs=fxobs, fxobsm=fxobsm)              
  
  # mcmc start !!!
  stime=proc.time()
  mcmc_out = get_mcmc_hbsarv7(
    data_all=data_all, id_group=id_group-1, xdelta=xdelta, xgrid=xgrid, xrange=xrange, 
    xinx=xinx-1, xover=xover, idn0=idn0-1, nobs1=nobs1, iptHB1=iptHB2, xmin=xmin, ztz=ztz, 
    vtv=vtv, wtw=wtw, phi2=phi2, 
    cpara=cpara, probit_para=probit_para, gamma_para=gamma_para, tau_para=tau_para, eta_para=eta_para,
    met_theta_para=met_theta_para, met_gamma0_para=met_gamma0_para,
    met_psi_para=met_psi_para, 
    theta0_para=theta0_para, squish_para1=squish_para1, squish_para2=squish_para2, 
    sigma_para=sigma_para, v_para=v_para, z_para=z_para, phi_para=phi_para, fx_para=fx_para, 
    basis_para=basis_para,
    alpha=alpha, betall=betall, phi=phi, lambda=lambda, lambdai=lambdai, thetall=thetall, 
    theta0=theta0, xiparall=xiparall, xipar0=xipar0, gammall=gammall, tau2all=tau2all, tauall=tauall, 
    sigma2=sigma2, sigma=sigma, iflagpsi=iflagpsi, iflagCenter=iflagCenter, iflaglm=iflaglm,
    iflagpn=iflagpn, iflagsc=iflagsc, bFlagZ=iflagZ, disp=disp)
  mcmctime=proc.time()-stime
  
  pmet          = mcmc_out$pmet
  alphag        = mcmc_out$alphag
  sigmag        = mcmc_out$sigmag
  sigma2_alphag = mcmc_out$sigma2_alphag
  sigma2_betag  = mcmc_out$sigma2_betag
  tauallg       = mcmc_out$tauallg
  eta0g         = mcmc_out$eta0g
  gamma_mugg    = mcmc_out$gamma_mugg
  gamma_alphag  = mcmc_out$gamma_alphag
  gamma_betag   = mcmc_out$gamma_betag
  gammallg      = mcmc_out$gammallg
  theta0g       = mcmc_out$theta0g
  phig          = mcmc_out$phig
  lambdag       = mcmc_out$lambdag
  zetag         = mcmc_out$squishg$zetag
  omegag        = mcmc_out$squishg$omegag
  zeta0g        = mcmc_out$squishg$zeta0g
  omega0g       = mcmc_out$squishg$omega0g
  psiallg       = mcmc_out$squishg$psiallg
  psi0g         = mcmc_out$squishg$psi0g
  betag         = mcmc_out$betag
  fgridg        = mcmc_out$fgridg
  fxobsm        = mcmc_out$fxobsg[,1, drop=F]
  fxobss        = mcmc_out$fxobsg[,2, drop=F]
  thetam        = mcmc_out$thetam
  thetas        = mcmc_out$thetas
  
  fxgridm   = fgridg[,1:ngroup, drop=F]                        
  fxgrids   = fgridg[,(ngroup+1):(2*ngroup), drop=F]
  f0xgridm  = fgridg[,2*ngroup+1, drop=F]              
  f0xgrids  = fgridg[,2*ngroup+2, drop=F]           
  f0Bxgridm = fgridg[,2*ngroup+3, drop=F]             
  f0Bxgrids = fgridg[,2*ngroup+4, drop=F]              
  
  betam = betag[,1:nparw, drop=F]
  betas = betag[,-(1:nparw), drop=F]
  
  mcmc.draws = list()
  mcmc.draws$pmet          = pmet
  mcmc.draws$alphag        = alphag
  mcmc.draws$sigmag        = sigmag      
  mcmc.draws$sigma2_alphag = sigma2_alphag
  mcmc.draws$sigma2_betag  = sigma2_betag
  mcmc.draws$tauallg       = tauallg     
  mcmc.draws$eta0g         = eta0g       
  mcmc.draws$gamma_mugg    = gamma_mugg  
  mcmc.draws$gamma_alphag  = gamma_alphag
  mcmc.draws$gamma_betag   = gamma_betag 
  mcmc.draws$gammallg      = gammallg    
  mcmc.draws$theta0g       = theta0g     
  mcmc.draws$phig          = phig        
  mcmc.draws$lambdag       = lambdag     
  mcmc.draws$zetag         = zetag       
  mcmc.draws$omegag        = omegag      
  mcmc.draws$zeta0g        = zeta0g      
  mcmc.draws$omega0g       = omega0g     
  mcmc.draws$psiallg       = psiallg     
  mcmc.draws$psi0g         = psi0g       
  mcmc.draws$betag         = betag       
  mcmc.draws$fgridg        = fgridg      
  
  
  # Compute summary statistics
  pmetm       = pmet/nmcmc
  
  fxgridm     = fxgridm/smcmc
  fxgrids     = sqrt(abs(fxgrids - smcmc*fxgridm^2)/smcmc)
  f0xgridm    = f0xgridm/smcmc
  f0xgrids    = sqrt(abs(f0xgrids-smcmc*f0xgridm^2)/smcmc)
  f0Bxgridm   = f0Bxgridm/smcmc
  f0Bxgrids   = sqrt(abs(f0Bxgrids-smcmc*f0Bxgridm^2)/smcmc)
  fxobsm      = fxobsm/smcmc
  theta0m     = apply(theta0g,2,mean)
  theta0s     = apply(theta0g,2,sd)
  thetam      = thetam/smcmc
  thetas      = sqrt(abs(thetas-smcmc*thetam^2)/smcmc)
  eta0m       = apply(eta0g,2,mean)
  eta0s       = apply(eta0g,2,sd)
  
  
  fit.draws = list()
  fit.draws$fxobs    = fxobsm
  fit.draws$fxobss   = fxobss
  fit.draws$fxgrid   = fxgridm
  fit.draws$f0xgrid  = f0xgridm
  fit.draws$f0xgrids = f0xgrids
  
  post.est = list()
  post.est$fxgridm   = fxgridm
  post.est$fxgrids   = fxgrids
  post.est$f0xgridm  = f0xgridm
  post.est$f0xgrids  = f0xgrids
  post.est$f0Bxgridm = f0Bxgridm
  post.est$f0Bxgrids = f0Bxgrids
  post.est$fxobsm    = fxobsm
  post.est$theta0m   = theta0m
  post.est$theta0s   = theta0s
  post.est$thetam    = thetam
  post.est$thetas    = thetas
  post.est$eta0m     = eta0m
  post.est$eta0s     = eta0s
  
  if(nparv > 0){
    post.est$alpham    = matrix(apply(alphag,2,mean))
    post.est$alphas    = matrix(apply(alphag,2,sd))
    post.est$alpha_pg0 = matrix(apply(alphag>0,2,mean))
    post.est$valpham   = vdata %*% alpham
  }else{
    post.est$alpham    = 0
    post.est$alphas    = 1
    post.est$valpham   = matrix(0,ntot,1)
  }
  post.est$sigmam      = apply(sigmag,2,mean)
  post.est$sigmas      = apply(sigmag,2,sd)
  post.est$tauallm     = apply(tauallg,2,mean)
  post.est$taualls     = apply(tauallg,2,sd)
  
  post.est$gammallm      = apply(gammallg,2,mean) 
  post.est$gammalls      = apply(gammallg,2,sd) 
  post.est$gamma_mugm    = mean(gamma_mugg)
  post.est$gamma_mugs    = sd(gamma_mugg)
  post.est$gamma_alpham  = mean(gamma_alphag)
  post.est$gamma_alphas  = sd(gamma_alphag)
  post.est$gamma_betam   = mean(gamma_betag)
  post.est$gamma_betas   = sd(gamma_betag)
  
  post.est$phim       = matrix(apply(phig,2,mean),nparz,nparw)
  post.est$phis       = matrix(apply(phig,2,sd),  nparz,nparw)
  post.est$phi_pg0    = matrix(apply(phig>0,2,mean),nparz,nparw)
  post.est$lambdam    = matrix(apply(lambdag,2,mean),nparw,nparw)
  post.est$lambdas    = matrix(apply(lambdag,2,sd),  nparw,nparw)
  post.est$betam      = betam/smcmc
  post.est$betas      = sqrt(abs(betas - smcmc*betam^2)/smcmc)
  
  if(iflagsc==4){
    zetam    = apply(zetag,2,mean)
    zetas    = apply(zetag,2,sd)
    zeta0m    = mean(zeta0g)
    zeta0s    = sd(zeta0g)
    omegam    = apply(omegag,2,mean)
    omegas    = apply(omegag,2,sd)
    omega0m   = mean(omegag)
    omega0s   = sd(omegag)
    if(iflagpsi>0){
      psim    = apply(psiallg,2,mean)
      psis    = apply(psiallg,2,sd)
      psi0m   = mean(psi0g)
      psi0s   = sd(psi0g)
    }
  }
  
  if (nparw > 1) {
    wbetam = apply(wdata*post.est$betam[id_beta,],1,sum)
  } else {
    wbetam = wdata*(post.est$betam/smcmc)[id_beta]
  }
  post.est$wbetam = wbetam
  
  yresid = ydata - post.est$valpham - wbetam 
  
  res.out=list()
  res.out$model='HBSARv7'
  #res.out$family=family
  #res.out$link=link
  res.out$y=ydata
  res.out$yresid=yresid
  res.out$x=as.matrix(xdata)
  res.out$group=id_group.org
  if(!missing(wdata)) res.out$w=wdata
  if(!missing(zdata)) res.out$z=zdata
  if(!missing(vdata)) res.out$v=vdata
  res.out$n      = ntot
  res.out$ngroup = ngroup
  res.out$nparw  = nparw
  res.out$nparv  = nparv
  res.out$nparz  = nparz
  res.out$nint   = nint
  res.out$nbasis = nbasis
  res.out$xgrid  = xgrid    
  res.out$yname  = yname
  res.out$zname  = zname
  res.out$wnames = wname
  res.out$xname  = xname
  
  res.out$shape=shape
  
  res.out$prior=prior
  
  res.out$mcmctime=mcmctime
  res.out$mcmc=mcvals
  res.out$pmet=pmet
  #res.out$imodmet=imodmetg
  
  res.out$mcmc.draws = mcmc.draws
  res.out$fit.draws  = fit.draws
  res.out$post.est   = post.est
  return(res.out)
}
