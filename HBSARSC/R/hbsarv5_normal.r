"hbsarv5" <- function(y, w, x, v, z, nbasis,id_group, 
                      nint=200, mcmc = list(), prior = list(),xmin,xmax, 
                      shape = "Increasing", iflagCenter=FALSE, disp=FALSE) {
  
  # check arguments
  if (missing(y)) stop("y is required")
  if (missing(x)) stop("x is required")
  
  if (missing(w)) w = NULL
  if (missing(v)) v = NULL
  if (missing(z)) z = NULL
  
  
  shape_type = c("Increasing", "Decreasing")
  shape <- match.arg(shape,shape_type,several.ok=TRUE)
  
  
  ydata = y
  wdata = w
  vdata = v
  zdata = z
  xdata = x
  
  # type of shape restriction
  fshape=function.shape(shape)
  fmodel=fshape$fmodel
  fpm   =fshape$fpm
  nfun  =fshape$nfun
  
  # response
  if(!is.matrix(ydata)) ydata=as.matrix(ydata)
  yname=colnames(ydata)
  if(is.null(yname)) yname='y'
  colnames(ydata)=yname
  ntot=length(ydata)
  
  #group
  if(missing(id_group))id_group=rep(1,ntot)
  
  # group infromation
  ugroups  = sort(unique(id_group))  # unique groups. Note: unique does not sort
  ngroup   = length(ugroups)
  if(ngroup != length(ugroups)){
    cat("\nError with group ids in Zdata and 2nd column of datall\n")
  }
  
  # W: design matrix for random parametric effect
  # W includes intercept
  if(is.null(wdata)) {
    wdata=matrix(1.0,nrow=ntot,ncol=1)
    wnames='CNST'
    colnames(wdata)=wnames
    ndimw=0 # number of W variables, excluding the intercept
  }else{
    if(!is.matrix(wdata)) wdata=as.matrix(wdata)
    wnames=colnames(wdata)
    if(is.null(wnames)) wnames=paste('w',1:ncol(wdata),sep='')
    wdata=cbind(1.0,wdata)
    wnames=c('CNST',wnames)
    colnames(wdata)=wnames
    ndimw=ncol(wdata)-1
  }
  
  nparw=ndimw+1
  
  #-----------------------------------
  # ver5. has a fixed effect V
  # V does not include intercept term 
  #-----------------------------------
  
  if(is.null(vdata)){
    vdata  = matrix(0, ntot,1)
    vnames = "None"
    vtv    = matrix(1, 1,1)
    nparv  = 0
  } else{
    if(!is.matrix(vdata)) vdata  = as.matrix(vdata)
    vnames = colnames(vdata)
    if(is.null(vnames)) vnames=paste('v',1:ncol(vdata),sep="")
    vtv    = crossprod(vdata)
    nparv=ncol(vdata)
  }
  
  #-------------------------------------------------------------
  # ver5. has Z matrix for explaning parametric coefficient beta
  # beta_j = z_j*Phi + delta, delta ~ N(0,Lambda)
  # Z matrix includes intercept term
  #------------------------------------------------------------ 
  
  if(is.null(zdata)) {
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
  nparz=ncol(zdata)
  
  # nonparametric part
  if(!is.matrix(xdata)) xdata=as.matrix(xdata)
  if(nfun!=ncol(xdata)) stop('The number of shape and columns of x should be same.')
  xname=colnames(xdata)
  if(is.null(xname)) xname=paste('x',1:nfun,sep='')
  colnames(xdata)=xname
  
  #--------------------------------------------------------------
  # * MCMC parameters
  #--------------------------------------------------------------
  mcvals=list(nblow0=1000,nblow=10000,smcmc=1000,nskip=10,ndisp=1000,maxmodmet=5)
  mcvals[names(mcmc)]=mcmc
  nblow0    = mcvals$nblow0
  nblow     = mcvals$nblow
  smcmc     = mcvals$smcmc
  nskip     = mcvals$nskip
  ndisp     = mcvals$ndisp
  maxmodmet = mcvals$maxmodmet
  
  nmcmc     = nblow + nskip*smcmc
  
  nmcmcall  = maxmodmet*nblow0 + nmcmc
  
  # id_group may not be sequential from 1 ... ngroup
  # Need an index to get correct betam
  id_beta  = matrix(0,ntot,1)
  for(j in 1:ngroup){
    id = which(id_group == ugroups[j])
    id_beta[id] = matrix(j,length(id),1)
  }
  
  outpt  = MakeHBPointer(id_group)
  nobs1  = outpt$nobs
  iptHB1 = outpt$iptd
  
  iptHB2 = vector("list", length = length(iptHB1)) # index for c++
  for (i in 1:length(iptHB1)) {
    iptHB2[[i]] = iptHB1[[i]]-1  
  }
  
  # Pointer into w'w
  iptw2   = as.matrix(cumsum(rep(nparw,ngroup)))
  iptw2   = cbind(c(1,(iptw2[1:(ngroup-1),]+1)),iptw2)
  
  # Compute w'w for each population
  wtw     = matrix(0,nparw*ngroup,nparw)
  for(j in 1:ngroup){
    wj    = wdata[iptHB1[[j]],]
    wj    = as.matrix(wj)
    wtw[iptw2[j,1]:iptw2[j,2],] = crossprod(wj)
  }
  
  # Compute z'z 
  ztz = crossprod(zdata)
  
  # Grid on [xmim,xmax] for plotting functions
  # Also, grid for integration
  if (missing(xmin)) xmin = min(xdata); 
  if (missing(xmax)) xmax = max(xdata)
  a      = floor(xmin)
  if(a<xmin) xmin = a
  b      = ceiling(xmax)
  if(b>xmax) xmax = b
  xdelta = (xmax-xmin)/nint
  xgrid  = seq(xmin,xmax,xdelta)
  xgrid  = matrix(xgrid,nint+1,1)
  xrange = xmax-xmin
  
  #------------------------------------------------
  # Computational strategy
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
  met_alpha = rep(3,ngroup);        #  Keep alpha fixed
  
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
  metm_g0      = 0.00001               # mean of metropolis step size
  mets_g0      = 0.01                  # std dev of metroplis step size
  met_alpha_g0 = 2 + (metm_g0/mets_g0)^2
  met_beta_g0  = metm_g0*(met_alpha_g0 -1)
  
  # Initialize smoothing parameters
  kall     = matrix(1:nbasis,nbasis,1)   # number of cosine functions, including intercept
  kall0    = matrix(0:nbasis,nbasis+1,1)
  
  
  # Compute Cosine functions on xgrid and xdata 
  phi0xgrid = matrix(1/sqrt(xrange),nint+1,1)
  phixgrid  = sqrt(2/xrange)*cos(pi*xgrid%*%(t(kall)-xmin)/xrange)
  phixgrid  = cbind(phi0xgrid,phixgrid)
  phi2xgrid = phixgrid^2    # Cosine squared functions on xgrid.  Used for bias
  
  
  #-----------------------------------------------------
  # Initialize priors
  privals=list(iflagprior=0,tau2_r0=5,tau2_s0=10,w0=1,#alpha_m0=3, alpha_s0=50, 
               iflagpsi=1,psifixed=100,omega_m0=(xmin+xmax)/2,omega_s0=xrange/8, # for shape restriction 
               alpha_m0=0,alpha_v0=1000,phi_b0=0,lambda_f0=5, lambda_g0=10,eta_r0=5,eta_s0=10,gamma_lmu_m0=0,gamma_lmu_v0=100,gamma_lsigma_m0=0,gamma_lsigma_v0=100,# from v5
               sigma2_r0=5,sigma2_s0=10,   # normal
               kappa_m0=1,kappa_v0=100,    # negative binomial 
               varphi_m0=1,varphi_v0=100)  # varphi for beta family

  privals[names(prior)]=prior
  
  
  # Fixed effect alpha ~ N(m0,V0)
  if(nparv>0){
    alpha_m0    = matrix(privals$alpha_m0,nparv,1)
    alpha_v0    = privals$alpha_v0*diag(nparv)
    alpha_v0i   = diag(nparv)/1000
    alpha_v0im0 = alpha_v0i %*% alpha_m0
  }else{  # no fixed effects: use placeholders
    alpha_m0    = matrix(0,      1,1)
    alpha_v0    = matrix(1000,   1,1)
    alpha_v0i   = matrix(1/1000, 1,1)
    alpha_v0im0 = matrix(0,      1,1)
  }
  
  # When Normal
  # Error variance sigma2 ~ IG(r0/2,s0/2)
  sigma2_r0  = privals$sigma2_r0
  sigma2_s0  = privals$sigma2_s0
  sigma2_rn  = sigma2_r0 + ntot
  
  # HB Regression beta = zdata*Phi + delta
  phidim    = nparw*nparz
  phi_b0    = matrix(privals$phi_b0,phidim,1)    # prior mean 
  phi_v0i   = diag(phidim)/10000                 # prior precision
  phi_v0ib0 = matrix(0,phidim,1)
  
  # Var(delta) = Lambda ~ IW
  lambda_f0  = privals$lambda_f0
  lambda_fn  = lambda_f0 + ngroup
  lambda_g0  = privals$lambda_g0*diag(nparw)
  lambda_g0i = diag(nparw)/10
  
  # theta_jk ~ N(theta_0k, tau_k^2*exp(-k*gamma))
  # tau_k^2  ~ IG(r0/2,s0/2)
  tau2_r0    = privals$tau2_r0
  tau2_s0    = privals$tau2_s0
  tau2_rn    = tau2_r0 + ngroup
  
  #  theta_0k ~ N(0,eta*exp(-k*gamma))
  #  eta ~
  eta_r0    = privals$eta_r0
  eta_s0    = privals$eta_s0
  eta_rn    = eta_r0 + nbasis+1
  
  # gamma_j ~ Gamma(alpha,beta) and gamma_0 = alpha/beta
  # Changer parameterization
  # mu = alpha/beta and sigma = sqrt(alpha)/beta
  # ln(mu) and ln(sigma) have normal distributions
  gamma_lmu_m0    = privals$gamma_lmu_m0
  gamma_lmu_v0    = privals$gamma_lmu_v0
  gamma_lsigma_m0 = privals$gamma_lsigma_m0
  gamma_lsigma_v0 = privals$gamma_lsigma_v0
  # gamma_j ~ G(alpha,beta).  Set alpha > abot
  abot             = 2   # Keep alpha > abot
  labot            = log(sqrt(abot))
  
  w0   = privals$w0
  wk   = sum(kall)/2;
  
  #------------------------------------------------------
  # Initialize parameters
  alpha           = alpha_m0
  betall          = matrix(0,ngroup,nparw)
  phi             = matrix(0,nparz,nparw)
  phidim          = nparw*nparz
  lambda          = diag(1,nparw,nparw)
  lambdai         = diag(1,nparw,nparw)
  thetall         = matrix(0,nbasis+1,ngroup)
  xiparll         = matrix(-.5,ngroup,1)                # xipar_j0 = log(theta_j0)
  thetall[1,]     = matrix(exp(xiparll),1,ngroup)       # lower level spectral coefficients
  
  theta0          = matrix(0,nbasis+1,1)     # upper level spectral coefficients
  xipar0          = 0                        # xipar_00 = log(theta_00)
  theta0[1,]      = 1
  
  gammall         = matrix(1,ngroup,1)       # lower level smoothing
  
  # Initialize parameters for gamma_j ~ G(alpha,beta)
  gamma_alpha     = 3            # Bound gamma_alpha > abot = 2, say
  gamma_beta      = gamma_alpha
  gamma_mu        = gamma_alpha/gamma_beta
  gamma_sigma2    = gamma_mu/gamma_beta
  gamma_sigma     = sqrt(gamma_sigma2)
  gamma_lmu       = log(gamma_mu)
  gamma_lsigma    = log(gamma_sigma)
  gamma0          = gamma_mu                                       # upper level smoothing
  
  sigma2          = 1
  sigma           = sqrt(sigma2)
  tau2all         = matrix(1,nbasis+1,1)   # upper level variances for theta_jk
  tauall          = sqrt(tau2all)
  eta2            = 1                      # variance parameter for theta_0k, k > 0
  eta             = sqrt(eta2)
  
  
  wbeta           = matrix(0,ntot,1)
  for(j in 1:ngroup){   # Compute w*beta
    wj    = as.matrix(wdata[iptHB1[[j]],])
    bj    = as.matrix((betall[j,]))
    wbeta[iptHB1[[j]],1] = wj %*% bj
  }
  
  if(nparv>0){
    valpha = vdata %*% alpha
  } else {
    valpha = matrix(0,ntot,1)
  }
  
  #------------------------------------------------------
  # Matrices for saving MCMC iterations
  if(nparv>0){
    alphag  = matrix(0,smcmc,nparv)
  }else{
    alphag  = matrix(0,smcmc,1)
  }
  
  sigmag    = matrix(0,smcmc,1)           # Error variance
  tauallg   = matrix(0,smcmc,nbasis+1)    # HB stdev for spectral coefficients
  etag      = matrix(0,smcmc,1)           # STDEV for upper level spectral coefficients
  gamma0g   = matrix(0,smcmc,1)           # Smoothing parameter for upper level model
  gammallg  = matrix(0,smcmc,ngroup)      # Smoothing parameters for lower level model
  thetam    = matrix(0,nbasis+1,ngroup)  
  thetas    = thetam
  theta0g   = matrix(0,smcmc,nbasis+1)
  fxgridm   = matrix(0,nint+1,ngroup)
  fxgrids   = matrix(0,nint+1,ngroup)
  dfxgridm  = matrix(0,nint+1,ngroup)
  f0xgridm  = matrix(0,nint+1,1)
  f0xgrids  = matrix(0,nint+1,1)
  
  # Added 19.11.25
  f0Bxgridm  = matrix(0,nint+1,1)
  f0Bxgrids  = matrix(0,nint+1,1)
  
  phig      = matrix(0,smcmc,phidim)
  nparw2    = nparw*nparw
  lambdag   = matrix(0,smcmc,nparw2)
  betam     = betall
  betas     = betall
  
  gamma_alphag  = matrix(0,smcmc,1)
  gamma_betag   = matrix(0,smcmc,1)
  gamma_lmug    = matrix(0,smcmc,1)
  gamma_lsigmag = matrix(0,smcmc,1)
  
  outfx           = rcpp_GetUpfxgrid(thetall[,1],phixgrid,xdelta,xrange,iflagCenter)
  
  fxgridj         = outfx$fx      # f computed on xgrid
  dfxgridj        = outfx$dfx     # derivative of f computed on xgrid
  fxgridall       = matrix(fxgridj,nint+1,ngroup)
  dfxgridall      = matrix(dfxgridj,nint+1,ngroup)
  f0xgrid         = matrix(0,nint+1,1)
  f0Bxgrid        = f0xgrid
  fxobs           = matrix(0,ntot,1)
 
  # No metropolis for theta if free f. Set maxmodemet=0
  if(max(fmodel)==0) maxmodmet=0
  
  data_all = cbind(ydata, valpha, wbeta, fxobs)
  cpara = c(ntot, ngroup, nparv, nparw, nparw2, nparz, phidim, nbasis, nblow, nblow0, maxmodmet, nskip,
            nmcmcall, smcmc, nint)
  met_vec_para = cbind(metw, metm, met_alpha, mets, metv, met_mAM, met_vAM, met_beta, met_beta_AM, met_var_all)
  met_para = c(iflagAM, met_alpha_g0, met_beta_g0)
  gamma_para =  c(gamma_alpha,gamma_beta,gamma_sigma,gamma_lmu,gamma_lsigma,
                  gamma0,gamma_lmu_m0,gamma_lmu_v0,gamma_lsigma_m0,gamma_lsigma_v0)
  sigma_para = c(sigma2_s0, sigma2_rn, sigma2, sigma)
  tau_para = c(tau2_s0, tau2_rn)
  eta_para = c(eta_s0, eta_rn, eta2, eta)
  
  fndfxgridall = cbind(fxgridall, dfxgridall)
  f0 = cbind(f0xgrid, f0Bxgrid)
  
  # mcmc start !!!
  stime=proc.time()
  mcmc_out = get_mcmc_hbsarv5_nromal( 
    data_all=data_all, vdata=vdata, wdata=wdata, zdata=zdata, id_group=id_group-1, iptw2=iptw2-1, 
    wtw=wtw, ztz=ztz, xdelta=xdelta, xgrid=xgrid, xrange=xrange, xinx=xinx-1, xover=xover, idn0=idn0-1,
    met_vec_para=met_vec_para, met_para=met_para, phi0xgrid=phi0xgrid, phixgrid=phixgrid, 
    phi2xgrid=phi2xgrid, fndfxgridall=fndfxgridall, f0=f0, vtv=vtv, betall=betall, phi=phi, 
    lambda=lambda, lambdai=lambdai, thetall=thetall, xiparll=xiparll, theta0=theta0, xipar0=xipar0, 
    gammall=gammall, gamma_para=gamma_para, tau2all=tau2all, tauall=tauall, alpha_m0=alpha_m0, 
    alpha_v0=alpha_v0, alpha_v0i=alpha_v0i, alpha_v0im0=alpha_v0im0, sigma_para=sigma_para, 
    phi_b0=phi_b0, phi_v0i=phi_v0i, phi_v0ib0=phi_v0ib0, lambda_f0=lambda_f0, lambda_fn=lambda_fn, 
    lambda_g0=lambda_g0, lambda_g0i=lambda_g0i, tau_para=tau_para, eta_para=eta_para, abot=abot, 
    labot=labot, nobs1=nobs1, iptHB1=iptHB2, cpara=cpara, iflagCenter=iflagCenter, iflagpn=fpm, disp=disp)
  mcmctime=proc.time()-stime
  
  pmetg   = mcmc_out$pmet   
  alphag  = mcmc_out$alphag 
  sigmag  = mcmc_out$sigmag         
  tauallg = mcmc_out$tauallg
  etag    = mcmc_out$etag   
  gammag  = mcmc_out$gammag 
  thetam  = mcmc_out$thetam 
  thetas  = mcmc_out$thetas 
  theta0g = mcmc_out$theta0g
  fgridg  = mcmc_out$fgridg 
  fxobsm  = mcmc_out$fxobsg[,1]
  fxobss  = mcmc_out$fxobsg[,2]
  phig    = mcmc_out$phig   
  lambdag = mcmc_out$lambdag
  betag   = mcmc_out$betag
  
  gamma0g       = gammag[,1]
  gamma_alphag  = gammag[,2]
  gamma_betag   = gammag[,3]
  gamma_lmug    = gammag[,4]
  gamma_mug     = gammag[,5]
  gamma_lsigmag = gammag[,6]
  gammallg      = gammag[,-(1:6)]
  
  fxgridm   = fgridg[,1:ngroup]                        
  fxgrids   = fgridg[,(ngroup+1):(2*ngroup)]  
  dfxgridm  = fgridg[,(2*ngroup+1):(3*ngroup)]
  f0xgridm  = fgridg[,3*ngroup+1]              
  f0xgrids  = fgridg[,3*ngroup+2]           
  f0Bxgridm = fgridg[,3*ngroup+3]             
  f0Bxgrids = fgridg[,3*ngroup+4]              
  
  betam = betag[,1:nparw]
  betas = betag[,-(1:nparw)]
  
  mcmc.draws = list()
  if (nparv>0) mcmc.draws$alpha = alphag
  mcmc.draws$tauall         = tauallg
  mcmc.draws$eta            = etag
  mcmc.draws$gamma0         = gamma0g
  mcmc.draws$gamma_alpha    = gamma_alphag
  mcmc.draws$gamma_beta     = gamma_betag
  mcmc.draws$gamma_lmu      = gamma_lmug
  mcmc.draws$gamma_lsigma   = gamma_lsigmag
  mcmc.draws$gammall        = gammallg
  mcmc.draws$theta0         = theta0g
  mcmc.draws$phi            = phig
  mcmc.draws$lambda         = lambdag
  mcmc.draws$beta           = betam
  mcmc.draws$theta          = thetam
  mcmc.draws$thetas         = thetas
  mcmc.draws$phi            = phig
  mcmc.draws$lambda         = lambdag
  mcmc.draws$sigma          = sigmag
  
  fit.draws = list()
  fit.draws$fxobs    = fxobsm
  fit.draws$fxobss   = fxobss
  fit.draws$dfxgrid  = dfxgridm
  fit.draws$fxgrid   = fxgridm
  fit.draws$f0xgrid  = f0xgridm
  fit.draws$f0xgrids = f0xgrids
  
  post.est = list()
  post.est$pmetm       = pmet/nmcmc
  post.est$fxgridm     = fxgridm/smcmc
  post.est$fxgrids     = sqrt(abs(fxgrids - smcmc*fxgridm^2)/smcmc)
  post.est$dfxgridm    = dfxgridm/smcmc
  post.est$f0xgridm    = f0xgridm/smcmc
  post.est$f0xgrids    = sqrt(abs(f0xgrids-smcmc*f0xgridm^2)/smcmc)
  post.est$f0Bxgridm   = f0Bxgridm/smcmc
  post.est$f0Bxgrids   = sqrt(abs(f0Bxgrids-smcmc*f0Bxgridm^2)/smcmc)
  post.est$theta0m     = apply(theta0g,2,mean)
  post.est$theta0s     = apply(theta0g,2,sd)
  post.est$thetam      = thetam/smcmc
  post.est$thetamm     = apply(post.est$thetam,1,mean)
  post.est$thetas      = sqrt(abs(thetas-smcmc*thetam^2)/smcmc)
  
  
  if(nparv > 0){
    alpham    = matrix(apply(alphag,2,mean))
    alphas    = matrix(apply(alphag,2,sd))
    alpha_pg0 = matrix(apply(alphag>0,2,mean))
    valpham   = vdata %*% alpham
  }else{
    alpham    = 0
    alphas    = 1
    valpham   = matrix(0,ntot,1)
  }
  
  if (nparw > 1) {
    wbetam     = apply(wdata*(betam/smcmc)[id_beta,],1,sum)
  } else {
    wbetam     = wdata*(betam/smcmc)[id_beta]
  }
  
  post.est$alpham  = alpham
  post.est$alphas  = alphas
  post.est$valpham = valpham
  post.est$sigmam  = mean(sigmag)
  post.est$sigmas  = sd(sigmag)
  
  post.est$tauallm       = apply(tauallg,2,mean)
  post.est$taualls       = apply(tauallg,2,sd)
  post.est$etam          = mean(etag)
  post.est$etas          = sd(etag)
  post.est$gammallm      = apply(gammallg,2,mean) 
  post.est$gammalls      = apply(gammallg,2,sd) 
  post.est$gamma0m       = mean(gamma0g)
  post.est$gamma0s       = sd(gamma0g)
  post.est$gamma_alpham  = mean(gamma_alphag)
  post.est$gamma_alphas  = sd(gamma_alphag)
  post.est$gamma_betam   = mean(gamma_betag)
  post.est$gamma_betas   = sd(gamma_betag)
  post.est$phim       = matrix(apply(phig,2,mean),    nparz,nparw)
  post.est$phis       = matrix(apply(phig,2,sd),      nparz,nparw)
  post.est$phi_pg0    = matrix(apply(phig>0,2,mean) , nparz,nparw)
  post.est$lambdam    = matrix(apply(lambdag,2,mean), nparw,nparw)
  post.est$lambdas    = matrix(apply(lambdag,2,sd),   nparw,nparw)
  post.est$betam      = betam/smcmc
  post.est$betas      = sqrt(abs(betas - smcmc*betam^2)/smcmc)
  post.est$wbetam     = wbetam
  
  res.out=list()
  res.out$model='HBSARv5'
  #res.out$family=family
  #res.out$link=link
  res.out$y=ydata
  res.out$yresid=ydata - valpham - wbetam     
  res.out$x=as.matrix(xdata)
  res.out$group=id_group
  if(!missing(wdata)) res.out$w=wdata
  if(!missing(zdata)) res.out$z=zdata
  if(!missing(vdata)) res.out$v=vdata
  res.out$n=ntot
  res.out$ngroup= ngroup
  res.out$nparw = nparw
  res.out$nparv = nparv
  res.out$nparz = nparz
  res.out$nint  = nint
  res.out$nbasis= nbasis
  res.out$xgrid = xgrid    
  res.out$yname = yname
  res.out$wnames= wnames
  res.out$xname = xname
  
  res.out$shape=shape
  
  res.out$prior=prior
  
  res.out$mcmctime=mcmctime
  res.out$mcmc=mcvals
  res.out$pmet=pmetg
  #res.out$imodmet=imodmetg
  
  res.out$mcmc.draws = mcmc.draws
  res.out$fit.draws  = fit.draws
  res.out$post.est   = post.est
  
  return(res.out)
}
