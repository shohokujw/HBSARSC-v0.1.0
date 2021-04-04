"gHBSARSC" <- function(shape="Increasing", ntot=2000, ngroup=100, nbasis=20, nparv=0, nparw=2, nparz=3, nparx=0,
                       xmin=0, xmax=10, nint=500, iflagZ=0, iflagHBsigma=0, iflagpsi=1, iflaglm=0, 
                       iflagCenter=1, iflagSpanX=0, iflagRandn=0, disp=F) {
  
  if (shape == "Free") {
    iflagsc = 0
    iflagpn = 1
  } else if (shape == "Increasing") {
    iflagsc = 1
    iflagpn = 1
  } else if (shape == "Decreasing") {
    iflagsc = 1
    iflagpn = -1
  } else if (shape == "IncreasingConvex") { # increasing and convex or decreasing and concave
    iflagsc = 2
    iflagpn = 1
  } else if (shape == "DecreasingConcave") {
    iflagsc = 2
    iflagpn = -1
  } else if (shape == "IncreasingConcave") { # increasing and concave or decreasing and convex
    iflagsc = 3
    iflagpn = 1
  } else if (shape == "DecreasingConvex") {
    iflagsc = 3
    iflagpn = -1
  } else if (shape == "InvertedU") {
    iflagsc  = 4
    iflagpn = 1
  } else if (shape == "Ushape") {
    iflagsc = 4
    iflagpn = -1
  } else if (shape == "IncreasingS") { # increasing convex-to-concave or decreasing concave-to-convex
    iflagsc = 5
    iflagpn = 1
  } else if (shape == "DecreasingS") { 
    iflagsc = 5
    iflagpn = -1
  } else if (shape == "IncreasingRotatedS") { # increasing concave-to-convex or decreasing convex-to-concave 
    iflagsc = 6
    iflagpn = 1
  } else if (shape == "DecreasingRotatedS") {
    iflagsc = 6
    iflagpn = -1
  }  else {
    stop("Not available")
  }
  
  
  ## local functions
  MakeHBPointer <- function(id){
    uid = sort(unique(id))   # unique populations
    nf  = length(uid)        # number of unique groups
    
    nobs = matrix(0,nf,1)          # number of observations in each group
    for(j in 1:nf){                # loop over groups
      jrow = which(id == uid[j])   # identify rows for group j
      nobs[j] = length(jrow)       # number of observations in group j
      if(j == 1){
        iptd = list(jrow)          # start the list iptd
      }else{
        iptd[[j]] = jrow           # add new member to list
      }
    }  # end loop over groups
    return(list("nobs"=nobs,"iptd"=iptd))
  }
  
  trapint <- function(f,xdelta){
    nr = length(f)
    a  = sum(f)
    a  = a - (f[1] + f[nr])/2
    return(a*xdelta)
  }
  
  CumTrap = function(f,delta){
    nr    = NROW(f)
    fcum	= cumsum(f[1:nr-1]) + cumsum(f[2:nr])
    fcum	= fcum*delta/2
    fcum	= matrix(c(0,fcum),nr,1)
    return(fcum)
  }
  
  GetSquish = function(omega,psi,xgrid){
    c       = psi*(xgrid-omega)
    # Worry about underflow and over flow
    # if c < -100, then hfun = 1
    # if c >  100, then hfun = -1
    idz     = which(c< -100)
    if(length(idz>0)) c[idz] = -100
    idz     = which(c>100)
    if(length(idz>0)) c[idz] = 100
    c       = exp(c)
    return((1-c)/(1+c))
  }
  
  GetMonofxgrid  = function(fpars){
    z     = fpars$phix %*% fpars$theta
    if(iflagZ == 1){
      dfx  = exp(z)
    }else{
      dfx  = z^2
    }
    
    fx         = fpars$iflagpn*CumTrap(dfx,fpars$delta)
    if(fpars$iflagCenter==0){  # do not center f
      return(fx)
    }else{  # center f
      # Compute mean
      c         = trapint(fx,fpars$delta)
      return(fx-c/fpars$range)
    }
  }  # End GetMonofxgrid
  
  GetMonoSamefxgrid  = function(fpars){
    z     = fpars$phix %*% fpars$theta
    if(iflagZ == 1){
      d2fx  = exp(z)
    }else{
      d2fx  = z^2
    }
    
    dfx         = fpars$iflagpn*CumTrap(d2fx,fpars$delta)
    fx          = CumTrap(dfx,fpars$delta)
    if(fpars$iflagCenter==0){  # do not center f
      return(fx)
    }else{  # center f
      # Compute mean
      c         = trapint(fx,fpars$delta)
      return(fx-c/fpars$range)
    }
  }  # End GetMonoSamefxgrid
  
  GetMonoDiffxgrid  = function(fpars){
    z     = fpars$phix %*% fpars$theta
    if(iflagZ == 1){
      d2fx  = exp(z)
    }else{
      d2fx  = z^2
    }
    dfx         = fpars$iflagpn*CumTrap(d2fx,fpars$delta)
    fx          = CumTrap(rev(dfx),fpars$delta)
    if(fpars$iflagCenter==0){  # do not center f
      return(fx)
    }else{  # center f
      # Compute mean
      c         = trapint(fx,fpars$delta)
      return(fx-c/fpars$range)
    }
  }  # End GetMonoDiffxgrid
  
  GetUfxgrid  = function(fpars){
    fx         = 0
    z         = fpars$phix %*% fpars$theta
    if(iflagZ == 1){
      dfx  = exp(z)
    }else{
      dfx  = z^2
    }
    dfx   = fpars$hfun*dfx
    fx         = fpars$iflagpn*CumTrap(dfx,fpars$delta)
    if(fpars$iflagCenter==0){  # do not center f
      return(fx)
    }else{  # center f
      # Compute mean
      c         = trapint(fx,fpars$delta)
      return(fx-c/fpars$range)
    }
  }  # End GetUfxgrid
  
  GetSfxgrid  = function(fpars,d3){
    z     = fpars$phix %*% fpars$theta
    if(iflagZ == 1){
      d2fx  = exp(z)
    }else{
      d2fx  = z^2
    }
    d2fx        = d3*fpars$hfun*d2fx
    # d3 =  1 -> hfun is "squish down" for "S"
    # d3 = -1 -> hfun is "squish up" for "rotated S"
    dfx         = CumTrap(d2fx,fpars$delta)
    c2          = min(c(0,dfx))
    fx          = CumTrap(dfx,fpars$delta)
    fx          = fpars$iflagpn*(fx - c2*(xgrid-xmin))
    if(fpars$iflagCenter==0){  # do not center f
      return(fx)
    }else{  # center f
      # Compute mean
      c         = trapint(fx,fpars$delta)
      fx        = fx - c/fpars$range
      return(fx)
    }
  }  # End GetSfxgrid
  
  Getfx = function(fpars){
    fx     = 0
    # free function
    if(fpars$iflagsc == 0)  return(fpars$phix %*% fpars$theta)
    #------------------------------------------
    # monotone function
    if(fpars$iflagsc == 1){
      fx = GetMonofxgrid(fpars)
      return(fx)
    } 
    
    # Increasing convex or decreasing and concave
    if(fpars$iflagsc == 2){
      fx = GetMonoSamefxgrid(fpars)
      return(fx)
    }  
    # Increasing concave or decreasing and convex
    if(fpars$iflagsc == 3){
      fx = GetMonoDiffxgrid(fpars)
      return(fx)
    }  
    #-------------------------------------------
    # U-shaped function
    if(fpars$iflagsc == 4){
      fx = GetUfxgrid(fpars)
      return(fx)
    }  
    # S-shaped (increasing,convex-to-concave) or (decreasing,concave-to-convex)
    if(fpars$iflagsc == 5){
      fx  = GetSfxgrid(fpars,1)
      return(fx)
    }
    # S-shaped (increasing,concave-to-convex) or (decreasing,convex-to-concave)
    if(fpars$iflagsc == 6){
      fx  = GetSfxgrid(fpars,-1)
      return(fx)
    }
    
  }  # End Getfx
  
  GetMonof0xgrid  = function(f0pars){
    
    z          = f0pars$phix %*% f0pars$theta
    v2         = f0pars$phix2 %*% f0pars$vpara
    if(iflagZ == 1){
      dfx  = exp(z + v2/2)
    }else{
      dfx  = z^2 + v2
    }
    fx         = f0pars$iflagpn*CumTrap(dfx,f0pars$delta)
    if(iflagCenter==0){  # do not center f
      return(fx)
    }else{  # center f
      # Compute mean
      c         = trapint(fx,f0pars$delta)
      return(fx-c/f0pars$range)
    }
  }  # End GetMonof0xgrid
  
  GetUf0xgrid  = function(f0pars){
    z          = f0pars$phix %*% f0pars$theta
    v2         = f0pars$phix2 %*% f0pars$vpara
    if(iflagZ == 1){
      dfx  = exp(z + v2/2)
    }else{
      dfx  = z^2 + v2
    }
    dfx        = f0pars$hfun*dfx
    fx         = f0pars$iflagpn*CumTrap(dfx,f0pars$delta)
    if(f0pars$iflagCenter==0){  # do not center f
      return(fx)
    }else{  # center f
      # Compute mean
      c         = trapint(fx,f0pars$delta)
      return(fx-c/f0pars$range)
    }
  }  # End GetUf0xgrid
  
  GetMonoSamef0xgrid  = function(f0pars){
    
    z          = f0pars$phix %*% f0pars$theta
    v2         = f0pars$phix2 %*% f0pars$vpara
    if(iflagZ == 1){
      d2fx  = exp(z + v2/2)
    }else{
      d2fx  = z^2 + v2
    }
    dfx         = f0pars$iflagpn*CumTrap(d2fx,f0pars$delta)
    fx          = CumTrap(dfx,f0pars$delta)
    if(iflagCenter==0){  # do not center f
      return(fx)
    }else{  # center f
      # Compute mean
      c         = trapint(fx,f0pars$delta)
      return(fx-c/f0pars$range)
    }
  }  # End GetMonoSamef0xgrid
  
  GetMonoDiff0xgrid  = function(f0pars){
    
    z          = f0pars$phix %*% f0pars$theta
    v2         = f0pars$phix2 %*% f0pars$vpara
    if(iflagZ == 1){
      d2fx  = exp(z + v2/2)
    }else{
      d2fx  = z^2 + v2
    }
    dfx         = f0pars$iflagpn*CumTrap(d2fx,f0pars$delta)
    fx          = CumTrap(rev(dfx),f0pars$delta)
    if(iflagCenter==0){  # do not center f
      return(fx)
    }else{  # center f
      # Compute mean
      c         = trapint(fx,f0pars$delta)
      return(fx-c/f0pars$range)
    }
  }  # End GetMonoSamef0xgrid
  
  GetSf0xgrid  = function(fpars,d3){
    z          = fpars$phix  %*% fpars$theta
    v2         = fpars$phix2 %*% fpars$vpara
    z     = fpars$phix %*% fpars$theta
    if(iflagZ == 1){
      d2fx  = exp(z + v2/2)
    }else{
      d2fx  = z^2 + v2
    }
    d2fx        = d3*fpars$hfun*d2fx
    # d3 =  1 -> hfun is "squish down"
    # d3 = -1 -> hfun is "squish up"
    dfx         = CumTrap(d2fx,fpars$delta)
    c2          = min(c(0,dfx))
    fx          = CumTrap(dfx,fpars$delta)
    fx          = fpars$iflagpn*(fx - c2*(xgrid-xmin))
    if(fpars$iflagCenter==0){  # do not center f
      return(fx)
    }else{  # center f
      # Compute mean
      c         = trapint(fx,fpars$delta)
      fx        = fx - c/fpars$range
      return(fx)
    }
  }  # End GetSf0xgrid
  
  Getf0x = function(f0pars){
    fx     = 0
    # free function
    if(f0pars$iflagsc == 0)  return(f0pars$phix %*% f0pars$theta)
    # monotone function
    if(f0pars$iflagsc == 1){
      fx         = GetMonof0xgrid(f0pars)
      return(fx)
    } 
    # Increasing and convex or decreasing and concave
    if(f0pars$iflagsc == 2){ 
      fx         = GetMonoSamef0xgrid(f0pars)
      return(fx)
    } # End Increasing and convex or decreasing and concave
    # Increasing and concave or decreasing and convex
    if(f0pars$iflagsc == 3){ 
      fx         = GetMonoDiff0xgrid(f0pars)
      return(fx)
    } # End Increasing and concave or decreasing and convex
    # U-shaped function
    if(f0pars$iflagsc == 4){ 
      fx         = GetUf0xgrid(f0pars)
      return(fx)
    } # End u-shaped funtion
    # S-shaped (increasing,convex-to-concave) or (decreasing,concave-to-convex)
    if(f0pars$iflagsc == 5){
      fx  = GetSf0xgrid(f0pars,1)
      return(fx)
    }
    # S-shaped (increasing,concave-to-convex) or (decreasing,convex-to-concave)
    if(f0pars$iflagsc == 6){
      fx  = GetSf0xgrid(f0pars,-1)
      return(fx)
    }
  } #End Getf0x
  
  GetSCfxobs = function(xgrid,xinx,xover,fxgrid){
    xdelta = xgrid[2] - xgrid[1]
    n      = length(xinx)
    fxobs  = matrix(0,n,1)
    idn0   = which(xover>0)
    if(sum(idn0)<n){  # Got data on grid points
      fxobs[-idn0] = fxgrid[xinx[-idn0]]
    }else{  # overshoot
      y0     = fxgrid[xinx[idn0]]
      y1     = fxgrid[(xinx[idn0]+1)]
      fxobs[idn0] = y0 + (y1-y0)*xover[idn0]/xdelta  # Linearly interpolate
    }
    fxob   = iflagpn*fxobs
    return(fxobs)
  }
  
  ############### End of Local Functions ###############################################
  
  nparx       = 0      # E[Z(x)] = mu(x) = sum_{j=1}^nparx theta_{j} h(x)
  # If free, no constant, and nparx >= 0
  if(iflagsc>0 & nparx==0) nparx = 1  # make sure you get the constant
  
  
  # Generate error variance sigma: homogeneous or heterogeneous model
  if(iflagHBsigma == 0){
    # Homogenous error variance
    sigmat  = .5
  }else{
    # HB error variance
    # sigma_j^2 is IG(alpha/2,beta/2)
    sigma2_mu     = 1
    sigma2_v      = .5
    sigma2_v2     = sigma2_v^2
    sigma2_alphat = 2*(sigma2_mu^2/sigma2_v2 + 1)
    sigma2_betat  = 2*sigma2_mu*(sigma2_alphat/2 -1 )
    sigma2t       = sigma2_betat/(2*rgamma(ngroup,shape=sigma2_alphat/2,rate=1))
    sigmat        = sqrt(sigma2t)     # error std dev
  }
  
  
  
  # gamma_j ~ G(gamma_alpha, gamma_beta)
  gamma_mu     = .8   # upper-level model smoothing parameter
  gamma_sigma  = .3
  gamma_sigma2 = gamma_sigma^2
  gamma_alphat = gamma_mu/gamma_sigma2  
  gamma_betat  = gamma_alphat/gamma_mu
  # Smoothing parameters for lower-level spectral model
  gammallt    = rgamma(ngroup,shape=gamma_alphat,rate=gamma_betat)
  # Upper level smoother
  gamma0t    = gamma_mu
  
  
  
  # tau_k^2 ~ IG(tau2_r/2, tau2_s/2)
  tau2_mu    = 1     # Mean
  tau2_sd    = .1     # Standard deviation
  tau2_r     = 2*(2+(tau2_mu/tau2_sd)^2)
  tau2_s     = (tau2_r-4)*tau2_mu
  
  
  # Upper level: theta_0k ~ N(0,eta02*exp(-gamma0*log(1+k)))
  eta02t       = 1
  eta0t        = sqrt(eta02t)
  
  
  
  # Number of basis functions + parametric model
  ntheta = nbasis + nparx  # Total number of Z functions, including parametric model
  
  # Make groups for HB
  if(iflagRandn == 1){  # random group membership
    idHB1     = 1 + floor(ngroup*runif(ntot))   # idHB1 gives group membership
  }else{   # same number of observations in each group
    nobs    = floor(ntot/ngroup)
    ntot    = nobs*ngroup
    idHB1   = rep(1:ngroup,nobs)
  }
  
  uid       = sort(unique(idHB1))
  if(ngroup != length(uid)) print("Something wrong with HB groups")
  idout     = MakeHBPointer(idHB1)             # Pointers into groups
  nobs      = idout$nobs
  iptHB1    = idout$iptd
  
  iptHB2 = vector("list", length = length(iptHB1)) # index for c++
  for (i in 1:length(iptHB1)) {
    iptHB2[[i]] = iptHB1[[i]]-1  
  }
  
  # Got fixed effects
  if(nparv>0){
    # Generate v data
    vdata          = matrix(rnorm(ntot*nparv),ntot,nparv)
    vnames         = paste("V",1:nparv,sep="")
    alphat         = matrix(floor(.1*rnorm(nparv)*1000)/1000)
    valphat        = vdata %*% alphat
  }else{
    valphat        = 0
    vnames         = "Error!"
  }
  
  
  # Generate Z data
  if (nparz > 0 && nparw > 0) {
    zdata             = matrix(1,ngroup,nparz)
    if (nparz > 1) {
      zdata[,2:nparz] = rnorm(ngroup*(nparz-1))
      znames          = paste("Z",1:(nparz-1),sep="")
      znames          = c("CNST",znames)  
    } else { znames="CNST" }
    colnames(zdata) = znames
    # Generate Wdata 
    wdata    = matrix(1,ntot,nparw)
    if(nparw>1){
      for(j in 1:(nparw-1)){
        wdata[,(j+1)] = rnorm(ntot)
        wnames    = paste("W",1:(nparw-1),sep="")
        wnames    = c("CNST",wnames)
        colnames(wdata) = wnames
      }
    }else{
      wnames           = "CNST"
      colnames(wdata) = wnames
    } 
    
    # beta_j = phi'z_j + delta_j and delta_j ~ N(0,lambda)
    phit = rnorm(nparw*nparz)/3
    phit = floor(phit*1000)/1000
    phit = matrix(phit,nparz,nparw)
    colnames(phit)  = wnames
    row.names(phit) = znames
    # Generate error variance
    if(nparw==1){
      lambdat   = .1
      lambdat12 = sqrt(lambdat)
    }else{
      v        = (.01 + .3*runif(nparw))^2
      lambdat  = diag(v)
      lambdat12 = chol(lambdat)  # Upper diagonal: t(lambdat12)  %*% lambdat12 = lambdat
      colnames(lambdat) = wnames
    }
    
    
    # Generate beta
    berrors = matrix(rnorm(ngroup*nparw),ngroup,nparw) %*% lambdat12
    betat   = zdata%*%phit + berrors
    colnames(betat) = wnames
  } else {
    wdata    = matrix(0,ntot,1)
  }
  
  
  
  
  # Generate X data on [xmin, xmax]
  if(iflagSpanX > 0){   # X is uniform across the range of x
    xdata    = xmin + (xmax-xmin)*matrix(runif(ntot),ntot,1)
  }else{   # Data for different groups are concentrated at different areas
    xdata     = matrix(0,ntot,1)
    xrangeall = matrix(0,ngroup,1)
    xmid      = (xmin+xmax)/2   # midpoint of xgrid
    xq1       = (3*xmin + xmax)/4
    xq3       = (xmin + 3*xmax)/4
    p1        = .25  # p1 from xmin to xmid
    p2        = .25  # p2 from xmid to xmax
    p3        = .25  # p3 from xq1 to xq3
    # p4 from xmin to xq1 and xq3 to xmax
    for(j in 1:ngroup){
      u      = runif(1)
      if(u < p1){  # xmin to xmid
        xj           = xmin + (xmid-xmin)*matrix(runif(nobs[j]))
        xrangeall[j] = 1
        
      }else if(u < p1 + p2){ # xmid to xmax
        xj           = xmid + (xmax-xmid)*matrix(runif(nobs[j]))
        xrangeall[j] = 2
      }else if(u< p1 + p2 + p3){  # xq1 to xq3
        xj           = xq1 + (xq3-xq1)*matrix(runif(nobs[j]))
        xrangeall[j] = 3
      }else{  # xmin to xq1 or xq3 to xmax
        nj           = nobs[j]
        n1           = floor(nj/2)
        n2           = nj - n1
        xj1          = xmin + (xq1-xmin)*matrix(runif(n1))
        xj2          = xq3  + (xmax-xq3)*matrix(runif(n2))
        xj           = c(xj1,xj2)
        xrangeall[j] = 4
      }
      xdata[iptHB1[[j]]] = xj
    } # End loop over groups
  }
  
  #-----------------------------------------------------
  # Grid on [xmim,xmax] for plotting functions
  # Also, grid for numerical integration
  
  xdelta = (xmax-xmin)/nint
  xgrid  = seq(xmin,xmax,xdelta)
  xgrid  = matrix(xgrid,nint+1,1)
  xrange = xmax-xmin
  xmid   = (xmin+xmax)/2
  zgrid  = (xgrid-xmid)/xrange  # Standardized
  
  #-------------------------------------------------------
  # Use f on xgrid + interpolation to compute f at xdata
  # It is much faster than computing the intergram from xmin to x at each observation
  # Find location of xdata in xgrid
  # xgrid[xinx[i]] < xi <= xgrid[xinx[i]+1]
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
  
  # Get excess of xdata over boundry of xgrid based on xinx
  # Used in numerical integration
  xover = xdata - xgrid[xinx]
  idn0  = which(xover>0)  # which xdata is not at grid
  
  #------------------------------------------------------------
  # Vector used in computing variance of spectral coefficients
  # Var(theta_jk) = tau_k^2*exp(-k*gamma_j)
  # Population j and frequency k
  # Trick
  # Note that the if nparx > 0, then the first nparx parameters have variance = tau_k^2
  # So the first pad the first nparx elements of ktheta with zeros.
  kall     = matrix(1:nbasis)  # number of cosine functions, including intercept
  lkall    = log(1+kall)
  ktheta   = kall
  if(nparx > 0){
    ktheta   = matrix(c(rep(0,nparx),kall))
  }
  
  #---------------------------------------------------
  # Compute Cosine functions on xgrid and xdata 
  phi0xgrid = matrix(1/sqrt(xrange),nint+1,1)
  phixgrid  = sqrt(2/xrange)*cos(pi*xgrid%*%(t(kall)-xmin)/xrange)
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
    
  }else{
    # Free.  No constant, but may have parametric functions for mu(x)
    if(nparx > 0)  
      for(k in 1:nparx){
        phixgrid = cbind(zgrid^k,phixgrid)
      }
  }
  
  # If free function, compute phi at xdata
  # Compute Cosine functions on xgrid and xdata 
  if(iflagsc==0){
    zds     = (xdata-xmid)/xrange
    phixdata  = sqrt(2/xrange)*cos(pi*xdata%*%(t(ktheta)-xmin)/xrange)
    if(nparx > 0){
      for(k in 1:nparx){
        phixdata = cbind(zds^k,phixdata)
      }
    }
  }else{
    phixdata = 0
  }
  
  hfunall      = matrix(1,nint+1,ngroup)
  #Squish function's parameters for U shaped
  if(iflagsc>=4){
    # U-shaped f
    zeta0t     = log(.5/.5)
    zeta_s0     = .1
    zetallt    = zeta0t + zeta_s0*rnorm(ngroup)
    ez         = exp(zetallt)
    ez0        = exp(zeta0t)
    omegallt   = xmin + xrange*ez/(1+ez)
    omega0t    = xmin + xrange*ez0/(1+ez0)
    psiallt    = matrix(75,ngroup,1)
    if(iflagpsi>0){
      psit_log_m0 = log(75)
      psit_log_s0 = 1
      psiallt_log = psit_log_m0 + psit_log_s0*rnorm(ngroup)
      psiallt     = exp(psiallt_log)
    }
    # Compute squish functions
    hfunall      = matrix(0,nint+1,ngroup)
    for(j in 1:ngroup){
      hfunall[,j] = GetSquish(omegallt[j],psiallt[j],xgrid)
    }
  }
  
  
  # Get variances for HB theta
  # tau_k^2 ~ IG(r/2,s/2)
  tau2allt    = tau2_s/(2*rgamma(ntheta,shape=tau2_r/2,rate=1))
  
  # Need to make tau_0 smaller than other tau because it theta_0 is log normal
  if(iflagsc > 0) tau2allt[1] = .1   # variance for log(theta_j0) = w_j0
  tauallt     = sqrt(tau2allt)
  
  
  
  
  # Generate Fourier coefficients for upper level model
  vark0         = (1+kall/gamma_betat)^(-gamma_alphat)
  #vark0         = exp(-lkall*gamma0t)
  if(nparx>0) vark0 = c(rep(1,nparx),vark0)
  vark0         = eta02t*vark0
  theta0t       = sqrt(vark0)*rnorm(ntheta) 
  theta0t[1]    = abs(theta0t[1])
  xipar0        = log(theta0t[1])
  
  
  theta0t[1]      = 3
  theta0t[2]      = 1
  theta0t[3]      = -1
  theta0t[4]      = -.5
  theta0t[5]      = 0
  if(iflagsc>0){
    xipar0t       = 1
    theta0t[1]    = exp(xipar0t)
    xiparallt     = xipar0t + .1*rnorm(ngroup,1)
  }
  
  
  # Generate Fourier coefficients for lower level model
  a          = matrix(tau2allt,ntheta,ngroup)
  b          = exp(-(matrix(ktheta) %*% t(matrix(gammallt))))
  
  thv        = a * b    # Element by element multiplication
  thsd       = sqrt(thv)
  z          = matrix( rnorm(ngroup*ntheta),ntheta,ngroup)
  thetallt   = matrix(theta0t,ntheta,ngroup) + thsd*z
  if(iflagsc>0)  thetallt[1,] = matrix(exp(xiparallt),1,ngroup)
  
  
  # Generate functions and f at observations
  fxgridt = matrix(0,nint+1,ngroup)
  fxobst  = matrix(0,ntot,1)
  # Loop over groups
  for(j in 1:ngroup){
    # Get function for group j
    fpars  = list(iflagsc=iflagsc,iflagpn=iflagpn,iflagCenter=iflagCenter,
                  theta=thetallt[,j],phix=phixgrid,
                  delta=xdelta, range = xrange, hfun=hfunall[,j])
    fxj         = Getfx(fpars)
    fxgridt[,j] = fxj
    # Generate observations
    a   = matrix(xinx[iptHB1[[j]]])
    b   = matrix(xover[iptHB1[[j]]])
    c   = matrix(fxj)
    # Compute f at observations
    if(iflagsc>0){  # Shape constraint: use interpolation based of fxgrid
      fxobsj = GetSCfxobs(xgrid,a,b,c)
    }else{  # free model
      fxobsj = phixdata[iptHB1[[j]],] %*% matrix(thetallt[,j])
    }
    fxobst[iptHB1[[j]]] = matrix(fxobsj)
  }
  
  # Compute upper-level functions
  hfunm = 1
  # Compute E[hfun]
  if(iflagsc >= 4){
    hfunm = apply(hfunall,1,mean)
  }
  if(iflagsc>0){
    
    fpara            = theta0t    # Upper level spectral coefficients
    vpar             = (1/(1+kall/gamma_betat))^gamma_alphat
    if(nparx>0) vpar = c(rep(1,nparx),vpar)
    thvb             = tau2allt*vpar
    # If Z2, theta_0 = exp(xipar)  Need to adjust theta and Var(theta)
    if(iflagZ==0){
      fpara[1]         = fpara[1]*exp(tau2allt[1]/2)  # log normal theta_{0,j}
      etv0             = exp(tau2allt[1])
      thvb[1]          = exp(2*xipar0t  + tau2allt[1])*(etv0-1)  # lognormal variance V(theta_j0)
    }
    #thvb             = apply(thetallt,1,sd)^2
    thvb             = matrix(thvb,ntheta,1)
    
    f0pars = list(iflagsc = iflagsc, iflagpn = iflagpn, iflagCenter = iflagCenter,
                  theta = fpara,vpara = vpar, phix = phixgrid, phix2 = phi2xgrid, 
                  range = xrange, delta = xdelta, hfun = hfunm)
    f0xgridt = Getf0x(f0pars)
    
    
  }else{
    f0xgridt = phixgrid %*% theta0t
  } # End generate upper level f0 on xgrid
  #---------------------------------------------------------
  
  #----------------------------------------------
  # Compute Wi*beta_i
  
  wbetat = matrix(0,ntot,1)
  rerr   = matrix(0,ntot,1)
  for(j in 1:ngroup){ # Loop over populations
    bi   = matrix(betat[j,])
    wi   = wdata[iptHB1[[j]],]
    wibi = wi %*% bi
    wbetat[iptHB1[[j]],] = wibi
    # Homogeneous or HB error variance
    if(iflagHBsigma == 0){
      rerr[iptHB1[[j]]]  = sigmat*rnorm(nobs[j])
    }else{
      rerr[iptHB1[[j]]]  = sigmat[j]*rnorm(nobs[j])
    }
  }
  
  ydata  = valphat + wbetat + fxobst + rerr
  yresid = fxobst + rerr
  
  
  if(iflaglm == 0){
    # Normal model
    if(nparv>0){
      datall = data.frame(cbind(idHB1,ydata,xdata,vdata,wdata))
      dnames = c("Group","Y","X",vnames,wnames)
      colnames(datall) = dnames
    }else{
      datall = data.frame(cbind(idHB1,ydata,xdata,wdata))
      dnames = c("Group","Y","X",wnames)
      colnames(datall) = dnames
    }
    
  }
  
  
  if(iflaglm == 1){
    # Probit model
    # Convert latent ydata to 0/1 and save in ydata0
    ydata0 = ifelse(ydata >= 0,1,0)   # Make 0/1 data
    if(nparv>0){
      datall = data.frame(cbind(idHB1,ydata,xdata,vdata,wdata))
      dnames = c("Group","Y","X",vnames,wnames)
      colnames(datall) = dnames
    }else{
      datall = data.frame(cbind(idHB1,ydata,xdata,wdata))
      dnames = c("Group","Y","X",wnames)
      colnames(datall) = dnames
    }
  }
  
  if (disp==T) {
    
    fmin = apply(fxgridt,1,min)
    fmax = apply(fxgridt,1,max)
    favg = apply(fxgridt,1,mean)
    matplot(xgrid,cbind(f0xgridt,favg,fmin,fmax),type="l",
            main="f0, Mean f, Min f, Max f",
            xlab="X",ylab="f")
    
    matplot(xgrid,cbind(f0xgridt,fxgridt),type="l",
            main="f0 and fj for all j",
            xlab="X",ylab="f")
    
    
    ymin  = floor(min(c(f0xgridt,fxgridt,yresid)))
    ymax  = ceiling(max(c(f0xgridt,fxgridt,yresid)))
    matplot(xgrid,cbind(f0xgridt,fxgridt),type="l",ylim=c(ymin,ymax),
            xlab="X",ylab="Parametric Residuals")
    points(xdata,yresid)
    
    matplot(ktheta,cbind(theta0t,thetallt),type="l",xlab="Frequency",ylab="Theta")
    
    thm = matrix(apply(thetallt,1,mean))
    matplot(ktheta,cbind(theta0t,thm),type="l",xlab="Frequency",ylab="Theta",
            main="theta_0 & Mean theta")
    
    matplot(xgrid,cbind(f0xgridt,favg),type="l",
            main="f0 & Mean f",
            xlab="X",ylab="f")
  }
  
  #---------------------------------------------------------------
  out = list()
  
  out$datall = datall
  out$zdata  = zdata
  
  ftall           = data.frame(cbind(xgrid,f0xgridt,fxgridt))
  ftnames         = paste("ft",0:ngroup,sep="")
  colnames(ftall) = c("X",ftnames)
  out$TrueFun = ftall
  
  tallt                = data.frame(cbind(ktheta,theta0t,thetallt))
  tnames               = paste("Theta",1:ngroup,sep="_")
  colnames(tallt)      = c("K","Theta_Pop",tnames)
  out$TrueFourierCoeff = tallt
  
  out$phit=phit
  out$lambdat=lambdat
  out$sigmat=sigmat
  if(nparv>0){
    out$alphat=alphat
  }
  out$betat=betat
  out$thetallt=thetallt
  
  if(iflagsc>=4){
    out$omegallt=omegallt
    out$psiallt=psiallt
  }
  
  out$valphat = valphat
  out$wbetat  = wbetat
  out$fxobst  = fxobst
  out$rerr    = rerr
  
  
  parall = c(nparx,nparv,nparw,xmin,xmax,iflagsc,
             iflagpn,iflagpsi,iflaglm,iflagCenter,iflagZ,iflagHBsigma,nbasis)
  out$parall=parall
  
  return(out)
}