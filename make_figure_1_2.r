rm(list = ls())
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
library(devtools)
library(hbsamGP)

# 100 Observations per Group
set.seed(1)
sim = gHBSARSC(shape="Increasing", ntot = 10000, ngroup = 100, iflagSpanX = 1,
               nparv = 2, nparw = 3, nparz = 2, iflagCenter=1)

nparx        = sim$parall[1]              # Number of function in E[Z(x)], including intercept
nparv        = sim$parall[2]              # Number of fixed effects (no constant)
nparw        = sim$parall[3]              # Number of random effects (with constant)
xmin         = sim$parall[4]              # minimum range of x
xmax         = sim$parall[5]              # maximum range of x
iflagsc      = sim$parall[6]              # flag for shape constraint
iflagpn      = sim$parall[7]              # flag for increasing or decreasing
iflagpsi     = sim$parall[8]              # flag for fixed or hb model of psi
iflaglm      = sim$parall[9]              # Normal, Probit, ....
iflagCenter  = sim$parall[10]             # 1 mean int f = 0, 0 mean f(xmin) = 0
iflagZ       = sim$parall[11]             # 1 means df/dx = exp(Z) and 0 means Z^2
iflagHBsigma = sim$parall[12]
nbasis       = sim$parall[13]             # # of basis fuctions, exludes constant

datall   = sim$datall
zdata    = sim$zdata[,-1]                # drop intercept
ntot     = NROW(datall)
nc       = NCOL(datall)
id_group = as.matrix(datall[,1])
ydata    = as.matrix(datall[,2])
xdata    = as.matrix(datall[,3])
vdata    = as.matrix(datall[,4:(3+nparv)])
wdata    = as.matrix(datall[,(nc-nparw+2):nc])
set.seed(1)
##################################################################################################
fout = hbsarv10(y=ydata, w=wdata, x=xdata, z=zdata, nbasis=nbasis, id_group=id_group, nint=500, 
                mcmc = list(nblow=10000,smcmc=10000), prior = list(),
                shape = 'Increasing', iflagCenter=iflagCenter, iflagHBsigma=iflagHBsigma)
fout$mcmctime

'intsim' <- function(f,delta)
{# integration using Simpson's rule
  r=length(f)
  if(r==2*floor(r/2)){
    stop('ERROR: Even number of rows for simpson\'s integration')
  }else if(r==3){	# Simple Simpson's integration
    t=c(1,4,1)
    fint=sum(t*f)*delta/3
  }else{	# Composite Simpson's integration
    t=c(1,rep(c(4,2),times=(r-3)/2),4,1)
    fint=sum(t*f)*delta/3
  }
  fint
}

truef  = as.matrix(sim$TrueFun[,-c(1,2)])
estif  = fout$post.est$fxgridm
dx     = sim$TrueFun[2,1] - sim$TrueFun[1,1]
ngroup = length(unique(id_group))

RMISE = numeric(ngroup)
CORR  = numeric(ngroup)
for (j in 1:ngroup) {
  tmp = (truef[,j] - estif[,j])^2
  RMISE[j] = sqrt(intsim(tmp, dx))
  CORR[j]  = cor(truef[,j], estif[,j])
}

mean(RMISE)
sd(RMISE)

mean(CORR)
sd(CORR)


###############################################


ugroups = unique(id_group)
ngroup  = length(ugroups)

kall0  = 0:nbasis

## lower level spectral coefficient
pdf(paste0("fig1_a.pdf"))
matplot(kall0, fout$post.est$thetam, type="l",xlab="Frequency", ylab="theta", lwd=1)
dev.off()

## upper level spectral coefficient
pdf(paste0("fig1_b.pdf"))
plot(kall0, fout$post.est$theta0m, type="l", lwd=2, xlab="Frequency", ylab="theta")
lines(kall0, apply(fout$post.est$thetam,1,mean), lwd=2, col="red", lty=2)
dev.off()

## Lower-Level Functions
pdf(paste0("fig1_c.pdf"))
matplot(fout$xgrid, fout$post.est$fxgridm, type="l", xlab="xobs", ylab="f(x)", lwd=1)
dev.off()

## Upper-Level Functions
pdf(paste0("fig1_d.pdf"))
minfj = apply(fout$post.est$fxgridm, 1, min)
maxfj = apply(fout$post.est$fxgridm, 1, max)
plot(fout$xgrid, fout$post.est$f0xgridm, type="l", xlab="x", ylab="f", lwd=2,
     ylim = range(c(minfj, maxfj)))
lines(fout$xgrid, apply(fout$post.est$fxgridm, 1, mean), lwd=2, col="red", lty=2)
lines(fout$xgrid, minfj, lty=3, lwd=2, col="green")
lines(fout$xgrid, maxfj, lty=4, lwd=2, col="blue")
dev.off()

# id - Resudual
datatrue = sim$TrueFun
xgridt   = datatrue[,1]
ftrue    = datatrue[,-(1:2)]
mf = apply(fout$post.est$fxgridm, 1, mean)
xrange = range(xgridt)
resid  = ydata - fout$post.est$wbetam

# figue 2
fname = c("a", "b", "c", "d")
fig4_idx = c(1, 25, 50, 100)
for(i in 1:length(fig4_idx)){
  j = fig4_idx[i]
  bool   = (id_group == j)
  residj = resid[bool]
  xj     = xdata[bool]
  
  fj    = fout$post.est$fxgridm[,j]
  f0j   = fout$post.est$f0xgridm
  pname = paste("f ",ugroups[j],sep="")
  ymax  = ceiling(max(c(residj,fj,f0j)))
  ymin  = floor(min(c(residj,fj,f0j)))
  pdf(paste0("fig2_",fname[i],".pdf"))
  plot(fout$xgrid, fj,type="l",ylim=c(ymin,ymax), xlab="x",ylab=pname, lwd=2, col="black", lty=1)
  lines(fout$xgrid,f0j,col="green",lwd=2,lty=3)
  lines(xgridt,ftrue[,j],col="red",lwd=2,lty=2)
  points(xj,residj,pch=19)
  dev.off()
}


save.image("make_figure_1_2.rdata")
