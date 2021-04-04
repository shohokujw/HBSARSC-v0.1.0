rm(list = ls())
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

library(devtools)
library(hbsamGP)

# Social Medical Survey of Children Attending Child 
# Health Clinics (SMOCC (Herngreen et al., 1994))
library(brokenstick)
library(dplyr)

# ga : gestation time.
# bw : birth weight

data = brokenstick::smocc.hgtwgt
data <- data %>% select("id", "wgt", "age", "sex", "ga", "bw") %>% 
  na.omit() %>%
  mutate(sex = ifelse(sex=="male",1,-1)) %>%
  mutate(bw  = bw/1000)


ydata  = data$wgt
xdata  = data$age
wdata  = data[,c("sex","ga","bw")]
id_group = data[,"id"]

##################################################################################################
set.seed(1)
fout1 = hbsarv10(y=ydata, w=wdata, x=xdata, v=NULL, z=NULL,
                 nbasis=20, id_group=id_group, nint=500, mcmc = list(nblow=100000,smcmc=50000), prior = list(),
                 shape = 'Free', iflagCenter=0, iflagHBsigma = TRUE)
fout1$mcmctime


##################################################################################################
set.seed(1)
fout2 = hbsarv10(y=ydata, w=wdata, x=xdata, v=NULL, z=NULL,
                nbasis=20, id_group=id_group, nint=500, mcmc = list(nblow=100000,smcmc=50000), prior = list(),
                shape = 'Increasing', iflagCenter=1)
fout2$mcmctime

save.image("make_figure_7_8.rdata")

ugroups = unique(id_group)
ngroup  = length(ugroups)

# Fig5-a all
pdf(paste0("fig7_a.pdf"))
matplot(fout1$xgrid, fout1$post.est$fxgridm,type="l", xlab="Age", ylab="f(x)", lwd=1)
dev.off()

# Fig5-b
pdf(paste0("fig7_b.pdf"))
minfj = apply(fout1$post.est$fxgridm, 1, min)
maxfj = apply(fout1$post.est$fxgridm, 1, max)
plot(fout1$xgrid, fout1$post.est$f0xgridm, type="l", xlab="x", ylab="f", lwd=2,
     ylim = range(c(minfj, maxfj)))
lines(fout1$xgrid, apply(fout1$post.est$fxgridm, 1, mean), lwd=2, col="red", lty=2)
lines(fout1$xgrid, minfj, lty=3, lwd=2, col="green")
lines(fout1$xgrid, maxfj, lty=4, lwd=2, col="blue")
dev.off()



# Fig5-c
pdf(paste0("fig7_c.pdf"))
matplot(fout2$xgrid, fout2$post.est$fxgridm,type="l", xlab="Age", ylab="f(x)", lwd=1)
dev.off()

# Fig5-d 
pdf(paste0("fig7_d.pdf"))
minfj = apply(fout2$post.est$fxgridm, 1, min)
maxfj = apply(fout2$post.est$fxgridm, 1, max)
plot(fout2$xgrid, fout2$post.est$f0xgridm, type="l", xlab="x", ylab="f", lwd=2,
     ylim = range(c(minfj, maxfj)))
lines(fout2$xgrid, apply(fout2$post.est$fxgridm, 1, mean), lwd=2, col="red", lty=2)
lines(fout2$xgrid, minfj, lty=3, lwd=2, col="green")
lines(fout2$xgrid, maxfj, lty=4, lwd=2, col="blue")
dev.off()



# Fig6
# Fig6-a, b
mf = apply(fout1$post.est$fxgridm, 1, mean)
resid  = ydata - fout1$post.est$wbetam
fname = c("a", "b")
fig6_idx = c(11121, 11122)
for(i in 1:length(fig6_idx)){
  id = fig6_idx[i]
  j  = which(ugroups == id)
  bool = (id_group == id)
  rj = resid[bool]
  xj = xdata[bool]
  
  fxgridj = fout1$post.est$fxgridm[,j]
  pname   = paste("f ",ugroups[j],sep="")
  ymax    = ceiling(max(c(rj,fxgridj,mf)))
  ymin    = floor(min(c(rj,fxgridj,mf)))
  bout    = cbind(fxgridj,mf)
  pdf(paste0("fig8_",fname[i],".pdf"))
  matplot(fout1$xgrid, bout,type="l",ylim=c(ymin,ymax), xlab="Age",ylab=paste0("f",id), lwd=1,lty=c(1,2))
  points(xj,rj,pch=19)
  dev.off()
}

# Fig6-c, d
mf = apply(fout2$post.est$fxgridm, 1, mean)
resid  = ydata - fout2$post.est$wbetam
fname = c("c", "d")
fig6_idx = c(11121, 11122)
for(i in 1:length(fig6_idx)){
  id = fig6_idx[i]
  j  = which(ugroups == id)
  bool = (id_group == id)
  rj = resid[bool]
  xj = xdata[bool]
  
  fxgridj = fout2$post.est$fxgridm[,j]
  pname   = paste("f ",ugroups[j],sep="")
  ymax    = ceiling(max(c(rj,fxgridj,mf)))
  ymin    = floor(min(c(rj,fxgridj,mf)))
  bout    = cbind(fxgridj,mf)
  pdf(paste0("fig8_",fname[i],".pdf"))
  matplot(fout2$xgrid, bout,type="l",ylim=c(ymin,ymax), xlab="Age",ylab=paste0("f",id), lwd=1,lty=c(1,2))
  points(xj,rj,pch=19)
  dev.off()
}


