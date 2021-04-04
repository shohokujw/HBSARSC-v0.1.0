rm(list = ls())
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

library(devtools)
library(hbsamGP)
data("comp")

head(comp)

ydata    = comp$`Log(Sales/EMP) 2018`
xdata    = comp$`Log(XAD/EMP) 2017`
wdata    = comp$`Log EMP 2018`
group    = comp$Group
id_group = c(factor(group, labels = 1:length(unique(group))))
##################################################################################################
set.seed(1)
fout1 = hbsarv10(y=ydata, w=wdata, x=xdata, nbasis=20, id_group=id_group, nint=500, 
                mcmc = list(nblow=100000,smcmc=50000), prior = list(), shape = 'Free', iflagCenter=0)
fout1$mcmctime

##################################################################################################
set.seed(1)
fout2 = hbsarv10(y=ydata, w=wdata, x=xdata, nbasis=20, id_group=id_group, nint=500, 
                 mcmc = list(nblow=100000,smcmc=50000), prior = list(),
                 shape = 'Increasing', iflagCenter=1)
fout2$mcmctime
save.image("make_figure_9_10.rdata")

ugroups = unique(id_group)
ngroup  = length(ugroups)
gname   = unique(comp$`Group Name`)


# Fig9-a
pdf("fig9_a.pdf")
matplot(fout1$xgrid, fout1$post.est$fxgridm,type="l", ylim=c(min(fout1$post.est$fxgridm),3),
        xlab="Log(XAD/EMP) 2017", ylab="Adjusted Log(Sales/EMP) 2018", lwd=1)
dev.off()

# Fig9-b
pdf("fig9_b.pdf")
minfj = apply(fout1$post.est$fxgridm, 1, min)
maxfj = apply(fout1$post.est$fxgridm, 1, max)
plot(fout1$xgrid, fout1$post.est$f0xgridm, type="l", xlab="x", ylab="f", lwd=2,
     ylim = range(c(minfj, maxfj)))
lines(fout1$xgrid, apply(fout1$post.est$fxgridm, 1, mean), lwd=2, col="red", lty=2)
lines(fout1$xgrid, minfj, lty=3, lwd=2, col="green")
lines(fout1$xgrid, maxfj, lty=4, lwd=2, col="blue")
dev.off()


# Fig9-a
pdf("fig9_c.pdf")
matplot(fout2$xgrid, fout2$post.est$fxgridm,type="l", ylim=c(min(fout2$post.est$fxgridm),3),
        xlab="Log(XAD/EMP) 2017", ylab="Adjusted Log(Sales/EMP) 2018", lwd=1)
dev.off()

# Fig9-b
pdf("fig9_d.pdf")
minfj = apply(fout2$post.est$fxgridm, 1, min)
maxfj = apply(fout2$post.est$fxgridm, 1, max)
plot(fout2$xgrid, fout2$post.est$f0xgridm, type="l", xlab="x", ylab="f", lwd=2,
     ylim = range(c(minfj, maxfj)))
lines(fout2$xgrid, apply(fout2$post.est$fxgridm, 1, mean), lwd=2, col="red", lty=2)
lines(fout2$xgrid, minfj, lty=3, lwd=2, col="green")
lines(fout2$xgrid, maxfj, lty=4, lwd=2, col="blue")
dev.off()


# Fig10
# Fig10-a, b
mf = apply(fout1$post.est$fxgridm, 1, mean)
resid  = ydata - fout1$post.est$wbetam
fname = c("a", "b")
fig10_idx = c(3, 31) # Capital Goods, Pharmaceuticals
for(i in 1:length(fig10_idx)){
  j = fig10_idx[i]
  bool = (id_group == j)
  rj = resid[bool]
  xj = xdata[bool]
  
  fxgridj = fout1$post.est$fxgridm[,j]
  pname   = paste("f ",ugroups[j],sep="")
  ymax    = ceiling(max(c(rj,fxgridj,mf)))
  ymin    = floor(min(c(rj,fxgridj,mf)))
  bout    = cbind(fxgridj,mf)
  pdf(paste0("fig10_",fname[i],".pdf"))
  matplot(fout1$xgrid, bout,type="l",ylim=c(ymin,ymax), xlab="Age",ylab=gname[j], lwd=1,lty=c(1,2))
  points(xj,rj,pch=19)
  dev.off()
}

# Fig10-c, d
mf = apply(fout2$post.est$fxgridm, 1, mean)
resid  = ydata - fout2$post.est$wbetam
fname = c("c", "d")
fig10_idx = c(3, 31) # Capital Goods, Pharmaceuticals
for(i in 1:length(fig10_idx)){
  j = fig10_idx[i]
  bool = (id_group == j)
  rj = resid[bool]
  xj = xdata[bool]
  
  fxgridj = fout2$post.est$fxgridm[,j]
  pname   = paste("f ",ugroups[j],sep="")
  ymax    = ceiling(max(c(rj,fxgridj,mf)))
  ymin    = floor(min(c(rj,fxgridj,mf)))
  bout    = cbind(fxgridj,mf)
  pdf(paste0("fig10_",fname[i],".pdf"))
  matplot(fout2$xgrid, bout,type="l",ylim=c(ymin,ymax), xlab="Age",ylab=gname[j], lwd=1,lty=c(1,2))
  points(xj,rj,pch=19)
  dev.off()
}


