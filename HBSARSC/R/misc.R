"vech"=function (mat) 
{# Reshapes the lower triangular portion of a symmetric matrix
  # into a column vector, 'RSGHB' R package
  nr <- dim(mat)[1]
  nc <- dim(mat)[2]
  rV <- rep(0, nr * (nc + 1)/2)
  k <- 1
  for (i in 1:nr) {
    for (j in 1:i) {
      if (j <= i) {
        rV[k] = mat[i, j]
        k <- k + 1
      }
    }
  }
  return(as.matrix(rV))
}

"CosFun"=function(x,k,xmin,xrange) 
{
  tmp=sqrt(2/xrange)*cos(matrix(k,ncol=1)%*%matrix((pi*(x-xmin)/xrange),nrow=1))
  return(t(tmp))
}

#---------------------------------------------------
# MakeHBPointer
#    Make a pointer into the matrics for groups
#    Input
#        id   = identifies groups
#    Output
#        nobs = number of observations in each group
#        iptd = list of pointers into group
#    Usage
#        xj   = xdata[iptd[[j]],]  gives you the x data for group j
#    Call
#        outpt = MakeHBPointer(id)
#---------------------------------------------------
"MakeHBPointer" = function(id){
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


"GetMonof"<-function(theta,phixobs,quadfacts, fpm)
{
  return(fpm*QuadMult(theta, phixobs, quadfacts))
}

"QuadMult"<-function(x,qvech,quadfacts)
{
  colSums(quadfacts[,1]*x[quadfacts[,2]]*qvech*x[quadfacts[,3]]);
}


'function.shape'=function(shape=c('Free','Increasing','Decreasing','IncreasingConvex','DecreasingConcave',
                                  'IncreasingConcave','DecreasingConvex','IncreasingS','DecreasingS',
                                  'IncreasingRotatedS','DecreasingRotatedS','InvertedU','Ushape'))
{
  choices=c('Free','Increasing','Decreasing','IncreasingConvex','DecreasingConcave',
            'IncreasingConcave','DecreasingConvex','InvertedU','Ushape','IncreasingS','DecreasingS',
            'IncreasingRotatedS','DecreasingRotatedS')
  shape=match.arg(shape,choices,several.ok=TRUE)
  nfun=length(shape)
  fmodel=numeric(nfun)
  fpm=numeric(nfun)
  for(i in 1:nfun){
    switch(shape[i],
           Free={fmodel[i]=0;fpm[i]=1},
           Increasing={fmodel[i]=1;fpm[i]=1},
           Decreasing={fmodel[i]=1;fpm[i]=-1},
           IncreasingConvex={fmodel[i]=2;fpm[i]=1},
           DecreasingConcave={fmodel[i]=2;fpm[i]=-1},
           IncreasingConcave={fmodel[i]=3;fpm[i]=1},
           DecreasingConvex={fmodel[i]=3;fpm[i]=-1},
           InvertedU={fmodel[i]=4;fpm[i]=1},
           Ushape={fmodel[i]=4;fpm[i]=-1},
           IncreasingS={fmodel[i]=5;fpm[i]=1},
           DecreasingS={fmodel[i]=5;fpm[i]=-1},
           IncreasingRotatedS={fmodel[i]=6;fpm[i]=1},
           DecreasingRotatedS={fmodel[i]=6;fpm[i]=-1})
  }
  list(fmodel=fmodel,fpm=fpm,nfun=nfun)
}