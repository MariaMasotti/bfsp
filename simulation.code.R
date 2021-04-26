
library(Rfast);library(fields);library(FastGP);library(sparseMVN);library(extraDistr);library(Matrix)

job = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

cat(paste0("Processing Job #", job, "\n"))

################### functions #####################

a2c<-function(x,y,center){
  shift.x<-x-center[1]
  shift.y<-y-center[2]
  a2c<-atan2(shift.x,shift.y)
  a2c<-ifelse(a2c<0,a2c+2*pi,a2c)
  return(a2c)
}


d2c<-function(x,y,center){
  d2c<-sqrt((x-center[1])^2+(y-center[2])^2)
  return(d2c)
}

likelihood=function(param,param0,value,x,y,Sigma,X){
  XB<-X%*%param[5:15]
  cluster<-+(d2c(x,y,param[1:2])<=XB)
  log.lik<-sapply(0:1,calc.lik,value=value,Sigma,param0=param0,param=param,cluster=cluster)
  return(log.lik)
}

likelihood0=function(param,param0,value,x,y,Sigma,X){
  XB<-X%*%param[5:15]
  cluster<-+(d2c(x,y,param[1:2])<=XB)
  log.lik<-calc.lik(0,value,Sigma,param0,param,cluster)
  return(log.lik)
}

likelihood1=function(param,param0,value,x,y,Sigma,X){
  cluster<-+(d2c(x,y,param[1:2])<=X%*%param[5:15])
  log.lik<-calc.lik(1,value,Sigma,param0,param,cluster)
  return(log.lik)
}




calc.lik<-function(n,value,Sigma,param0,param,cluster){
  if(n==0){
    ind0<-which(cluster==0)
    Y0<-value[ind0]
    n0<-length(Y0)
    cov<-param0[3]*Sigma[ind0,ind0]+param0[2]*diag(n0)
    log.lik<- dmvn.sparse(t(Y0), rep(param0[1],length(Y0)), Cholesky(as(cov,"dsCMatrix")), prec = F, log = TRUE)
  }else{
    ind1<-which(cluster==1)
    Y1<-value[ind1]
    n1<-length(Y1)
    if(length(Y1)<2) return(-Inf)
    cov<-param[4]*Sigma[ind1,ind1]+param0[2]*diag(n1)
    log.lik<- dmvn.sparse(t(Y1), rep(param[3],length(Y1)), Cholesky(as(cov,"dsCMatrix")), prec = F, log = TRUE)
  }
  return(as.numeric(log.lik))
}


prior_param=function(param,param0,x,y,X,xlim,ylim){
  XB<-X%*%param[5:15]
  bp<-cbind(x=XB*sin(a2c(x,y,param[1:2]))+param[1],y=XB*cos(a2c(x,y,param[1:2]))+param[2])
  if(any(bp[,1]<xlim[1])|any(bp[,1]>xlim[2])|any(bp[,2]<ylim[1])|any(bp[,2]>ylim[2])){return(-Inf)}
  prior=sum(dgamma(1/param[4],.1,.1,log=T),na.rm=T)+sum(dhnorm(param[3]-param0[1], sigma = 10, log = TRUE),na.rm=T)+sum(dunif(param[5],min=0,max=.5,log=T),na.rm=T)+sum(dnorm(param[6:15],mean=0,sd=1,log=T),na.rm=T)
  return(prior)
}

prior_param0=function(param0){
  prior=(dnorm(param0[1],0,100,log=T)+sum(dgamma(1/param0[2:3],.1,.1,log=T)))
  return(prior)
}



param_update=function(param,param0,la,la1,la2,gam,gam1,gam2,i,value,x,y,ll,Sigma,X,xlim,ylim){
  
  if(i<=50) {
    param_star<-propose_new(param[3:4],(.01/2)*diag(2))
  }else{
    param_star<-.99*propose_new(param[3:4],exp(la1)*gam1)+.01*propose_new(param[3:4],(.01/2)*diag(2))}
  parameters_star<-param
  parameters_star[3:4]<-param_star
  prior_ratio=prior_param(parameters_star,param0,x,y,X,xlim,ylim)-prior_param(param,param0,x,y,X,xlim,ylim)
  if(prior_ratio>-Inf){
    ll_star<-likelihood1(parameters_star,param0,value,x,y,Sigma,X)
    like_ratio=(ll_star)-ll[2]
    alpha=exp(prior_ratio+like_ratio)
    if(runif(1)<alpha){
      param[3:4]<-param_star
      ll[2]<-ll_star
    } 
  }
  if(i<=50) {
    param_star<-propose_new(param[5:10],(.01/6)*diag(6))
  }else{
    param_star<-.99*propose_new(param[5:10],exp(la)*gam)+.01*propose_new(param[5:10],(.01/6)*diag(6))}
  parameters_star<-param
  parameters_star[5:10]<-param_star
  prior_ratio=prior_param(parameters_star,param0,x,y,X,xlim,ylim)-prior_param(param,param0,x,y,X,xlim,ylim)
  if(prior_ratio>-Inf){
    ll_star<-likelihood(parameters_star,param0,value,x,y,Sigma,X)
    like_ratio=sum(ll_star)-sum(ll)
    alpha=exp(prior_ratio+like_ratio)
    if(runif(1)<alpha){
      param[5:10]<-param_star
      ll<-ll_star
    }
  }
  if(i<=50) {
    param_star<-propose_new(param[11:15],(.01/5)*diag(5))
  }else{
    param_star<-.99*propose_new(param[11:15],exp(la2)*gam2)+.01*propose_new(param[11:15],(.01/5)*diag(5))}
  parameters_star<-param
  parameters_star[11:15]<-param_star
  prior_ratio=prior_param(parameters_star,param0,x,y,X,xlim,ylim)-prior_param(param,param0,x,y,X,xlim,ylim)
  if(prior_ratio>-Inf){
    ll_star<-likelihood(parameters_star,param0,value,x,y,Sigma,X)
    like_ratio=sum(ll_star)-sum(ll)
    alpha=exp(prior_ratio+like_ratio)
    if(runif(1)<alpha){
      param[11:15]<-param_star
      ll<-ll_star
    }
  }
  #}
  
  return(list(param,ll))
}


propose_new<-function(param,cov){
  prop<-rcpp_rmvnorm(1,cov,param)
  return(prop)
}

param0_update=function(param0,la0,gam0,i,param,value,x,y,ll,Sigma,X,xlim,ylim){
  if(i<=50) {
    param0_star<-c(propose_new(param0,(.01/3)*diag(3)))
  }else{
    param0_star<-c(.99*propose_new(param0,exp(la0)*gam0)+.01*propose_new(param0,(.01/3)*diag(3)))}
  
  prior_ratio=prior_param0(param0_star)+prior_param(param,param0_star,x,y,X,xlim,ylim)-prior_param0(param0)-prior_param(param,param0,x,y,X,xlim,ylim)
  if(prior_ratio>-Inf){ 
    ll_star<-likelihood0(param,param0_star,value,x,y,Sigma,X)
    like_ratio=ll_star-ll[1]
    alpha=exp(prior_ratio+like_ratio)
    if(runif(1)<alpha){
      param0<-param0_star
      ll[1]<-ll_star
    } 
  }
  return(list(param0,ll))
}

bfsp<-function(loc,value,iterations,N,inits,xlim,ylim){
  
  #calc Sigma 
  dist<-dist(loc, diag=T, upper=T)
  wend<-Wendland(dist, theta = .1, dimension=2, k=2)
  dist_mat<-as.matrix(dist)
  wend_mat<-as.matrix(wend)
  diag(wend_mat)<-1
  Dinv <- 1/dist_mat
  diag(Dinv)<-0
  Dinv <- Dinv/rowSums(Dinv)
  prec<-diag(dim(loc)[1])-.99*Dinv
  sig<-spdinv(prec)
  Sigma<-sig*wend_mat
  
  
  #set values 
  x<-loc[,1]
  y<-loc[,2]
  param0<-inits[[1]]
  param<-inits[[2]]
  X<-cbind2(cbind2(1,sin(outer(a2c(x,y,param[1:2]),(1:N)))),cos(outer(a2c(x,y,param[1:2]),(1:N))))
  ll<-likelihood(param,param0,value,x,y,Sigma,X)
  mu0<-param0
  mu<-param[5:10]
  mu1<-param[3:4]
  mu2<-param[11:15]
  gam0<-diag(3)
  gam=diag(6)
  gam1=diag(2)
  gam2=diag(5)
  param0_keeps<-array(dim=c(iterations,3))
  param_keeps<-array(dim=c(iterations,15))
  la0<--1
  la<--1
  la1<--1
  la2<--1
  
  
  
  
  start.time <- Sys.time()
  for(i in 1:iterations){
    
    
    #update centroid 
    XB<-X%*%param[5:15]
    clust<-+(d2c(x,y,param[1:2])<=XB)
    center_star=c(mean(x[clust==1]),mean(y[clust==1]))
    if((abs(center_star[1]-param[1])>.01|abs(center_star[2]-param[2])>.01)){
      bp<-cbind(x=XB*sin(a2c(x,y,param[1:2]))+param[1],y=XB*cos(a2c(x,y,param[1:2]))+param[2])
      theta_star<-a2c(bp[,1],bp[,2],center_star)
      f_star<-d2c(bp[,1],bp[,2],center_star)
      X_star<-cbind2(cbind2(1,sin(outer(theta_star,(1:N)))),cos(outer(theta_star,(1:N))))
      #if(!has_error(solve(t(X_star)%*%X_star))){
      beta_star<-solve(t(X_star)%*%X_star)%*%t(X_star)%*%f_star
      if(beta_star[1]<0){
        beta_star[1]<-0
      }
      param_star<-c(center_star,param[3:4],beta_star)
      X_star<-cbind2(cbind2(1,sin(outer(a2c(x,y,param_star[1:2]),(1:N)))),cos(outer(a2c(x,y,param_star[1:2]),(1:N))))
      prior_star<-prior_param(param_star,param0,x,y,X_star,xlim,ylim)
      ll_star<-likelihood( param_star,param0,value,x,y,Sigma,X_star)
      if(sum(ll_star)+sum(prior_star)>-Inf){
        param[5:15]<-beta_star
        param[1:2]<-center_star
        ll<-ll_star
        X<-X_star
      }
      #}
    }
    
    
    param0_move<-param0_update(param0,la0,gam0,i,param,value,x,y,ll,Sigma,X,xlim,ylim)
    param0<-param0_move[[1]]
    ll<-param0_move[[2]]
    param_move<-param_update(param,param0,la,la1,la2,gam,gam1,gam2,i,value,x,y,ll,Sigma,X,xlim,ylim)
    param<-param_move[[1]]
    ll<-param_move[[2]]
    
    
    
    # #plot to check progress 
    
    if(i %% 100== 0) cat(paste(100*i/iterations,"% complete","\n"))
    
    
    # #update covariance of parameters
    if(i>=50){
      mu0<-mu0+(1/(i+1))*(param0-mu0)
      gam0<-gam0+(1/(i+1))*((param0-mu0)%*%t(param0-mu0)-gam0)
      
      mu<-mu+(1/(i+1))*(param[5:10]-mu)
      mu1<-mu1+(1/(i+1))*(param[3:4]-mu1)
      mu2<-mu2+(1/(i+1))*(param[11:15]-mu2)
      gam<-gam+(1/(i+1))*((param[5:10]-mu)%*%t(param[5:10]-mu)-gam)
      gam1<-gam1+(1/(i+1))*((param[3:4]-mu1)%*%t(param[3:4]-mu1)-gam1)
      gam2<-gam2+(1/(i+1))*((param[11:15]-mu2)%*%t(param[11:15]-mu2)-gam2)
      
    }
    
    #keep things and calc average boundaries and centers 
    param0_keeps[i,]<-param0
    param_keeps[i,]<-param
    
    
    
    #update step  size 
    if(i%%20==0){
      last20<-param0_keeps[(i-20):i,]
      accept<-1-mean(duplicated(round(last20,8)))
      if(accept<.4) la0<-la0-.1
      if(accept>.4) la0<-la0+.1
      
      last20<-param_keeps[(i-20):i,5:10]
      accept<-1-mean(duplicated(round(last20,8)))
      if(accept<.2) la<-la-.1
      if(accept>.2) la<-la+.1
      
      last20<-param_keeps[(i-20):i,3:4]
      accept<-1-mean(duplicated(round(last20,8)))
      if(accept<.4) la1<-la1-.1
      if(accept>.4) la1<-la1+.1
      
      last20<-param_keeps[(i-20):i,11:15]
      accept<-1-mean(duplicated(round(last20,8)))
      if(accept<.2) la2<-la2-.1
      if(accept>.2) la2<-la2+.1
      
    }
    
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  list1<-list(param_keeps,param0_keeps,time.taken)
  
  list1
  
}


myfunc<-function(m){
  
  ##########generate simulated data############
  loc<-as.matrix(expand.grid(seq(0,1,length=50),seq(0,1,length=50)))
  xx<-loc[,1]
  yy<-loc[,2]
  
  
  theta1<-.05
  theta0<-.05
  smooths<-c(.5,1.5)
  nugs<-c(.1,.5)
  shape<-c("square","heart","triangle")
  grid<-expand.grid(1:50,smooths,nugs,shape)
  smooth1<-grid[m,2]
  smooth0<-1
  nug<-grid[m,3]
  shape<-grid[m,4]
  center<-c(.6,.6)
  
  if(shape=="triangle"){
    tru_cluster<-ifelse(yy<xx&yy<(1-xx)&yy>=.2,1,0)
    center=c(mean(xx[tru_cluster==1]),mean(yy[tru_cluster==1]))
  }
  
  if(shape=="square"){
    tru_cluster<-ifelse((xx<=(center[1]+.2)&xx>=(center[1]-.2)&yy>=(center[2]-.2)&yy<=(center[2]+.2)),1,0)
  }
  
  
  if(shape=="heart"){
    loc<-as.data.frame(expand.grid(seq(-center[1],1-center[1],length=50),seq(-center[2],1-center[2],length=50)))
    names(loc)<-c("x","y")
    xhrt <- function(t) .24*sin(t)^3
    yhrt <- function(t) .195*cos(t)-.075*cos(2*t)-.03*cos(3*t)-.015*cos(4*t)
    dat<- data.frame(t=seq(0, 2*pi, by=0.01) )
    dat$y=yhrt(dat$t)
    dat$x=xhrt(dat$t)
    cen<-c(0,0)
    thet<-a2c(dat$x,dat$y,cen)
    f<-d2c(dat$x,dat$y,cen)
    X<-cbind2(cbind2(1,sin(outer(thet,(1:100)))),cos(outer(thet,(1:100))))
    beta<-solve(t(X)%*%X)%*%t(X)%*%f
    theta<-a2c(loc[,1],loc[,2],cen)
    X<-cbind2(cbind2(1,sin(outer(theta,(1:100)))),cos(outer(theta,(1:100))))
    f<-X%*%beta
    bp.shape<-cbind.data.frame(x=f*sin(a2c(loc[,1],loc[,2],cen))+cen[1],y=f*cos(a2c(loc[,1],loc[,2],cen))+cen[2])
    dist.loc<-d2c(loc[,1],loc[,2],c(0,0))
    dist.h<-d2c(bp.shape[,1],bp.shape[,2],c(0,0))
    tru_cluster<-ifelse(dist.loc<=dist.h,1,0)
    range01 <- function(x){(x-min(x))/(max(x)-min(x))}
    xx<-range01(loc$x)
    yy<-range01(loc$y)
  }
  
  
  
  x<-xx
  y<-yy
  loc<-cbind.data.frame(x,y)
  loc0<-loc[tru_cluster==0,]
  loc1<-loc[tru_cluster==1,]
  n0=dim(loc0)[1]
  n1=dim(loc1)[1]
  
  Cov1<-stationary.cov(loc1, x2=NULL, Covariance = "Matern", theta=theta1,smoothness=smooth1)
  Cov0<-stationary.cov(loc0, x2=NULL, Covariance = "Matern", theta=theta0,smoothness=smooth0)
  
  Y0     <- rcpp_rmvnorm(1,Cov0+nug*diag(n0),rep(-1,n0))
  Y1     <- rcpp_rmvnorm(1,Cov1+nug*diag(n1),rep(2,n1))
  
  value<-array(dim=c(length(x),1))
  value[tru_cluster==0]<-Y0
  value[tru_cluster==1]<-Y1
  
  N<-5
  iterations<-50000
  inits<-list(param0=c(-1,1,1),param=c(.5,.5,1,1,.4,rep(0,2*N)))
  xlim<-c(0,1)
  ylim<-c(0,1)
  
  #run bfsp 
  bfsp_out<-bfsp(loc,value,iterations,N,inits,xlim,ylim)
  return(list(bfsp_out,x,y,value,tru_cluster))
}


set.seed(job)
result<- myfunc(job)
saveRDS(result , paste0("paper1.simulation.CAR.", job, ".rds"))



