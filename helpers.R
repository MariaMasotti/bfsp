

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

# Generate simiulated data 
simulated.matern.data<-function(smooth,nugget,shape){
  loc<-as.matrix(expand.grid(seq(0,1,length=50),seq(0,1,length=50)))
  xx<-loc[,1]
  yy<-loc[,2]
  smooth1<-smooth
  smooth0<-1
  
  
  if(shape=="triangle"){
    tru_cluster<-ifelse(yy<xx&yy<(1-xx)&yy>=.2,1,0)
  }
  
  center<-c(.6,.6)
  if(shape=="square"){
    tru_cluster<-ifelse((xx<=(center[1]+.2)&xx>=(center[1]-.2)&yy>=(center[2]-.2)&yy<=(center[2]+.2)),1,0)
  }
  
  if(shape=="circle"|shape=="Circle"){
    tru_cluster<-ifelse((xx-center[1])^2+(yy-center[2])^2<.25^2,1,0)
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
  loc0<-loc[tru_cluster==0,]
  loc1<-loc[tru_cluster==1,]
  n0=dim(loc0)[1]
  n1=dim(loc1)[1]
  
  Cov1<-stationary.cov(loc1, x2=NULL, Covariance = "Matern", theta=theta1,smoothness=smooth1)
  Cov0<-stationary.cov(loc0, x2=NULL, Covariance = "Matern", theta=theta0,smoothness=smooth0)
  
  Y0     <- rcpp_rmvnorm(1,Cov0+nugget*diag(n0),rep(-1,n0))
  Y1     <- rcpp_rmvnorm(1,Cov1+nugget*diag(n1),rep(2,n1))
  
  value<-array(dim=c(length(x),1))
  value[tru_cluster==0]<-Y0
  value[tru_cluster==1]<-Y1
  
  return(list(loc,value))
  
}




likelihood=function(param,param0,value,x,y,Sigma,X,N){
  XB<-X%*%param[5:(5+2*N)]
  cluster<-+(d2c(x,y,param[1:2])<=XB)
  log.lik<-sapply(0:1,calc.lik,value=value,Sigma,param0=param0,param=param,cluster=cluster)
  return(log.lik)
}

likelihood0=function(param,param0,value,x,y,Sigma,X,N){
  XB<-X%*%param[5:(5+2*N)]
  cluster<-+(d2c(x,y,param[1:2])<=XB)
  log.lik<-calc.lik(0,value,Sigma,param0,param,cluster)
  return(log.lik)
}

likelihood1=function(param,param0,value,x,y,Sigma,X,N){
  cluster<-+(d2c(x,y,param[1:2])<=X%*%param[5:(5+2*N)])
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


prior_param=function(param,param0,x,y,X,xlim,ylim,N){
  XB<-X%*%param[5:(5+2*N)]
  bp<-cbind(x=XB*sin(a2c(x,y,param[1:2]))+param[1],y=XB*cos(a2c(x,y,param[1:2]))+param[2])
  if(any(bp[,1]<xlim[1])|any(bp[,1]>xlim[2])|any(bp[,2]<ylim[1])|any(bp[,2]>ylim[2])){return(-Inf)}
  prior=sum(dgamma(1/param[4],.1,.1,log=T),na.rm=T)+sum(dhnorm(param[3]-param0[1], sigma = 10, log = TRUE),na.rm=T)+sum(dunif(param[5],min=0,max=.5,log=T),na.rm=T)+sum(dnorm(param[6:15],mean=0,sd=1,log=T),na.rm=T)
  return(prior)
}

prior_param0=function(param0){
  prior=(dnorm(param0[1],0,100,log=T)+sum(dgamma(1/param0[2:3],.1,.1,log=T)))
  return(prior)
}



param_update=function(param,param0,la,la1,la2,gam,gam1,gam2,i,value,x,y,ll,Sigma,X,xlim,ylim,N){
  
  if(i<=50) {
    param_star<-propose_new(param[3:4],(.01/2)*diag(2))
  }else{
    param_star<-.99*propose_new(param[3:4],exp(la1)*gam1)+.01*propose_new(param[3:4],(.01/2)*diag(2))}
  parameters_star<-param
  parameters_star[3:4]<-param_star
  prior_ratio=prior_param(parameters_star,param0,x,y,X,xlim,ylim,N)-prior_param(param,param0,x,y,X,xlim,ylim,N)
  if(prior_ratio>-Inf){
    ll_star<-likelihood1(parameters_star,param0,value,x,y,Sigma,X,N)
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
  prior_ratio=prior_param(parameters_star,param0,x,y,X,xlim,ylim,N)-prior_param(param,param0,x,y,X,xlim,ylim,N)
  if(prior_ratio>-Inf){
    ll_star<-likelihood(parameters_star,param0,value,x,y,Sigma,X,N)
    like_ratio=sum(ll_star)-sum(ll)
    alpha=exp(prior_ratio+like_ratio)
    if(runif(1)<alpha){
      param[5:10]<-param_star
      ll<-ll_star
    }
  }
  
  for(p in 2:(ceiling(2*N/5))){
    index<-(5*p+1):min((5*p+5),length(param))
  if(i<=50) {
    param_star<-propose_new(param[index],(.01/length(index))*diag(length(index)))
  }else{
    param_star<-.99*propose_new(param[index],exp(la2[[p-1]])*gam2[[p-1]])+.01*propose_new(param[index],(.01/length(index))*diag(length(index)))}
  parameters_star<-param
  parameters_star[index]<-param_star
  prior_ratio=prior_param(parameters_star,param0,x,y,X,xlim,ylim,N)-prior_param(param,param0,x,y,X,xlim,ylim,N)
  if(prior_ratio>-Inf){
    ll_star<-likelihood(parameters_star,param0,value,x,y,Sigma,X,N)
    like_ratio=sum(ll_star)-sum(ll)
    alpha=exp(prior_ratio+like_ratio)
    if(runif(1)<alpha){
      param[index]<-param_star
      ll<-ll_star
    }
  }
  }

  
  return(list(param,ll))
}


propose_new<-function(param,cov){
  prop<-rcpp_rmvnorm(1,cov,param)
  return(prop)
}

param0_update=function(param0,la0,gam0,i,param,value,x,y,ll,Sigma,X,xlim,ylim,N){
  if(i<=50) {
    param0_star<-c(propose_new(param0,(.01/3)*diag(3)))
  }else{
    param0_star<-c(.99*propose_new(param0,exp(la0)*gam0)+.01*propose_new(param0,(.01/3)*diag(3)))}
  
  prior_ratio=prior_param0(param0_star)+prior_param(param,param0_star,x,y,X,xlim,ylim,N)-prior_param0(param0)-prior_param(param,param0,x,y,X,xlim,ylim,N)
  if(prior_ratio>-Inf){ 
    ll_star<-likelihood0(param,param0_star,value,x,y,Sigma,X,N)
    like_ratio=ll_star-ll[1]
    alpha=exp(prior_ratio+like_ratio)
    if(runif(1)<alpha){
      param0<-param0_star
      ll[1]<-ll_star
    } 
  }
  return(list(param0,ll))
}

bfsp<-function(loc,value,iterations,burn,N,inits,xlim,ylim){
  
  #calc tapered car covariance matrix 
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
  
  
  #set up inits 
  x<-loc[,1]
  y<-loc[,2]
  param0<-inits[[1]]
  param<-inits[[2]]
  X<-cbind2(cbind2(1,sin(outer(a2c(x,y,param[1:2]),(1:N)))),cos(outer(a2c(x,y,param[1:2]),(1:N))))
  ll<-likelihood(param,param0,value,x,y,Sigma,X,N)
  mu0<-param0
  mu<-param[5:10]
  mu1<-param[3:4]
  mu2<-list()
  for(p in 2:(ceiling(2*N/5))){
    mu2[[p-1]]<-param[(5*p+1):min((5*p+5),length(param))]
  }
  gam0<-diag(3)
  gam=diag(6)
  gam1=diag(2)
  gam2=list()
  for(p in 2:(ceiling(2*N/5))){
    gam2[[p-1]]<-diag(length(mu2[[p-1]]))
  }
  param0_keeps<-array(dim=c(iterations,3))
  param_keeps<-array(dim=c(iterations,(5+2*N)))
  la0<--1
  
  la<--1
  la1<--1
  la2<-rep(list(-1),length(mu2))
  
  
  
  
  start.time <- Sys.time()
  #run mcmc 
  for(i in 1:iterations){
    
    
    #update centroid 
    XB<-X%*%param[5:(5+2*N)]
    clust<-+(d2c(x,y,param[1:2])<=XB)
    center_star=c(mean(x[clust==1]),mean(y[clust==1]))
    if((abs(center_star[1]-param[1])>.01|abs(center_star[2]-param[2])>.01)){
      bp<-cbind(x=XB*sin(a2c(x,y,param[1:2]))+param[1],y=XB*cos(a2c(x,y,param[1:2]))+param[2])
      theta_star<-a2c(bp[,1],bp[,2],center_star)
      f_star<-d2c(bp[,1],bp[,2],center_star)
      X_star<-cbind2(cbind2(1,sin(outer(theta_star,(1:N)))),cos(outer(theta_star,(1:N))))
      
      beta_star<-solve(t(X_star)%*%X_star)%*%t(X_star)%*%f_star
      if(beta_star[1]<0){
        beta_star[1]<-0
      }
      param_star<-c(center_star,param[3:4],beta_star)
      X_star<-cbind2(cbind2(1,sin(outer(a2c(x,y,param_star[1:2]),(1:N)))),cos(outer(a2c(x,y,param_star[1:2]),(1:N))))
      prior_star<-prior_param(param_star,param0,x,y,X_star,xlim,ylim,N)
      ll_star<-likelihood( param_star,param0,value,x,y,Sigma,X_star,N)
      if(sum(ll_star)+sum(prior_star)>-Inf){
        param[5:(5+2*N)]<-beta_star
        param[1:2]<-center_star
        ll<-ll_star
        X<-X_star
      }
      
    }
    
    
    param0_move<-param0_update(param0,la0,gam0,i,param,value,x,y,ll,Sigma,X,xlim,ylim,N)
    param0<-param0_move[[1]]
    ll<-param0_move[[2]]
    param_move<-param_update(param,param0,la,la1,la2,gam,gam1,gam2,i,value,x,y,ll,Sigma,X,xlim,ylim,N)
    param<-param_move[[1]]
    ll<-param_move[[2]]
    
    
    
    # #update covariance of parameters
    if(i>=50){
      mu0<-mu0+(1/(i+1))*(param0-mu0)
      gam0<-gam0+(1/(i+1))*((param0-mu0)%*%t(param0-mu0)-gam0)
      
      mu<-mu+(1/(i+1))*(param[5:10]-mu)
      mu1<-mu1+(1/(i+1))*(param[3:4]-mu1)
      for(p in 2:(ceiling(2*N/5))){
        mu2[[p-1]]<-mu2[[p-1]]+(1/(i+1))*(param[(5*p+1):min((5*p+5),length(param))]-mu2[[p-1]])
      }
      
      gam<-gam+(1/(i+1))*((param[5:10]-mu)%*%t(param[5:10]-mu)-gam)
      gam1<-gam1+(1/(i+1))*((param[3:4]-mu1)%*%t(param[3:4]-mu1)-gam1)
      for(p in 2:(ceiling(2*N/5))){
      gam2[[p-1]]<-gam2[[p-1]]+(1/(i+1))*((param[(5*p+1):min((5*p+5),length(param))]-mu2[[p-1]])%*%t(param[(5*p+1):min((5*p+5),length(param))]-mu2[[p-1]])-gam2[[p-1]])
      }
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
      
      for(p in 2:(ceiling(2*N/5))){
      last20<-param_keeps[(i-20):i,(5*p+1):min((5*p+5),length(param))]
      accept<-1-mean(duplicated(round(last20,8)))
      if(accept<.2) la2[[p-1]]<-la2[[p-1]]-.1
      if(accept>.2) la2[[p-1]]<-la2[[p-1]]+.1
      }
    }
    
    if(i %% 100== 0) cat(paste(100*i/iterations,"% complete","\n"))
    
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  list1<-list(param_keeps[-c(1:burn),],param0_keeps[-c(1:burn),],time.taken)
  
  list1
  
}

plot_bfsp<-function(result,loc,value,credible.band){
  param<-colMeans(result[[1]])
  N<-(length(param)-5)/2
  x<-loc[,1]
  y<-loc[,2]
  X<-cbind2(cbind2(1,sin(outer(a2c(x,y,param[1:2]),(1:N)))),cos(outer(a2c(x,y,param[1:2]),(1:N))))
  XB<-X%*%param[5:(5+2*N)]
  bp<-cbind(x=XB*sin(a2c(x,y,param[1:2]))+param[1],y=XB*cos(a2c(x,y,param[1:2]))+param[2])
  plot<-ggplot(cbind.data.frame(x,y,value),aes(x,y))+geom_tile(aes(fill=value))+geom_point(data=as.data.frame(bp),aes(V1,V2),size=.5)+scale_fill_distiller(palette = "Spectral", direction = -1)
  
  if(credible.band==T){
    n<-dim(result[[1]])[1]
    f<-array(dim=c(n,200))
    mean_center<-colMeans(result[[1]][,1:2])
    for(i in 1:n){
      center<-result[[1]][i,1:2]
      X=cbind2(cbind2(1,sin(outer(a2c(x,y,center),(1:N)))),cos(outer(a2c(x,y,center),(1:N))))
      f_est<-X%*%result[[1]][i,5:(5+2*N)]
      bp<-cbind.data.frame(x=f_est*sin(a2c(x,y,center))+center[1],y=f_est*cos(a2c(x,y,center))+center[2])
      theta_star<-a2c(bp$x,bp$y,mean_center)
      f_star<-d2c(bp$x,bp$y,mean_center)
      X_star<-cbind2(cbind2(1,sin(outer(theta_star,(1:N)))),cos(outer(theta_star,(1:N))))
      beta_star<-solve(t(X_star)%*%X_star)%*%t(X_star)%*%f_star
      X<-cbind2(cbind2(1,sin(outer(seq(0,2*pi,length.out = 200),(1:N)))),cos(outer(seq(0,2*pi,length.out = 200),(1:N))))
      f[i,]<-X%*%beta_star
    }

    fhat<-colMeans(f)
    angle=seq(0,2*pi,length.out = 200)
    fsd<-colVars(f,std=T)
    dev<-array(dim=c(n,200))
    for(i in 1:n){
      dev[i,]<-abs(f[i,]-fhat)/fsd
    }
    dev1<-rowMaxs(dev,value=T)
    l0<-quantile(dev1,probs=c(.95))
    lower<-fhat-l0*fsd
    upper<-fhat+l0*fsd
    data<-cbind.data.frame(f=c(fhat,lower,upper),angle=seq(0,2*pi,length.out = 200),type=c(rep("estimate",200),rep("lower",200),rep("upper",200)))
    
    polar_band<-ggplot(cbind.data.frame(angle,fhat,lower,upper),aes(angle,fhat))+geom_line(color="black")+geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.4,fill="black")+coord_polar()+ylim(0,1)+theme_minimal()+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),panel.grid = element_blank(),panel.background = element_blank(),panel.grid.major.x=element_blank(), plot.margin = unit(c(-2, 0, -2, -.1), "cm"))
    
    
    point_plot<-ggplot(cbind.data.frame(theta=a2c(x,y,mean_center),dist=d2c(x,y,mean_center),value),aes(theta,dist))+geom_point(aes(color=value),shape=15,size=3)+coord_polar()+scale_color_distiller(palette = "Spectral", direction = -1)+theme_minimal()+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),panel.grid = element_blank(),panel.background = element_blank(),panel.grid.major.x=element_blank(),legend.position = "none", plot.margin = unit(c(-2, 0, -2, -.1), "cm"))+ylim(0,1)+xlim(0,2*pi)
    
    library(cowplot)
    
    plot<-ggdraw() +draw_plot(point_plot) +draw_plot(polar_band,scale=1)
    
  }
  
  print(plot)
}



