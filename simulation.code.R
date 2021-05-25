
# Source helper functions 
source("helpers.R")

# read in data 
data<-readRDS("data.simulations.rds")

# set job number 
job=1
loc<-data[,1:2]
value<-data[,job+2]

# Set parameters for mcmc 
iterations<-50000 # number of iterations 
burn<-10000 #number of samples to throw out as burn in 
N<-5 # number of basis functions 
inits<-list(c(-1,1,1),c(.5,.5,1,1,.4,rep(0,2*N))) #set initial values 
xlim<-c(-.1,1.1) #limits on x-axis 
ylim<-c(-.1,1.1) #limits on y-axis

#run mcmc 
result<- bfsp(loc,value,iterations,burn,N,inits,xlim,ylim)

#plot results 
plot_bfsp(result,loc,value,credible.band =F)
