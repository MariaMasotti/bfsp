
library(Rfast);library(fields);library(FastGP);library(sparseMVN);library(extraDistr);library(Matrix)


# Source helper functions 
source("helpers.R")

# Define parameters for all simulation settings 
theta1<-.05
theta0<-.05
smooths<-c(.5,1.5)
nuggets<-c(.1,.5)
shapes<-c("square","heart","triangle")
grid<-expand.grid(1:50,smooths,nuggets,shapes)

# set job number 
job=1
cat(paste0("Processing Job #", job, "\n"))

# Generate data for job 
smooth<-grid[job,2]
nugget<-grid[job,3]
shape<-grid[job,4]
set.seed(job)
data<-simulated.matern.data(smooth,nugget,shape)
loc<-data[[1]]
value<-data[[2]]

# Set parameters for mcmc 
N<-5 # number of basis functions 
iterations<-10000# number of iterations 
inits<-list(c(-1,1,1),c(.5,.5,1,1,.4,rep(0,2*N))) #initial values 
xlim<-c(-.1,1.1) #limits on x-axis (boundary cannot go outside these limits)
ylim<-c(-.1,1.1) #limits on y-axis
burn<-1000

#run mcmc 
result<- bfsp(loc,value,iterations,burn,N,inits,xlim,ylim)

#plot results 
plot_bfsp(result,loc,value,credible.band = F)
