# BFSP 
## Bayesian Functional Spatial Partitioning 

Simulation data is included as "data.simulations.rds". The coordinates are saved as the first 2 columns of the matrix. The remaining columns contain the 600 simualted Gaussian process data with Matern covariance. 

The file "helpers.R" contains functions, libraries, and the 2 main functions bfsp() and plot_bfsp(). 

The file "simulation.code.R" generates and plots the boundary estimates for simulated data. 

### bfsp()
Estimates boundary via MCMC. 

#### Arguments

loc: n by 2 matrix of spatial coordinates 

value: n by 1 vector of values measured at each spatial coorinate 

iterations: total number of MCMC iterations 

burn: number of MCMC iterations to discard 

N: number of basis functions 

inits: list of initial values for outer and target regions respectively. First list item contains $(\mu_0,\tau^2,\sigma^2_0)$ and the second list item contains $(x_c,y_c,\mu_1,\sigma^2_1,b_0,b_1,b_2)$

xlim: Vector of lower and upper limit on x. The estimated boundary may not exceed these values of x. 

ylim: Vector of lower and upper limit on y. The estimated boundary may not exceed these values of y.  

#### Value

list containing MCMC samples of parameters of the target region, $(x_c,y_c,\mu_1,\sigma^2_1,b_0,b_1,b_2)$, parameters of the outer region, $(\mu_0,\tau^2,\sigma^2_0)$, and elapsed time. 

### plot_bfsp()
Plots the estimated boundary from bfsp().

#### Arguments

result: output from bfsp()

loc: n by 2 matrix of spatial coordinates 

value: n by 1 vector of values measured at each spatial coorinate 

credible.band: logical. If true, a $95\%$ credible band is plotted along with boundary. 

