
###########################################################################
##### Simulation study to check the performance of the Max--min design for
#####  UIT Tree graph versus the Balanced design versus Dunnett's design for K=4
###########################################################################

# load the following libraries

library('nloptr')
library(mvtnorm)
rm(list = ls())

# calculating the max--min design

UIT_tree_CDF_K4_avg_LFC_design <- function(x,mu,n,sigma,q){
  
  # Calculating mean and covariance matrix of Z_ij's
  
  y<- (1-x)/3
  
  de12 <- ((mu[1]-mu[2])*n)/sigma
  de13 <- ((mu[1]-mu[3])*n)/sigma
  de14 <- ((mu[1]-mu[4])*n)/sigma
  
  
  a <- y*y
  b <- (x+y)*(x+y)
  rho <- sqrt(a/b)
  
  mu1 <- (sqrt((x*y)/(x + y)))*de12
  mu2 <- (sqrt((x*y)/(x + y)))*de13
  mu3 <- (sqrt((x*y)/(x + y)))*de14
  
  m <- c(mu1,mu2,mu3)
  S <-  rbind(c(1,rho,rho),c(rho,1,rho),c(rho,rho,1))
  
  lw1 <- c(-q,-q,-q)
  up1 <- c(q,q,q)
  
  return(-(1-pmvnorm(mean=m, sigma=S, lower=lw1, upper=up1)))   
  
  
  
}
normal_probability_tree_K4_AVG <- function(x,del,n,sigma,q){
  
  
  #LFC's
  
  mu11<- c(-del/2,del/2,-del/2,-del/2)
  mu21<- c(-del/2,-del/2,del/2,-del/2)
  mu31<- c(-del/2,-del/2,-del/2,del/2)
  
  # Expected Power
  return((1/3)*UIT_tree_CDF_K4_avg_LFC_design(x,mu11,n,sigma,q)+(1/3)*UIT_tree_CDF_K4_avg_LFC_design(x,mu21,n,sigma,q)+(1/3)*UIT_tree_CDF_K4_avg_LFC_design(x,mu31,n,sigma,q))
  
}

# calculating the optimal design.

N<-48                # Total sample size
n <- sqrt(N)
sigma <- 8.71        # standard deviation
al <- 0.05/6          # significance level alpha
q <- qnorm(1-al)      # critical value
del<-5.8                     # delta



# Lower and upper bounds
lb <- c(0)
ub <- c(1)

#initial values
x0 <- c(0.3)

# Set optimization options.
local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )
opts <- list( "algorithm"= "NLOPT_GN_ISRES",
              "xtol_rel"= 1.0e-15,
              "maxeval"= 20000,
              "local_opts" = local_opts,
              "print_level" = 0 )


res <- nloptr ( x0 = x0,
                eval_f = normal_probability_tree_K4_AVG,
                del=del,
                n=n,
                sigma=sigma,
                q=q,
                lb = lb,
                ub = ub,
                opts = opts
)
design <- res$solution


# values of the parameters associated


m<-1e6             # number of simulations
K<-4

eta<-design             # max--min design for path IUT K=3

x_mm<-c(eta,(1-eta)/3,(1-eta)/3,(1-eta)/3)
x_d1<- 1/(1+sqrt(K-1))

x_d <- c(x_d1,(1-x_d1)/3,(1-x_d1)/3,(1-x_d1)/3)


true.mean<-c(97.6,101.6,102.2,103.4)   # true mean from real data
psd <- 8.71                   # pooled standard deviation of the real data



n_mm<- rep(0,K)
for (i in 1:K) {
  n_mm[i]<- ceiling(x_mm[i]*N)  # exact max--min design
}
mean_smm<-0
pwr_mm<- rep(0,m)

n_d <- rep(0,K)
for (i in 1:K) {
  n_d[i]<- ceiling(x_d[i]*N)  # exact Dunnett's design
}
mean_d<-0
pwr_d<- rep(0,m)

n_b<- rep(0,K)
for (i in 1:K) {
  n_b[i] <- ceiling( N/K)   # exact balanced design
}
mean_b<-0
pwr_b<- rep(0,m)

# generate from normal distribution

for (j in 1:m){
  
  for (i in 1:K) {
    
    # Simulations for Min--Max design
    sample_mm <- rep(0,n_mm[i])
    sample_mm <- rnorm(n_mm[i], mean=true.mean[i], psd)
    mean_smm[i] <- mean(sample_mm)
    
    
    # Simulations for Balanced design
    sample_b <- rep(0,n_b[i])
    sample_b <- rnorm(n_b[i], mean=true.mean[i], psd)
    mean_b[i] <- mean(sample_b)
    
    # Simulations for Dunnett design
    sample_d <- rep(0,n_d[i])
    sample_d <- rnorm(n_d[i], mean=true.mean[i], psd)
    mean_d[i] <- mean(sample_d)
  }
  
  pwr_mm[j] <- UIT_tree_CDF_K4_avg_LFC_design(eta,mean_smm,n,sigma,q)   # power vector for the simulated data for the min--max design
  pwr_b[j] <- UIT_tree_CDF_K4_avg_LFC_design(1/K,mean_b,n,sigma,q)      # power vector for the simulated data for the balanced design
  pwr_d[j] <- UIT_tree_CDF_K4_avg_LFC_design(x_d1,mean_d,n,sigma,q)     # power vector for the simulated data for the Dunnett's design
  
  
}



print('Median of power for the max-min design')

print(median(-pwr_mm))

print('Median of power for the balanced design')

print(median(-pwr_b))

print('Median of power for the Dunnett design')

print(median(-pwr_d))

