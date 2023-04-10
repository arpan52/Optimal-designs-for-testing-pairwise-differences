###########################################################################
##### Simulation study to check the performance of the Max--min design for
#####  IUT Tree graph versus the Balanced design versus Dunnett's design for K=4
###########################################################################

# load the following libraries

library('nloptr')
library(mvtnorm)
rm(list = ls())


IUT_tree_K4 <- function(x,mu,n,sigma,q){
  # Calculating mean and covariance matrix of Z_ij's
  
  y<- (1-x)/3
  
  mu12 <- (n*(mu[1]-mu[2]))/(sqrt((1/x)+(1/y)))
  mu13 <- (n*(mu[1]-mu[3]))/(sqrt((1/x)+(1/y)))
  mu14 <- (n*(mu[1]-mu[4]))/(sqrt((1/x)+(1/y)))
  
  
  p12<- pnorm(q,mu12)-pnorm(-q,mu12)
  p13<-pnorm(q,mu13)-pnorm(-q,mu13)
  p14<-pnorm(q,mu14)-pnorm(-q,mu14)
  
  
  
  lw1 <- c(-q,-q)
  up1 <- c(q,q)
  lw2 <- c(-q,-q,-q)
  up2 <- c(q,q,q)
  
  rho <- sqrt((y*y)/((x+y)*(x+y)))
  s<- rbind(c(1,rho),c(rho,1))
  mu1213<- c(mu12,mu13)
  p1213<- pmvnorm(mean=mu1213, sigma=s, lower=lw1, upper=up1)
  mu1214<- c(mu12,mu14)
  p1214<- pmvnorm(mean=mu1214, sigma=s, lower=lw1, upper=up1)
  mu1314<- c(mu13,mu14)
  p1314<- pmvnorm(mean=mu1314, sigma=s, lower=lw1, upper=up1)
  mu121314<- c(mu12,mu13,mu14)
  s121314<- rbind(c(1,rho,rho),c(rho,1,rho),c(rho,rho,1))
  p121314<- pmvnorm(mean=mu121314, sigma=s121314, lower=lw2, upper=up2)
  
  
  
  p<- -(1-p12-p13-p14+p1213+p1214+p1314-p121314)
  return(p)
}
IUT_tree_K4_AVG <- function(x,del,n,sigma,q){
  
  
  # LFC's
  mu1<- c(0,-del,del,del)
  mu2<- c(0,del,-del,del)
  mu3<- c(0,del,del,-del)
  
  
  return((1/3)*IUT_tree_K4(x,mu1,n,sigma,q)+(1/3)*IUT_tree_K4(x,mu2,n,sigma,q)+(1/3)*IUT_tree_K4(x,mu3,n,sigma,q))  #Expected Power
}

design <- res$solution

# calculating the optimal design.
N<-24                # Total sample size
n <- sqrt(N)
sigma <- 1.16        # standard deviation
al <- 0.05          # significance level alpha
q <- qnorm(1-al)      # critical value
del<-1.13                     # delta

# Lower and upper bounds
lb <- c(0)
ub <- c(1)

#initial values
x0 <- 0.3


# Set optimization options.
local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )
opts <- list( "algorithm"= "NLOPT_GN_ISRES",
              "xtol_rel"= 1.0e-15,
              "maxeval"= 20000,
              "local_opts" = local_opts,
              "print_level" = 0 )


res <- nloptr ( x0 = x0,
                eval_f = IUT_tree_K4_AVG,
                del=del,
                n=n,
                sigma=sigma,
                q=q,
                lb = lb,
                ub = ub,
                opts = opts
)

design <- res$solution
# Power function for UIT Tree graph K=3


# values of the parameters associated


m<-1e6             # number of simulations
K<-4

eta<-design            # max--min design for path IUT K=3



x_mm<-c(eta,(1-eta)/3,(1-eta)/3,(1-eta)/3)
x_d1<- 1/(1+sqrt(K-1))

x_d <- c(x_d1,(1-x_d1)/3,(1-x_d1)/3,(1-x_d1)/3)


true.mean<-c(6.70,14.12,9.02,7.83)   # true mean from real data
psd <- 1.16                  # pooled standard deviation of the real data



# Simulations for max--min design



n_mm<- rep(0,K)
for (i in 1:K) {
  n_mm[i]<- ceiling(x_mm[i]*N)  # exact max--min design
}
n_mm<-c(19,9,10,10)
mean_smm<-0
pwr_mm<- rep(0,m)

n_d <- rep(0,K)
for (i in 1:K) {
  n_d[i]<- ceiling(x_d[i]*N)  # exact Dunnett's design
}
mean_d<-0
pwr_d<- rep(0,m)

n_mm<- rep(0,K)
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
  
  pwr_mm[j] <- IUT_tree_K4(eta,mean_smm,n,sigma,q)    # power vector for the simulated data for the min--max design
  pwr_b[j] <- IUT_tree_K4(1/K,mean_b,n,sigma,q)       # power vector for the simulated data for the balanced design
  pwr_d[j] <- IUT_tree_K4(x_d1,mean_d,n,sigma,q)      # power vector for the simulated data for the Dunnett design
  
  
}



print('Median of power for the max-min design')

print(median(-pwr_mm))

print('Median of power for the balanced design')

print(median(-pwr_b))

print('Median of power for the Dunnett design')

print(median(-pwr_d))