###########################################################################
##### Simulation study to check the performance of the Max--min design for
#####  IUT Path graph versus the Balanced design for K=3
###########################################################################

# load the following libraries


library('nloptr')
library(mvtnorm)
rm(list = ls())


# for pairs (1,2),(2,3)

#install.packages("nloptr")
#install.packages("mvtnorm")
library(nloptr)
library(mvtnorm)

IUT_simple_K3 <- function(x,mu,n,sigma,q){
  
  # Calculating mean and covariance matrix of Z_ij 's
  
  t<-(1-2*x)/2
  mu12 <- (n*(mu[1]-mu[2]))/(sigma*sqrt((1/x)+(1/t)))
  mu23 <- (n*(mu[2]-mu[3]))/(sigma*sqrt((1/t)+(1/x)))
  
  
  p12<- pnorm(q,mu12)-pnorm(-q,mu12)
  p23<-pnorm(q,mu23)-pnorm(-q,mu23)
  
  lw1 <- c(-q,-q)
  up1 <- c(q,q)
  rho1223 <- -sqrt((x*x)/((x+t)*(t+x)))
  s1223<- rbind(c(1,rho1223),c(rho1223,1))
  mu1223<- c(mu12,mu23)
  p1223<- pmvnorm(mean=mu1223, sigma=s1223, lower=lw1, upper=up1)      
  
  p<- -(1-p12-p23+p1223)         # Power of the test
  return(p)
}

# calculating the optimal design.

N<-18                # Total sample size
n <- sqrt(N)
sigma <- 1.02        # standard deviation
al <- 0.05          # significance level alpha
q <- qnorm(1-al)      # critical value
del<-2.06                       # delta value
mu<- c(0,del,2*del)              # mu value


# Lower and upper bounds
lb <- 0
ub <- 0.5

#initial values
x0 <- 0.3

local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )
opts <- list( "algorithm"= "NLOPT_GN_ISRES",
              "xtol_rel"= 1.0e-15,
              "maxeval"= 20000,
              "local_opts" = local_opts,
              "print_level" = 0 )


res <- nloptr ( x0 = x0,
                eval_f = IUT_simple_K3,
                mu = mu,
                n = n,
                sigma = sigma,
                q = q,
                lb = lb,
                ub = ub,
                opts = opts
)
design <- res$solution



# values of the parameters associated


m<-1e6              # number of simulations


eta<-design             # max--min design for path IUT K=3
x_mm<-c(eta,1-2*eta,eta)
K<-3

true.mean<-c(6.7,9.02,6.96)   # true mean from real data
psd <- 1.02                   # pooled standard deviation of the real data


n_mm<- rep(0,K)
for (i in 1:K) {
  n_mm[i]<- ceiling(x_mm[i]*N)  # exact max--min design
}

mean_smm<-0
pwr_mm<- rep(0,m)

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
  }
  
  pwr_mm[j] <- IUT_simple_K3(eta,mean_smm,n,sigma,q) # power vector for the simulated data for the min--max design
  pwr_b[j] <- IUT_simple_K3(1/K,mean_b,n,sigma,q)    # power vector for the simulated data for the balanced design
  
}

print('Median of power for the max--min design')

print(median(-pwr_mm))

print('Median of power for the balanced design')

print(median(-pwr_b))


