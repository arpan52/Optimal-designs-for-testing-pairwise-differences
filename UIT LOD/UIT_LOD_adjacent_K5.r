# for pairs (1,2), (2,3),(3,4),(4,5)

#install.packages("nloptr")
#install.packages("mvtnorm"

library('nloptr')
library(mvtnorm)

UIT_Complete_CDF_K5_simple_power <- function(x,mu,del,n,q){
  
  # Calculating mean and covariance matrix of Z_ij's

  de12 <- ((mu[1]-mu[2])*n)/sigma
  de23 <- ((mu[2]-mu[3])*n)/sigma
  de34 <- ((mu[3]-mu[4])*n)/sigma
  de45 <- ((mu[4]-mu[5])*n)/sigma
  

  a12 <- x[1]*x[3]
  b12 <- (x[1]+x[2])*(x[2]+x[3])
  rho12 <- -sqrt(a12/b12)
  
  rho13 <- 0
  
  rho14<- 0
  
  a23 <- x[2]*x[4]
  b23 <- (x[2]+x[3])*(x[3]+x[4])
  rho23 <- -sqrt(a23/b23)
  rho24 <- 0
  
  a34 <- x[3]*x[5]
  b34 <- (x[3]+x[4])*(x[4]+x[5])
  rho34 <- -sqrt(a34/b34)
  
  
  mu1 <- (sqrt((x[1]*x[2])/(x[1] + x[2])))*de12
  mu2 <- (sqrt((x[2]*x[3])/(x[2] + x[3])))*de23
  mu3 <- (sqrt((x[3]*x[4])/(x[3] + x[4])))*de34
  mu4 <- (sqrt((x[4]*x[5])/(x[4] + x[5])))*de45
  
  m <- c(mu1,mu2,mu3,mu4)
  S <-  rbind(c(1,rho12,rho13,rho14),c(rho12,1,rho23,rho24),c(rho13,rho23,1,rho34),c(rho14,rho24,rho34,1))
  
  lw1 <- c(-q,-q,-q,-q)
  up1 <- c(q,q,q)
  
  
  return(-(1-pmvnorm(mean=m, sigma=S, lower=lw1, upper=up1)))  #Power
  
  
}


# calculating the optimal design. 

mu <- c(0,0,0,0.5,0)           # mu vector
sigma <- 1
N<- 180                       #Total number of subjects
n<- sqrt(N)                   
al <- 0.05/8                  # significance level
q <- qnorm(1-(al))            # critical value
del<-0.6                      # delta


# equality Constraints for x
const1 <- function(x,mu,del,n,q) {
  return( (x[1]+x[2]+x[3]+x[4]+x[5]-1))
}


# Lower and upper bounds
lb <- c(0,0,0,0,0)
ub <- c(1,1,1,1,1)

#initial values
x0 <- c(0.25,0.1,0.2,0.2,0.1)

# Set optimization options.
local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )
opts <- list( "algorithm"= "NLOPT_GN_ISRES",
              "xtol_rel"= 1.0e-15,
              "maxeval"= 20000,
              "local_opts" = local_opts,
              "print_level" = 0 )


res <- nloptr ( x0 = x0,
                eval_f = UIT_Complete_CDF_K5_simple_power,
                mu=mu,
                del=del,
                n=n,
                q=q,
                lb = lb,
                ub = ub,
                eval_g_eq = const1,
                opts = opts
)
print(res)

