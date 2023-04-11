# for tree pairs (1,2), (1,3), (1,4),(1,5)
#install.packages("nloptr")
#install.packages("mvtnorm")

library('nloptr')
library(mvtnorm)

UIT_Complete_CDF_K5_tree_design_given <- function(x,mu,del,n,q){
  
  # Calculating mean and covariance matrix of Z_ij's
  
  de12 <- ((mu[1]-mu[2])*n)/sigma
  de13 <- ((mu[1]-mu[3])*n)/sigma
  de14 <- ((mu[1]-mu[4])*n)/sigma
  de15 <- ((mu[1]-mu[5])*n)/sigma
  
  a1 <- x[2]*x[3]
  b1 <- (x[1]+x[2])*(x[1]+x[3])
  rho1 <- sqrt(a1/b1)
  
  a2 <- x[2]*x[4]
  b2 <- (x[1]+x[2])*(x[1]+x[4])
  rho2 <- sqrt(a2/b2)
  
  a3 <- x[2]*x[5]
  b3 <- (x[1]+x[2])*(x[1]+x[5])
  rho3 <- sqrt(a3/b3)
  
  a4 <- x[3]*x[4]
  b4 <- (x[1]+x[3])*(x[1]+x[4])
  rho4 <- sqrt(a4/b4)
  
  a5 <- x[3]*x[5]
  b5 <- (x[1]+x[3])*(x[1]+x[5])
  rho5 <- sqrt(a5/b5)
  
  a6 <- x[4]*x[5]
  b6 <- (x[1]+x[4])*(x[1]+x[5])
  rho6 <- sqrt(a6/b6)
  
  mu1 <- (sqrt((x[1]*x[2])/(x[1] + x[2])))*de12
  mu2 <- (sqrt((x[1]*x[3])/(x[1] + x[3])))*de13
  mu3 <- (sqrt((x[1]*x[4])/(x[1] + x[4])))*de14
  mu4 <- (sqrt((x[1]*x[5])/(x[1] + x[5])))*de15
  
  m <- c(mu1,mu2,mu3,mu4)
  S <-  rbind(c(1,rho1,rho2,rho3),c(rho1,1,rho4,rho5),c(rho2,rho4,1,rho6),c(rho3,rho5,rho6,1))  #Power
  
  lw1 <- c(-q,-q,-q,-q)
  up1 <- c(q,q,q,q)
  
  return(-(1-pmvnorm(mean=m, sigma=S, lower=lw1, upper=up1)))   
  
  
}

# calculating the optimal design. 


mu <- c(0,0,0,0.5,0)      # mu vector
sigma <- 1
N<- 180                       #Total number of subjects
n<- sqrt(N)                   
al <- 0.05/8                 # significance level
q <- qnorm(1-(al))            # critical value
del<-0.6   

# equality Constraints for x
const1 <- function(x,mu,del,n,q) {
  return( rbind(x[1]+x[2]+x[3]+x[4]+x[5]-1))
}

# Lower and upper bounds
lb <- c(0,0,0,0,0)
ub <- c(1,1,1,1,1)

#initial values
x0 <- c(0.3,0.1,0.2,0.2,0.1)
#x0 <- c(1/2,1/4,1/4)
# Optimization for tree k=4

# Set optimization options.
local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )
opts <- list( "algorithm"= "NLOPT_GN_ISRES",
              "xtol_rel"= 1.0e-15,
              "maxeval"= 20000,
              "local_opts" = local_opts,
              "print_level" = 0 )


res <- nloptr ( x0 = x0,
                eval_f = UIT_Complete_CDF_K5_tree_design_given,
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

