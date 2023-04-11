
# for pairs (1,2), (1,3),(1,4),(2,3),(2,4),(3,4)
#install.packages("nloptr")
#install.packages("mvtnorm")

library('nloptr')
library(mvtnorm)

UIT_Complete_CDF_K4_power_singular <- function(x,mu,del,n,q){
 
  # Calculating mean and covariance matrix of Z_ij's

  de12 <- ((mu[1]-mu[2])*n)/sigma
  de23 <- ((mu[2]-mu[3])*n)/sigma
  de34 <- ((mu[3]-mu[4])*n)/sigma
  de13 <- ((mu[1]-mu[3])*n)/sigma
  de14 <- ((mu[1]-mu[4])*n)/sigma
  de24 <- ((mu[2]-mu[4])*n)/sigma
  
  

  a12 <- x[1]*x[3]
  b12 <- (x[1]+x[2])*(x[2]+x[3])
  rho12 <- -sqrt(a12/b12)
  
  rho13 <- 0
  
  a14 <- x[2]*x[3]
  b14 <- (x[1]+x[2])*(x[1]+x[3])
  rho14 <- sqrt(a14/b14)
  
  a15 <- x[2]*x[4]
  b15 <- (x[1]+x[2])*(x[1]+x[4])
  rho15 <- sqrt(a15/b15)
  
  a16 <- x[1]*x[4]
  b16 <- (x[1]+x[2])*(x[2]+x[4])
  rho16 <- -sqrt(a16/b16)
  
  a23 <- x[2]*x[4]
  b23 <- (x[2]+x[3])*(x[3]+x[4])
  rho23 <- -sqrt(a23/b23)
  
  a24 <- x[1]*x[2]
  b24 <- (x[2]+x[3])*(x[1]+x[3])
  rho24 <- sqrt(a24/b24)
  
  rho25 <- 0
  
  a26 <- x[3]*x[4]
  b26 <- (x[2]+x[3])*(x[2]+x[4])
  rho26 <- sqrt(a26/b26)
  
  a34 <- x[1]*x[4]
  b34 <- (x[1]+x[3])*(x[3]+x[4])
  rho34 <- -sqrt(a34/b34)
  
  a35 <- x[1]*x[3]
  b35 <- (x[1]+x[4])*(x[3]+x[4])
  rho35 <- sqrt(a35/b35)
  
  a36 <- x[2]*x[3]
  b36 <- (x[3]+x[4])*(x[2]+x[4])
  rho36 <- sqrt(a36/b36)
  
  a45 <- x[3]*x[4]
  b45 <- (x[1]+x[3])*(x[1]+x[4])
  rho45 <- sqrt(a45/b45)

  rho46 <- 0
  
  a56 <- x[1]*x[2]
  b56 <- (x[1]+x[4])*(x[2]+x[4])
  rho56 <- sqrt(a56/b56)
 
  
  mu1 <- (sqrt((x[1]*x[2])/(x[1] + x[2])))*de12
  mu2 <- (sqrt((x[2]*x[3])/(x[2] + x[3])))*de23
  mu3 <- (sqrt((x[3]*x[4])/(x[3] + x[4])))*de34
  mu4 <- (sqrt((x[1]*x[3])/(x[1] + x[3])))*de13
  mu5 <- (sqrt((x[1]*x[4])/(x[1] + x[4])))*de14
  mu6 <- (sqrt((x[2]*x[4])/(x[2] + x[4])))*de24
  
  m <- c(mu1,mu2,mu3,mu4,mu5,mu6)
  S <-  rbind(c(1,rho12,rho13,rho14,rho15,rho16),c(rho12,1,rho23,rho24,rho25,rho26),c(rho13,rho23,1,rho34,rho35,rho36),c(rho14,rho24,rho34,1,rho45,rho46),c(rho15,rho25,rho35,rho45,1,rho56),c(rho16,rho26,rho36,rho46,rho56,1))
  
  lw1 <- c(-q,-q,-q,-q,-q,-q)
  up1 <- c(q,q,q,q,q,q)
  
  return(-(1-pmvnorm(mean=m, sigma=S, lower=lw1, upper=up1)))   #Power 
  
  
}

# calculating the optimal design. 

mu <- c(0,-0.25,0.25,0)        # mu vector
sigma <- 1
N<- 180                       #Total number of subjects
n<- sqrt(N)                   
al <- 0.05/12                  # significance level
q <- qnorm(1-(al))            # critical value
del<-0.6                      # delta


# equality Constraints for x
const1 <- function(x,mu,del,n,q) {
  return( (x[1]+x[2]+x[3]+x[4]-1))
}

# Lower and upper bounds
lb <- c(0,0,0,0)
ub <- c(1,1,1,1)

#initial values
x0 <- c(0.1,0.2,0.3,0.2)


# Set optimization options.
local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )
opts <- list( "algorithm"= "NLOPT_GN_ISRES",
              "xtol_rel"= 1.0e-15,
              "maxeval"= 20000,
              "local_opts" = local_opts,
              "print_level" = 0 )


res <- nloptr ( x0 = x0,
                eval_f = UIT_Complete_CDF_K4_power_singular,
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

