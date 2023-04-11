# for tree pairs (1,2), (1,3)
#install.packages("nloptr")
#install.packages("mvtnorm")

library('nloptr')
library(mvtnorm)

UIT_Complete_CDF_K3_tree_power <- function(x,mu,del,n,q){
  
  # Calculating mean and covariance matrix of Z_ij's
  

  de12 <- ((mu[1]-mu[2])*n)/sigma
  de13 <- ((mu[1]-mu[3])*n)/sigma

 
  a1 <- x[2]*x[3]
  b1 <- (x[1]+x[2])*(x[1]+x[3])
  rho <- sqrt(a1/b1)
  
  mu1 <- (sqrt((x[1]*x[2])/(x[1] + x[2])))*de12
  mu2 <- (sqrt((x[1]*x[3])/(x[1] + x[3])))*de13
  
  m <- c(mu1,mu2)
  S <-  rbind(c(1,rho),c(rho,1))
  
  lw1 <- c(-q,-q)
  up1 <- c(q,q)
  
  return(-(1-pmvnorm(mean=m, sigma=S, lower=lw1, upper=up1)))   #POWER

  
}
# calculating the optimal design. 


mu <- c(0,0,0.5)        # mu vector
sigma <- 1
N<- 180                       #Total number of subjects
n<- sqrt(N)                   
al <- 0.05/4                 # significance level
q <- qnorm(1-(al))            # critical value
del<-0.6                      # delta

# equality Constraints for x
const1 <- function(x,mu,del,n,q) {
  return( (x[1]+x[2]+x[3]-1))
}

# Lower and upper bounds
lb <- c(0,0,0)
ub <- c(1,1,1)

#initial values
x0 <- c(0.33,0.33,0.33)

# Set optimization options.
local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )
opts <- list( "algorithm"= "NLOPT_GN_ISRES",
              "xtol_rel"= 1.0e-15,
              "maxeval"= 20000,
              "local_opts" = local_opts,
              "print_level" = 0 )


res <- nloptr ( x0 = x0,
                eval_f = UIT_Complete_CDF_K3_tree_power,
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

