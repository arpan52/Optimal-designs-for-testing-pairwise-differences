# for pairs (1,2),(1,3),(1,4),(1,5)

#install.packages("nloptr")
#install.packages("mvtnorm"
#install.packages("tictoc")

library(tictoc)
library('nloptr')
library(mvtnorm)



UIT_tree_CDF_K5_avg_LFC_design <- function(x,mu,n,sigma,q){
  
  # Calculating mean and covariance matrix of Z_ij's
  
  y<- (1-x)/4

  de12 <- ((mu[1]-mu[2])*n)/sigma
  de13 <- ((mu[1]-mu[3])*n)/sigma
  de14 <- ((mu[1]-mu[4])*n)/sigma
  de15 <- ((mu[1]-mu[5])*n)/sigma
  
  a <- y*y
  b <- (x+y)*(x+y)
  rho <- sqrt(a/b)

  mu1 <- (sqrt((x*y)/(x + y)))*de12
  mu2 <- (sqrt((x*y)/(x + y)))*de13
  mu3 <- (sqrt((x*y)/(x + y)))*de14
  mu4 <- (sqrt((x*y)/(x + y)))*de15

  m <- c(mu1,mu2,mu3,mu4)
  S <-  rbind(c(1,rho,rho,rho),c(rho,1,rho,rho),c(rho,rho,1,rho),c(rho,rho,rho,1))
  
  lw1 <- c(-q,-q,-q,-q)
  up1 <- c(q,q,q,q)
  
  return(-(1-pmvnorm(mean=m, sigma=S, lower=lw1, upper=up1)))   
  
  
}
normal_probability_tree_K5_AVG <- function(x,del,n,sigma,q){
 
  #LFC's
  
  mu11<- c(-del/2,del/2,-del/2,-del/2,-del/2)
  mu21<- c(-del/2,-del/2,del/2,-del/2,-del/2)
  mu31<- c(-del/2,-del/2,-del/2,del/2,-del/2)
  mu41<- c(-del/2,-del/2,-del/2,-del/2,del/2) 
  
  # Expected Power
  
  return((1/4)*UIT_tree_CDF_K5_avg_LFC_design(x,mu11,n,sigma,q)+(1/4)*UIT_tree_CDF_K5_avg_LFC_design(x,mu21,n,sigma,q)+(1/4)*UIT_tree_CDF_K5_avg_LFC_design(x,mu31,n,sigma,q)+(1/4)*UIT_tree_CDF_K5_avg_LFC_design(x,mu41,n,sigma,q))
  
}

# calculating the optimal design. 



N<- 250                       #Total number of subjects
n<- sqrt(N)                   
al <- 0.05/8                  # significance level
q <- qnorm(1-(al))            # critical value
del<-0.4                      # delta


# Lower and upper bounds
lb <- c(0)
ub <- c(1)

#initial values
x0 <- c(0.35)

# Set optimization options.
local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )
opts <- list( "algorithm"= "NLOPT_GN_ISRES",
              "xtol_rel"= 1.0e-15,
              "maxeval"= 20000,
              "local_opts" = local_opts,
              "print_level" = 0 )


res <- nloptr ( x0 = x0,
                eval_f = normal_probability_tree_K5_AVG,
                del=del,
                n=n,
                sigma=sigma,
                q=q,
                lb = lb,
                ub = ub,
                opts = opts
)
print(res)
toc()