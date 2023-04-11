# for pairs (1,2),(1,3)

#install.packages("nloptr")
#install.packages("mvtnorm")


library(nloptr)
library(mvtnorm)

normal_probability_tree_k3 <- function(x,mu,n,sigma,q){
  
  
  # Calculating mean and covariance matrix of Z_ij's
  
  de12 <- ((mu[1]-mu[2])*n)/sigma
  de13 <- ((mu[1]-mu[3])*n)/sigma
  
  y<- (1-x)/2
  a1 <- y*y
  
  b1 <- (x+y)*(x+y)
  
  rho <- sqrt(a1/b1)
  
  mu1 <- (sqrt((x*y)/(x + y)))*de12
  
  mu2 <- (sqrt((x*y)/(x + y)))*de13
  
  m <- c(mu1,mu2)
  S <-  rbind(c(1,rho),c(rho,1))
  
  lw1 <- c(-q,-q)
  up1 <- c(q,q)
  
  return(-(1-pmvnorm(mean=m, sigma=S, lower=lw1, upper=up1)))   
  
}
normal_probability_tree_K3_AVG <- function(x,del,n,sigma,q){

    #LFC's
  
  mu1<- c(-del/2,del/2,-del/2)
  mu2<- c(-del/2,-del/2,del/2)
  return((1/2)*normal_probability_tree_k3(x,mu1,n,sigma,q)+(1/2)*normal_probability_tree_k3(x,mu2,n,sigma,q))  #Expected Power
}

# calculating the optimal design.

N<- 60                       #Total number of subjects
n<- sqrt(N)                   
al <- 0.05/4                  # significance level
q <- qnorm(1-(al))            # critical value
del<-0.5                      # delta
sigma<-1


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
                eval_f = normal_probability_tree_K3_AVG,
                del=del,
                n=n,
                sigma=sigma,
                q=q,
                lb = lb,
                ub = ub,
                opts = opts
)
print(res)
