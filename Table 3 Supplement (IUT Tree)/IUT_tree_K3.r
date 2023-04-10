# for pairs (1,2),(1,3)

#install.packages("nloptr")
#install.packages("mvtnorm")
library(nloptr)
library(mvtnorm)

IUT_tree_K3 <- function(x,mu,n,sigmq,q){
  
  # Calculating mean and covariance matrix of Z_ij's

  y<- (1-x)/2

  mu12 <- (n*(mu[1]-mu[2]))/(sigma*sqrt((1/x)+(1/y)))
  mu13 <- (n*(mu[1]-mu[3]))/(sigma*sqrt((1/x)+(1/y)))
 
  p12<- pnorm(q,mu12)-pnorm(-q,mu12)
  p13<-pnorm(q,mu13)-pnorm(-q,mu13)


  lw1 <- c(-q,-q)
  up1 <- c(q,q)

  rho <- sqrt((y*y)/((x+y)*(x+y)))
  s<- rbind(c(1,rho),c(rho,1))
  mu1213<- c(mu12,mu13)
  p1213<- pmvnorm(mean=mu1213, sigma=s, lower=lw1, upper=up1)

  p<- -(1-p12-p13+p1213)
  return(p)
}
IUT_tree_K3_AVG <- function(x,del,n,sigma,q){
  
  # LFC's
  mu1<- c(0,-del,del)
  mu2<- c(0,del,-del)

   return((1/2)*IUT_tree_K3(x,mu1,n,sigmq,q)+(1/2)*IUT_tree_K3(x,mu2,n,sigmq,q))   # Expeceted power
}

# calculating the optimal design.

N<- 60                       #Total number of subjects
n<- sqrt(N)                   
al <- 0.05                  # significance level
q <- qnorm(1-(al))            # critical value
del<-0.5                      # delta
sigma <- 1


# Lower and upper bounds
lb <- 0.25
ub <- 0.5

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
                eval_f = IUT_tree_K3_AVG,
                del=del,
                n=n,
                sigma=sigma,
                q=q,
                lb = lb,
                ub = ub,
                opts = opts
)
print(res)
