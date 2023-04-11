# for pairs (1,2),(1,3),(1,4)

#install.packages("nloptr")
#install.packages("mvtnorm")
library(mvtnorm)

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

# calculating the optimal design.

N<-180                       #Total number of subjects
n<- sqrt(N)                   
al <- 0.05                  # significance level
q <- qnorm(1-(al))            # critical value
del<-0.6                      # delta
sigma<-1

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
print(res)
