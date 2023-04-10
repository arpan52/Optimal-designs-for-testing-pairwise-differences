# for pairs (1,2),(1,3),(1,4),(1,5)

#install.packages("nloptr")
#install.packages("mvtnorm"
#install.packages("tictoc")

library(tictoc)
tic()
library(nloptr)
library(mvtnorm)

IUT_tree_K5 <- function(x,mu,n,sigma,q){
  # Calculating mean and covariance matrix of Z_ij's
  
  y<- (1-x)/4

  mu12 <- (n*(mu[1]-mu[2]))/(sqrt((1/x)+(1/y)))
  mu13 <- (n*(mu[1]-mu[3]))/(sqrt((1/x)+(1/y)))
  mu14 <- (n*(mu[1]-mu[4]))/(sqrt((1/x)+(1/y)))
  mu15 <- (n*(mu[1]-mu[5]))/(sqrt((1/x)+(1/y)))
  
 
  
  p12<- pnorm(q,mu12)-pnorm(-q,mu12)
  p13<-pnorm(q,mu13)-pnorm(-q,mu13)
  p14<-pnorm(q,mu14)-pnorm(-q,mu14)
  p15<-pnorm(q,mu15)-pnorm(-q,mu15)
  
  
  
  lw1 <- c(-q,-q)
  up1 <- c(q,q)
  lw2 <- c(-q,-q,-q)
  up2 <- c(q,q,q)
  lw3 <- c(-q,-q,-q,-q)
  up3 <- c(q,q,q,q) 
  
  rho <- sqrt((y*y)/((x+y)*(x+y)))
  s<- rbind(c(1,rho),c(rho,1))
  mu1213<- c(mu12,mu13)
  p1213<- pmvnorm(mean=mu1213, sigma=s, lower=lw1, upper=up1)
  mu1214<- c(mu12,mu14)
  p1214<- pmvnorm(mean=mu1214, sigma=s, lower=lw1, upper=up1)
  mu1215<- c(mu12,mu15)
  p1215<- pmvnorm(mean=mu1215, sigma=s, lower=lw1, upper=up1)
  mu1314<- c(mu13,mu14)
  p1314<- pmvnorm(mean=mu1314, sigma=s, lower=lw1, upper=up1)
  mu1315<- c(mu13,mu15)
  p1315<- pmvnorm(mean=mu1315, sigma=s, lower=lw1, upper=up1)
  mu1415<- c(mu14,mu15)
  p1415<- pmvnorm(mean=mu1415, sigma=s, lower=lw1, upper=up1)
  mu121314<- c(mu12,mu13,mu14)
  s1<- rbind(c(1,rho,rho),c(rho,1,rho),c(rho,rho,1))
  p121314<- pmvnorm(mean=mu121314, sigma=s1, lower=lw2, upper=up2)
  mu121315<- c(mu12,mu13,mu15)
  p121315<- pmvnorm(mean=mu121315, sigma=s1, lower=lw2, upper=up2)
  mu131415<- c(mu13,mu14,mu15)
  p131415<- pmvnorm(mean=mu131415, sigma=s1, lower=lw2, upper=up2)
  mu12131415<- c(mu12,mu13,mu14,mu15)
  s2<- rbind(c(1,rho,rho,rho),c(rho,1,rho,rho),c(rho,rho,1,rho),c(rho,rho,rho,1))
  p12131415<- pmvnorm(mean=mu12131415, sigma=s2, lower=lw3, upper=up3)
  p<- -(1-p12-p13-p14-p15+p1213+p1214+p1215+p1314+p1315+p1415-p121314-p121315-p131415+p12131415)
  return(p)
}
IUT_tree_K5_AVG <- function(x,del,n,sigma,q){

  # LFC's
  mu1<- c(0,-del,-del,del,del)
  mu2<- c(0,-del,del,-del,del)
  mu3<- c(0,-del,del,del,-del)
  mu4<- c(0,del,-del,-del,del)
  mu5<- c(0,del,-del,del,-del)
  mu6<- c(0,del,del,-del,-del)
  
  return((1/6)*IUT_tree_K5(x,mu1,n,sigma,q)+(1/6)*IUT_tree_K5(x,mu2,n,sigma,q)+(1/6)*IUT_tree_K5(x,mu3,n,sigma,q)+(1/6)*IUT_tree_K5(x,mu4,n,sigma,q)+(1/6)*IUT_tree_K5(x,mu5,n,sigma,q)+(1/6)*IUT_tree_K5(x,mu6,n,sigma,q))
}

# calculating the optimal design. 



N<- 250                       #Total number of subjects
n<- sqrt(N)                   
al <- 0.05                  # significance level
q <- qnorm(1-(al))            # critical value
del<-0.55                     # delta
sigma<-1

# Lower and upper bounds
lb <- 0.1
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
                eval_f = IUT_tree_K5_AVG,
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
