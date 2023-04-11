# for pairs (1,2),(2,3),(3,4)

#install.packages("nloptr")
#install.packages("mvtnorm")
library(nloptr)
library(mvtnorm)

IUT_simple_K4 <- function(x,del,mu,n,sigma,q){


   # Calculating mean and covariance matrix of Z_ij 's
  
  t<-(1-2*x)/2
  mu12 <- (n*(mu[1]-mu[2]))/(sigma*sqrt((1/x)+(1/t)))
  mu23 <- (n*(mu[2]-mu[3]))/(sigma*sqrt((1/t)+(1/t)))
  mu34 <- (n*(mu[3]-mu[4]))/(sigma*sqrt((1/t)+(1/x)))
  
  p12<- pnorm(q,mu12)-pnorm(-q,mu12)
  p23<-pnorm(q,mu23)-pnorm(-q,mu23)
  p34<- pnorm(q,mu34)-pnorm(-q,mu34)
  
  lw1 <- c(-q,-q)
  up1 <- c(q,q)
  rho1223 <- -sqrt((x*t)/((x+t)*(t+t)))
  s1223<- rbind(c(1,rho1223),c(rho1223,1))
  mu1223<- c(mu12,mu23)
  p1223<- pmvnorm(mean=mu1223, sigma=s1223, lower=lw1, upper=up1)

  rho2334 <- -sqrt((t*x)/((t+t)*(t+x)))
  s2334<- rbind(c(1,rho2334),c(rho2334,1))
  mu2334<- c(mu23,mu34)
  p2334<- pmvnorm(mean=mu2334, sigma=s2334, lower=lw1, upper=up1)

  rho1234 <- 0
  s1234<- rbind(c(1,rho1234),c(rho1234,1))
  mu1234<- c(mu12,mu34)
  p1234<-  pmvnorm(mean=mu1234, sigma=s1234, lower=lw1, upper=up1)

  lw2 <- c(-q,-q,-q)
  up2 <- c(q,q,q)
 s122334<- rbind(c(1,rho1223,rho1234),c(rho1223,1,rho2334),c(rho1234,rho2334,1))
  mu122334<- c(mu12,mu23,mu34)
  p122334<-  pmvnorm(mean=mu122334, sigma=s122334, lower=lw2, upper=up2)
  

  p<- -(1-p12-p23-p34+p1223+p2334+p1234-p122334)        # Power of the test
  return(p)
}

 # calculating the optimal design

del<-0.6                       # delta value
mu<- c(0,del,2*del,3*del)       # mu value
N<- 180                          # Total number of subjects
n<- sqrt(N)
al <- 0.05                      # significance level
q <- qnorm(1-(al))                  # Critical value
sigma<-1
# Lower and upper bounds
lb <- 0
ub <- 0.5

#initial values
x0 <- 0.45


# Set optimization options.
local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )
opts <- list( "algorithm"= "NLOPT_GN_ISRES",
              "xtol_rel"= 1.0e-15,
              "maxeval"= 20000,
              "local_opts" = local_opts,
              "print_level" = 0 )


res <- nloptr ( x0 = x0,
                eval_f = IUT_simple_K4,
                del = del,
                mu = mu,
                n = n,
                sigma=sigma,
                q = q,
                lb = lb,
                ub = ub,
                opts = opts
)
print(res)
x2<- (1-2*res$solution)/2
x2